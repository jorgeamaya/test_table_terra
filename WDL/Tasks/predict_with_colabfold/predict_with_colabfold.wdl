version 1.0

task TestPredictWithColabfold {

    input {
        String bucket_name
        String screen_id
        File fasta_input_group_inventory
    }

    command <<<
        set -euxo pipefail

        # Prepare directory structure

        bucket_name="~{bucket_name}"
        screen_id="~{screen_id}"

        BUCKET_URI="gs://${bucket_name}"
        SCREEN_ROOT="${BUCKET_URI}/screens/${screen_id}_screen"

        # Get group value

        fasta_input_group_inventory="~{fasta_input_group_inventory}"
        
        GROUP_ID="$(awk -F '\t' 'NR==2 {print $4}' "${fasta_input_group_inventory}")"
        if [[ -z "${GROUP_ID}" ]]; then
            GROUP_ID=$(basename "${fasta_input_group_inventory}" .tsv | sed 's/fasta_input_group_inventory_//')
        fi

        # Make directories

        mkdir -p local/fasta_inputs local/predictions local/inventories local/logs

        # Set variables

        input_dir="local/fasta_inputs"
        predictions_dir="local/predictions"
        log_file="local/logs/predict_with_colabfold_${GROUP_ID}.log"
        inventory_file="local/inventories/colabfold_output_group_inventory_${GROUP_ID}.tsv"

        # Set environmental variables

        export TF_FORCE_UNIFIED_MEMORY="1"
        export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"
        export XLA_PYTHON_CLIENT_ALLOCATOR="platform"
        export TF_FORCE_GPU_ALLOW_GROWTH="true"

        # Download FASTA files listed in the group inventory
        
        awk -F '\t' 'NR > 1 {print $3}' "${fasta_input_group_inventory}" | while read -r fasta_relative_path; do
        
            if [[ -z "${fasta_relative_path}" ]]; then
                continue
            fi
        
            fasta_gs_path="${SCREEN_ROOT}/${fasta_relative_path}"
            fasta_local_path="${input_dir}/$(basename "${fasta_relative_path}")"
        
            echo "Downloading FASTA: ${fasta_gs_path}" | tee -a "${log_file}"
        
            if ! gcloud storage cp "${fasta_gs_path}" "${fasta_local_path}" >> "${log_file}" 2>&1; then
                echo "FAILED_DOWNLOAD: ${fasta_gs_path}" | tee -a "${log_file}"
                continue
            fi
        
            echo "DOWNLOADED_FASTA: ${fasta_local_path}" | tee -a "${log_file}"
        
        done

        # Verify that FASTA files were downloaded
        
        num_fastas=$(find "${input_dir}" -maxdepth 1 -name "*.fasta" | wc -l)
        
        if [[ "${num_fastas}" -eq 0 ]]; then
            echo "ERROR: No FASTA files were downloaded for group ${GROUP_ID}." | tee -a "${log_file}"
            exit 1
        fi

        # Predict with ColabFold
        
        echo "Starting ColabFold predictions for group ${GROUP_ID}" | tee -a "${log_file}"
        
        if ! colabfold_batch \
            --msa-mode mmseqs2_uniref_env \
            --pair-mode unpaired_paired \
            --pair-strategy greedy \
            --random-seed 42 \
            --num-models 5 \
            --num-recycle 3 \
            --model-type auto \
            --rank auto \
            --num-ensemble 1 \
            --num-seeds 1 \
            "${input_dir}" \
            "${predictions_dir}" >> "${log_file}" 2>&1; then
        
            echo "FAILED_PREDICTION_GROUP: ${GROUP_ID}" | tee -a "${log_file}"
        
        else
        
            echo "COMPLETED_PREDICTION_GROUP: ${GROUP_ID}" | tee -a "${log_file}"
        
        fi

        # Write ColabFold output group inventory
        
        echo "colabfold_output_file" > "${inventory_file}"
        
        find "${predictions_dir}" -type f | sort | while read -r file; do
            file_from_predictions_dir="${file#${predictions_dir}/}"
            echo "predictions/${file_from_predictions_dir}" >> "${inventory_file}"
        done

        # Make a local copy without group variable in the name for WDL output

        cp "${inventory_file}" "colabfold_output_group_inventory.tsv"

        # Copy data to Google Cloud bucket
        
        if [[ -d "${predictions_dir}" ]]; then
            gcloud storage cp --recursive "${predictions_dir}"/* "${SCREEN_ROOT}/predictions/" >> "${log_file}" 2>&1
        fi
        
        gcloud storage cp "${inventory_file}" "${SCREEN_ROOT}/inventories/" >> "${log_file}" 2>&1
        gcloud storage cp "${log_file}" "${SCREEN_ROOT}/logs/"
    >>>

    output {
        File colabfold_output_group_inventory = "colabfold_output_group_inventory.tsv"
    }

    runtime {
        gpuType: "nvidia-tesla-t4"
        gpuCount: 1
        cpu: 4
        disks: "local-disk 100 HDD"
        bootDiskSizeGb: 30
        preemptible: 3
        maxRetries: 1
        memory: "16 GB"
        docker: "ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2"
    }
}
