version 1.0

task TestPredictWithColabfold {

    input {
        String bucket_name
        String screen_id
        File fasta_input_group_inventory
    }

    command <<<
        set -uo pipefail

        BUCKET_URI="gs://~{bucket_name}"
        SCREEN_ROOT="${BUCKET_URI}/screens/~{screen_id}_screen"

        GROUP_ID=$(basename "~{fasta_input_group_inventory}" .tsv | sed 's/fasta_input_group_inventory_//')

        mkdir -p local/fasta_inputs local/colabfold_outputs local/logs

        LOG_FILE="local/logs/predict_with_colabfold_${GROUP_ID}.log"
        INVENTORY_FILE="colabfold_output_group_inventory.tsv"

        echo -e "output_file" > "${INVENTORY_FILE}"

        export TF_FORCE_UNIFIED_MEMORY="1"
        export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"
        export XLA_PYTHON_CLIENT_ALLOCATOR="platform"
        export TF_FORCE_GPU_ALLOW_GROWTH="true"

        awk -F '\t' 'NR > 1 {print $3}' "~{fasta_input_group_inventory}" | while read -r fasta_file; do

            if [ -z "${fasta_file}" ]; then
                continue
            fi

            fasta_base=$(basename "${fasta_file}" .fasta)
            local_fasta="local/fasta_inputs/$(basename "${fasta_file}")"
            gcs_fasta="${SCREEN_ROOT}/${fasta_file}"

            echo "==================================================" | tee -a "${LOG_FILE}"
            echo "Processing ${gcs_fasta}" | tee -a "${LOG_FILE}"

            if ! gcloud storage cp "${gcs_fasta}" "${local_fasta}" >> "${LOG_FILE}" 2>&1; then
                echo "FAILED_DOWNLOAD ${gcs_fasta}" | tee -a "${LOG_FILE}"
                continue
            fi

            before_file="local/logs/before_${fasta_base}.txt"
            after_file="local/logs/after_${fasta_base}.txt"

            find local/colabfold_outputs -type f | sort > "${before_file}"

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
                "${local_fasta}" \
                "local/colabfold_outputs" >> "${LOG_FILE}" 2>&1; then

                echo "FAILED_PREDICTION ${gcs_fasta}" | tee -a "${LOG_FILE}"
                continue
            fi

            find local/colabfold_outputs -type f | sort > "${after_file}"

            comm -13 "${before_file}" "${after_file}" | while read -r output_file; do
                rel_output="${output_file#local/colabfold_outputs/}"
                echo "${SCREEN_ROOT}/colabfold_outputs/${rel_output}" >> "${INVENTORY_FILE}"
            done

            echo "SUCCESS ${gcs_fasta}" | tee -a "${LOG_FILE}"

        done

        if compgen -G "local/colabfold_outputs/*" > /dev/null; then
            gcloud storage cp local/colabfold_outputs/* "${SCREEN_ROOT}/colabfold_outputs/" >> "${LOG_FILE}" 2>&1
        else
            echo "WARNING: no ColabFold outputs produced for group ${GROUP_ID}" | tee -a "${LOG_FILE}"
        fi

        gcloud storage cp "${INVENTORY_FILE}" "${SCREEN_ROOT}/inventories/colabfold_output_group_inventory_${GROUP_ID}.tsv" >> "${LOG_FILE}" 2>&1
        gcloud storage cp "${LOG_FILE}" "${SCREEN_ROOT}/logs/" >> "${LOG_FILE}" 2>&1
    >>>

    output {
        File colabfold_output_group_inventory = "colabfold_output_group_inventory.tsv"
    }

    runtime {
        gpuType: "nvidia-tesla-t4"
        gpuCount: 1
        cpu: 1
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 30
        preemptible: 3
        maxRetries: 1
        memory: "12 GB"
        docker: "ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2"
    }
}
