version 1.0

task AnalyzeScreen {
    input {
		String bucket_name
		String screen_id
		File fasta_input_inventory 
		File colabfold_output_inventory
		Array[File] subject_proteome_datasets 
		String analysis_id
		Int query_len 
		Array[String] aa_ranges_i
		Array[String] aa_ranges_j
		Int pae_threshold
    }

    command <<<
        set -euxo pipefail

        # Prepare directory structure

        bucket_name="~{bucket_name}"
        analysis_id="~{analysis_id}"
        BUCKET_URI="gs://${bucket_name}"
        ANALYSIS_ROOT="${BUCKET_URI}/analyses/${analysis_id}_screen"
    
        screen_id="~{screen_id}"
        SCREEN_ROOT="${BUCKET_URI}/screens/${screen_id}_screen"

        mkdir -p local/inputs local/outputs local/inventories local/logs

        # Set variables

        mode="analyze"
        query_len="~{query_len}"
        pae_threshold="~{pae_threshold}"

        input_dir="local/inputs"
        output_dir="local/outputs"
        log_dir="local/logs"
        log_file="${log_dir}/analyze_screen.log"

        fasta_input_inventory="~{fasta_input_inventory}"
        colabfold_output_inventory="~{colabfold_output_inventory}"
    
        # Build analysis_matrices from ranges

        analysis_matrices_file="${output_dir}/analysis_matrices.json"

        python <<PY
import json

aa_ranges_i = "~{sep=' ' aa_ranges_i}".split()
aa_ranges_j = "~{sep=' ' aa_ranges_j}".split()

if len(aa_ranges_i) != len(aa_ranges_j):
    raise ValueError(
        f"aa_ranges_i and aa_ranges_j must have the same length. "
        f"Got {len(aa_ranges_i)} and {len(aa_ranges_j)}."
    )

analysis_matrices = {}

for idx, (range_i, range_j) in enumerate(zip(aa_ranges_i, aa_ranges_j)):
    matrix_key = "required" if idx == 0 else f"optional_{idx}"
    matrix_name = f"matrix_{matrix_key}"

    analysis_matrices[matrix_name] = {
        "aa_ranges_i": range_i,
        "aa_ranges_j": range_j
    }

    if matrix_key != "required":
        analysis_matrices[matrix_name]["include"] = "true"

with open("${analysis_matrices_file}", "w") as f:
    json.dump(analysis_matrices, f, indent=4)
PY

        if [[ ! -s "${analysis_matrices_file}" ]]; then
            echo "ERROR: Failed to create analysis_matrices.json"
            exit 1
        fi

        # Validate input and write validation logs
        
        python -m protbindscreen.protbindscreen_runner \
            --mode "${mode}" \
            --query_len "${query_len}" \
            --analysis_matrices "${analysis_matrices_file}" \
            --pae_threshold "${pae_threshold}" \
            --log_dir "${log_dir}"


        # Copy to local VM the proteome dictionary

        subject_proteome_dictionary_file=""
        
        for file in ~{sep=' ' subject_proteome_datasets}; do
            base="$(basename "${file}")"
        
            if [[ "${base}" == *subject_proteome_dictionary.tsv ]]; then
                cp "${file}" "${input_dir}/${base}"
                subject_proteome_dictionary_file="${input_dir}/${base}"
            fi
        done
        
        if [[ -z "${subject_proteome_dictionary_file}" ]]; then
            echo "ERROR: No subject_proteome_dictionary.tsv found in subject_proteome_datasets."
            exit 1
        fi
        
        if [[ ! -s "${subject_proteome_dictionary_file}" ]]; then
            echo "ERROR: Subject proteome dictionary file is missing or empty: ${subject_proteome_dictionary_file}"
            exit 1
        fi

        # Download FASTA files listed in the inventory
              
        awk -F '\t' 'NR > 1 {print $3}' "${fasta_input_inventory}" | while read -r fasta_relative_path; do
        
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

        # Download ColabFold output files listed in the inventory
        # Keep only files directly under predictions/, skip nested MSA/pairgreedy directories.
        
        awk -F '\t' 'NR > 1 {print $1}' "${colabfold_output_inventory}" | while read -r colabfold_output_relative_path; do
        
            if [[ -z "${colabfold_output_relative_path}" ]]; then
                continue
            fi
        
            # Skip nested output files, e.g. predictions/*_env/... or predictions/*_pairgreedy/...
            if [[ "${colabfold_output_relative_path}" == */*/* ]]; then
                echo "SKIPPED_NESTED_COLABFOLD_OUTPUT: ${colabfold_output_relative_path}" | tee -a "${log_file}"
                continue
            fi
        
            colabfold_output_gs_path="${SCREEN_ROOT}/${colabfold_output_relative_path}"
            colabfold_output_local_path="${input_dir}/$(basename "${colabfold_output_relative_path}")"
        
            echo "Downloading ColabFold output: ${colabfold_output_gs_path}" | tee -a "${log_file}"
        
            if ! gcloud storage cp "${colabfold_output_gs_path}" "${colabfold_output_local_path}" >> "${log_file}" 2>&1; then
                echo "FAILED_DOWNLOAD: ${colabfold_output_gs_path}" | tee -a "${log_file}"
                continue
            fi
        
            echo "DOWNLOADED_COLABFOLD_OUTPUT: ${colabfold_output_local_path}" | tee -a "${log_file}"
        
        done

        # Verify that files were downloaded
        
        num_files=$(find "${input_dir}" -maxdepth 1 -type f | wc -l)
        
        if [[ "${num_files}" -eq 0 ]]; then
            echo "ERROR: No files were downloaded."
            exit 1
        fi

        # Run ProtBindScreen analyze module

        python -m protbindscreen.submission.screen_analysis \
            --input_dir "${input_dir}" \
            --output_dir "${output_dir}" \
            --analysis_name "${analysis_id}" \
            --query_len "${query_len}" \
            --analysis_matrices "${analysis_matrices_file}" \
            --subject_proteome_dictionary "${subject_proteome_dictionary_file}" \
            --pae_threshold "${pae_threshold}" \
            --log_dir "${log_dir}"

        # Copy to GCS
 
    output {
        Array[File] analysis_output_files = glob("Results/~{screen_name}/**/**/**")
    }

    runtime {
        cpu: 1 
        memory: "8 GiB" 
        disks: "local-disk 10 HDD" 
        bootDiskSizeGb: 10 
        preemptible: 3
        maxRetries: 1
        docker: 'us-central1-docker.pkg.dev/lithe-aileron-498218-r4/private-gar-protbindscreen-docker-images/protbindscreen:0.0.7'
    }
}
