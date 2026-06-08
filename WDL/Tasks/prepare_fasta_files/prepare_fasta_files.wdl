version 1.0

task TestPrepareFastaFiles {

	input {
		String bucket_name
		String screen_id
		String query_name
		String query_sequence
		Array[File] subject_proteome_datasets
	}

	command <<<
		set -e

    # Prepare directory structure

		bucket_name="~{bucket_name}"
		screen_id="~{screen_id}"
		BUCKET_URI="gs://${bucket_name}"
		SCREEN_ROOT="${BUCKET_URI}/screens/${screen_id}_screen"

		mkdir -p local/inputs local/predictions local/inventories local/logs

	# Set variables
	
		mode="submit"
		query_name="~{query_name}"
		query_sequence="~{query_sequence}"
		subject_proteome_datasets="~{sep=' ' subject_proteome_datasets}"
		log_dir="local/logs"
	
	# Validate input and write validation logs
	
		python -m protbindscreen.protbindscreen_runner \
			--mode "${mode}" \
			--query_name "${query_name}" \
			--query_sequence "${query_sequence}" \
			--subject_proteome_datasets ${subject_proteome_datasets} \
			--log_dir "${log_dir}"
	
    # Write validated inputs

		cat > local/inputs/screen_input.json << EOF
{
  "screen_id": "${screen_id}",
  "query_name": "${query_name}",
  "query_sequence": "${query_sequence}",
  "subject_proteome_datasets": "${subject_proteome_datasets}"
}
EOF

    # Copy subject proteome dataset files to local inputs

		for file in ${subject_proteome_datasets}; do
			cp "${file}" "local/inputs/$(basename "${file}")"
		done

    # Define local subject proteome dataset paths

		subject_proteome_dictionary_file="local/inputs/subject_proteome_dictionary.tsv"
		subject_native_sequences_file="local/inputs/subject_proteome_native_seq.tsv"
		subject_scrambled_sequences_file="local/inputs/subject_proteome_scrambled_seq.tsv"

		if [[ ! -f "${subject_proteome_dictionary_file}" ]]; then
			echo "ERROR: Missing local subject_proteome_dictionary.tsv." >&2
			exit 1
		fi

		if [[ ! -f "${subject_native_sequences_file}" ]]; then
			echo "ERROR: Missing local subject_proteome_native_seq.tsv." >&2
			exit 1
		fi

		if [[ ! -f "${subject_scrambled_sequences_file}" ]]; then
			echo "ERROR: Missing local subject_proteome_scrambled_seq.tsv." >&2
			exit 1
		fi

    # Write query TSV to local inputs

		query_tsv="local/inputs/query_${query_name}.tsv"

		printf "query_name\tquery_sequence\n" > "${query_tsv}"
		printf "%s\t%s\n" "${query_name}" "${query_sequence}" >> "${query_tsv}"

	# Generate FASTA prediction inputs and inventories

		python -m protbindscreen.submission.screen_submission \
			--query_protein_path "${query_tsv}" \
			--subject_native_sequences_path "${subject_native_sequences_file}" \
			--subject_scrambled_sequences_path "${subject_scrambled_sequences_file}" \
			--predictions_dir "local/predictions" \
			--inventory_dir "local/inventories" \
			--log_dir "${log_dir}"

		gcloud storage cp local/inputs/* "${SCREEN_ROOT}/inputs/"
		gcloud storage cp --recursive local/predictions/* "${SCREEN_ROOT}/predictions/"
		gcloud storage cp local/inventories/* "${SCREEN_ROOT}/inventories/"
		gcloud storage cp local/logs/* "${SCREEN_ROOT}/logs/"
	>>>

	output {
		File fasta_input_inventory = "local/inventories/fasta_input_inventory.tsv"
		Array[File] fasta_input_group_inventories = glob("local/inventories/fasta_input_group_inventory_*.tsv")
	}

	runtime {
        bootDiskSizeGb: 10
        disks: "local-disk 10 HDD" 
        cpu: 1
        memory: "16 GiB" 
        preemptible: 3
        maxRetries: 1
		docker: "google/cloud-sdk:slim"
	}
}
