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

		mkdir -p local/inputs local/fasta_inputs local/inventories local/logs

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

    # Copy subject proteome dictionary to local inputs
		subject_proteome_dictionary_file=""

		for file in ${subject_proteome_datasets}; do
			file_name="$(basename "${file}")"

			if [[ "${file_name}" == *subject_proteome_dictionary.tsv ]]; then
				subject_proteome_dictionary_file="${file}"
				break
			fi
		done

		if [[ -z "${subject_proteome_dictionary_file}" ]]; then
			echo "ERROR: Could not find subject_proteome_dictionary.tsv in subject_proteome_datasets." >&2
			exit 1
		fi

		cp "${subject_proteome_dictionary_file}" \
			local/inputs/subject_proteome_dictionary.tsv

    # Write query TSV to local inputs
		query_tsv="local/inputs/query_${query_name}.tsv"

		printf "query_name\tquery_sequence\n" > "${query_tsv}"
		printf "%s\t%s\n" "${query_name}" "${query_sequence}" >> "${query_tsv}"

	# Here the logic for fasta making and grouping to follow. 


		gcloud storage cp local/inputs/* "${SCREEN_ROOT}/inputs/"
		gcloud storage cp local/fasta_inputs/* "${SCREEN_ROOT}/predictions/"
		gcloud storage cp local/inventories/* "${SCREEN_ROOT}/inventories/"
		gcloud storage cp local/logs/* "${SCREEN_ROOT}/logs/"
	>>>

	output {
		File fasta_input_inventory = "local/inventories/fasta_input_inventory.tsv"
		Array[File] fasta_input_group_inventories = glob("local/inventories/fasta_input_group_inventory_*.tsv")
	}

	runtime {
		docker: "google/cloud-sdk:slim"
	}
}
