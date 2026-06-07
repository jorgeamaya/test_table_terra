version 1.0

task TestPredictWithColabfold {

	input {
		String bucket_name
		String screen_id
		File fasta_input_group_inventory
	}

	command <<<
		set -e

		BUCKET_URI="gs://~{bucket_name}"
		SCREEN_ROOT="${BUCKET_URI}/screens/~{screen_id}"

		GROUP_ID=$(basename "~{fasta_input_group_inventory}" .tsv | sed 's/fasta_input_group_inventory_//')

		mkdir -p local/colabfold_outputs local/logs

		echo -e "output_id\toutput_file" > colabfold_output_group_inventory.tsv

		while IFS=$'\t' read -r group_id fasta_file; do
			if [ "${group_id}" = "group_id" ]; then
				continue
			fi

			fasta_base=$(basename "${fasta_file}" .fasta)

			for rank in 001 002 003 004 005; do
				pdb_file="${fasta_base}_unrelaxed_rank_${rank}_model_1_seed_000.pdb"
				json_file="${fasta_base}_scores_rank_${rank}.json"

				echo "MOCK PDB ${fasta_base} rank ${rank}" > "local/colabfold_outputs/${pdb_file}"
				echo "{\"mock\": true, \"fasta_id\": \"${fasta_base}\", \"rank\": \"${rank}\"}" > "local/colabfold_outputs/${json_file}"

				echo -e "${fasta_base}_pdb_rank_${rank}\t${SCREEN_ROOT}/colabfold_outputs/${pdb_file}" >> colabfold_output_group_inventory.tsv
				echo -e "${fasta_base}_scores_rank_${rank}\t${SCREEN_ROOT}/colabfold_outputs/${json_file}" >> colabfold_output_group_inventory.tsv
			done

			echo "MOCK PAE ${fasta_base}" > "local/colabfold_outputs/${fasta_base}_pae.png"
			echo "MOCK COVERAGE ${fasta_base}" > "local/colabfold_outputs/${fasta_base}_coverage.png"
			echo "{\"mock_config\": true}" > "local/colabfold_outputs/${fasta_base}_config.json"

			echo -e "${fasta_base}_pae_png\t${SCREEN_ROOT}/colabfold_outputs/${fasta_base}_pae.png" >> colabfold_output_group_inventory.tsv
			echo -e "${fasta_base}_coverage_png\t${SCREEN_ROOT}/colabfold_outputs/${fasta_base}_coverage.png" >> colabfold_output_group_inventory.tsv
			echo -e "${fasta_base}_config_json\t${SCREEN_ROOT}/colabfold_outputs/${fasta_base}_config.json" >> colabfold_output_group_inventory.tsv

		done < "~{fasta_input_group_inventory}"

		echo "dummy predict with colabfold completed for group ${GROUP_ID}" > local/logs/predict_with_colabfold_${GROUP_ID}.log

		gcloud storage cp local/colabfold_outputs/* "${SCREEN_ROOT}/colabfold_outputs/"
		gcloud storage cp colabfold_output_group_inventory.tsv "${SCREEN_ROOT}/inventories/colabfold_output_group_inventory_${GROUP_ID}.tsv"
		gcloud storage cp local/logs/predict_with_colabfold_${GROUP_ID}.log "${SCREEN_ROOT}/logs/"
	>>>

	output {
		File colabfold_output_group_inventory = "colabfold_output_group_inventory.tsv"
	}

	runtime {
		docker: "google/cloud-sdk:slim"
	}
}
