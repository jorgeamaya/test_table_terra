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
		BUCKET_URI="gs://~{bucket_name}"
		SCREEN_ROOT="${BUCKET_URI}/screens/~{screen_id}"

		mkdir -p local/inputs local/fasta_inputs local/inventories local/logs

    # Validate input and proceed if validation checks are fully met

    ### we work here now

    # Write validated inputs
		echo '{"screen_id":"~{screen_id}","query_name":"~{query_name}","mode":"dummy"}' > local/inputs/screen_input.json

		cat > local/inputs/query.fasta <<EOF
>~{query_name}
~{query_sequence}
EOF

		echo "dummy native sequences" > local/inputs/subject_native_sequences.fasta
		echo "dummy scrambled sequences" > local/inputs/subject_scrambled_sequences.fasta
		echo -e "subject_id\tgene\tdescription" > local/inputs/subject_proteome_dictionary.tsv
		echo -e "P00001\tDUMMY1\tDummy subject 1" >> local/inputs/subject_proteome_dictionary.tsv
		echo -e "P00002\tDUMMY2\tDummy subject 2" >> local/inputs/subject_proteome_dictionary.tsv

		echo ">~{query_name}__P00001_native" > local/fasta_inputs/~{query_name}__P00001_native.fasta
		echo "MSEQUENCEONE" >> local/fasta_inputs/~{query_name}__P00001_native.fasta

		echo ">~{query_name}__P00001_scrambled" > local/fasta_inputs/~{query_name}__P00001_scrambled.fasta
		echo "MSEQUENCESCRAMBLEDONE" >> local/fasta_inputs/~{query_name}__P00001_scrambled.fasta

		echo ">~{query_name}__P00002_native" > local/fasta_inputs/~{query_name}__P00002_native.fasta
		echo "MSEQUENCETWO" >> local/fasta_inputs/~{query_name}__P00002_native.fasta

		echo ">~{query_name}__P00002_scrambled" > local/fasta_inputs/~{query_name}__P00002_scrambled.fasta
		echo "MSEQUENCESCRAMBLEDTWO" >> local/fasta_inputs/~{query_name}__P00002_scrambled.fasta

		echo -e "fasta_id\tfasta_file" > fasta_input_inventory.tsv
		for f in local/fasta_inputs/*.fasta; do
			b=$(basename "${f}")
			id="${b%.fasta}"
			echo -e "${id}\t${SCREEN_ROOT}/fasta_inputs/${b}" >> fasta_input_inventory.tsv
		done

		echo -e "group_id\tfasta_file" > fasta_input_group_inventory_001.tsv
		echo -e "group_001\t${SCREEN_ROOT}/fasta_inputs/~{query_name}__P00001_native.fasta" >> fasta_input_group_inventory_001.tsv
		echo -e "group_001\t${SCREEN_ROOT}/fasta_inputs/~{query_name}__P00001_scrambled.fasta" >> fasta_input_group_inventory_001.tsv

		echo -e "group_id\tfasta_file" > fasta_input_group_inventory_002.tsv
		echo -e "group_002\t${SCREEN_ROOT}/fasta_inputs/~{query_name}__P00002_native.fasta" >> fasta_input_group_inventory_002.tsv
		echo -e "group_002\t${SCREEN_ROOT}/fasta_inputs/~{query_name}__P00002_scrambled.fasta" >> fasta_input_group_inventory_002.tsv

		echo "dummy prepare fasta files completed" > local/logs/prepare_fasta_files.log

		gcloud storage cp local/inputs/* "${SCREEN_ROOT}/inputs/"
		gcloud storage cp local/fasta_inputs/* "${SCREEN_ROOT}/fasta_inputs/"
		gcloud storage cp fasta_input_inventory.tsv "${SCREEN_ROOT}/inventories/fasta_input_inventory.tsv"
		gcloud storage cp fasta_input_group_inventory_001.tsv "${SCREEN_ROOT}/inventories/fasta_input_group_inventory_001.tsv"
		gcloud storage cp fasta_input_group_inventory_002.tsv "${SCREEN_ROOT}/inventories/fasta_input_group_inventory_002.tsv"
		gcloud storage cp local/logs/prepare_fasta_files.log "${SCREEN_ROOT}/logs/prepare_fasta_files.log"
	>>>

	output {
		File fasta_input_inventory = "fasta_input_inventory.tsv"
		Array[File] fasta_input_group_inventories = ["fasta_input_group_inventory_001.tsv", "fasta_input_group_inventory_002.tsv"]
	}

	runtime {
		docker: "google/cloud-sdk:slim"
	}
}
