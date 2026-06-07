version 1.0

task TestInventoryColabfoldOutputs {

	input {
		String bucket_name
		String screen_id
		Array[File] colabfold_output_group_inventories
	}

	command <<<
		set -e

		BUCKET_URI="gs://~{bucket_name}"
		SCREEN_ROOT="${BUCKET_URI}/screens/~{screen_id}"

		echo -e "output_id\toutput_file" > colabfold_output_inventory.tsv

		for inventory in ~{sep=" " colabfold_output_group_inventories}; do
			tail -n +2 "${inventory}" >> colabfold_output_inventory.tsv
		done

		echo "dummy inventory colabfold outputs completed" > inventory_colabfold_outputs.log

		gcloud storage cp colabfold_output_inventory.tsv "${SCREEN_ROOT}/inventories/colabfold_output_inventory.tsv"
		gcloud storage cp inventory_colabfold_outputs.log "${SCREEN_ROOT}/logs/inventory_colabfold_outputs.log"
	>>>

	output {
		File colabfold_output_inventory = "colabfold_output_inventory.tsv"
	}

	runtime {
		docker: "google/cloud-sdk:slim"
	}
}
