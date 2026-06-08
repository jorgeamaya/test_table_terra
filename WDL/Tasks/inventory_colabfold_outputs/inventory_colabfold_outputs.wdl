version 1.0

task InventoryColabfoldOutputs {

    input {
        String bucket_name
        String screen_id
        Array[File] colabfold_output_group_inventories
    }

    command <<<
        set -euxo pipefail

        bucket_name="~{bucket_name}"
        screen_id="~{screen_id}"

        BUCKET_URI="gs://${bucket_name}"
        SCREEN_ROOT="${BUCKET_URI}/screens/${screen_id}_screen"

        mkdir -p local/inventories local/logs

        output_inventory="local/inventories/colabfold_output_inventory.tsv"
        log_file="local/logs/inventory_colabfold_outputs.log"

        echo "colabfold_output_file" > "${output_inventory}"

        for inventory in ~{sep=' ' colabfold_output_group_inventories}; do
            tail -n +2 "${inventory}" >> "${output_inventory}"
        done

        gcloud storage cp "${output_inventory}" "${SCREEN_ROOT}/inventories/" >> "${log_file}" 2>&1
        gcloud storage cp "${log_file}" "${SCREEN_ROOT}/logs/"
    >>>

    output {
        File colabfold_output_inventory = "local/inventories/colabfold_output_inventory.tsv"
    }

    runtime {
        bootDiskSizeGb: 10
        disks: "local-disk 10 HDD"
        cpu: 1
        memory: "4 GiB"
        preemptible: 3
        maxRetries: 1
        docker: "google/cloud-sdk:slim"
    }
}
