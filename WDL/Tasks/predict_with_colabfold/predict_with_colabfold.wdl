version 1.0

task PredictWithColabfold {
    input {
        Array[File] predictions_files
        Array[String] size_tags
    }

    command <<<
    #https://github.com/sokrypton/ColabFold/pkgs/container/colabfold
    #https://github.com/sokrypton/ColabFold/wiki/Running-ColabFold-in-Docker
    set -euxo pipefail

    # RUN COLABFOLD PREDICTIONS ON EACH OF THE FASTA FILE GROUPS for the input size_tag
    for path in $(find "~{predictions_dir}" -mindepth 1 -maxdepth 1 -type d -name "~{size_tag}*"); do
        colabfold \
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
        "${path}" \
        "${path}"
    done

    >>>
  
    output {
        Directory predictions_dir = "~{predictions_dir}"
    }
    runtime {
        gpu: 1 # because during fasta grouping the groups are filled sequentially
        memory: select_first([
            if size_tag == "Size1KB" then "15 GiB" else
            if size_tag == "Size2KB" then "30 GiB" else
            if size_tag == "Size3KB" then "45 GiB" else
            if size_tag == "Size4KB" then "60 GiB" else
            if size_tag == "Size5KB" then "75 GiB" else
            if size_tag == "Size6KB" then "90 GiB" else
            if size_tag == "Size7KB" then "105 GiB" else
            if size_tag == "Size8KB" then "120 GiB" else
            if size_tag == "Size9KB" then "135 GiB"
        ]) 
        disks: "local-disk 10 HDD" # never used this parameter, I don't know why is needed
        bootDiskSizeGb: 10 # same as above
        preemptible: 3
        maxRetries: 1
        docker: 'ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2'
    }
}
