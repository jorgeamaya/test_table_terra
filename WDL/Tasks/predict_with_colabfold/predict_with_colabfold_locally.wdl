version 1.0

task PredictWithColabfoldLocally {
    input {
        File fasta_file
    }

    command <<<
    #https://github.com/sokrypton/ColabFold/pkgs/container/colabfold
    #https://github.com/sokrypton/ColabFold/wiki/Running-ColabFold-in-Docker
    set -euxo pipefail

    # RUN COLABFOLD PREDICTIONS

    output_dir="colabfold_outputs"
    mkdir -p "${output_dir}"

    # Recommended ColabFold GPU environment variables
    export TF_FORCE_UNIFIED_MEMORY="1"
    export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"
    export XLA_PYTHON_CLIENT_ALLOCATOR="platform"
    export TF_FORCE_GPU_ALLOW_GROWTH="true"

    # Start continuous GPU logging in background
    nvidia-smi --query-gpu=timestamp,memory.used,memory.total,utilization.gpu \
           --format=csv,nounits,noheader --loop-ms=1000 > gpu_usage.csv &
    GPU_MON_PID=$!

    colabfold_batch \
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
        "~{fasta_file}" \
        "${output_dir}" 2>&1 | tee -a colabfold_log.txt
        kill $GPU_MON_PID || true
    >>>
  
    output {
        Array[File] colabfold_outputs = glob("colabfold_outputs/*")
    }
    runtime {
        gpu: 1
        cpu: 4
        disks: "local-disk 10 HDD" 
        bootDiskSizeGb: 10
        preemptible: 3
        maxRetries: 1
        memory: "12 GB"
        docker: 'ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2'
    }
}
