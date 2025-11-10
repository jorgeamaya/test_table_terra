version 1.0

task PredictWithColabfold {
    input {
        String query_name
        File a3m_file
    }

    command <<<
    set -euxo pipefail

    # RUN COLABFOLD PREDICTIONS
    ls -R
    output_dir="Results/~{query_name}_screen/predictions/colabfold_outputs"
    mkdir -p "${output_dir}"

    # Recommended ColabFold GPU environment variables
    export TF_FORCE_UNIFIED_MEMORY="1"
    export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"
    export XLA_PYTHON_CLIENT_ALLOCATOR="platform"
    export TF_FORCE_GPU_ALLOW_GROWTH="true"

    ls -R
    colabfold_batch \
        --random-seed 42 \
        --num-models 5 \
        --num-recycle 3 \
        --model-type auto \
        --rank auto \
        --num-ensemble 1 \
        --num-seeds 1 \
        "~{a3m_file}" \
        "${output_dir}" 

    ls -R
    >>>
  
    output {
        Array[File] colabfold_output_files = glob("colabfold_outputs/*")
    }
    runtime {
		gpuType: "nvidia-tesla-t4"
        gpuCount: 1
        nvidiaDriverVersion: "418.87.00"
        zones: "us-central1-c"
        cpu: 1
        disks: "local-disk 10 HDD" 
        bootDiskSizeGb: 10
        preemptible: 3
        maxRetries: 1
        memory: "12 GB"
        docker: 'us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/custom_build_cudabase_mmseqs2bin_colabfold:0.0.7'
    }
}
