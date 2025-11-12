version 1.0

task PredictWithColabfold {
    input {
        String query_name
        Array[File] a3m_files
    }

    command <<<
    set -euxo pipefail

    echo "GPU info:"
    nvidia-smi

    # RUN COLABFOLD PREDICTIONS
    ls -R
    output_dir="Results/~{query_name}_screen/predictions/colabfold_outputs"
    mkdir -p "${output_dir}"

    input_dir=$(dirname ~{a3m_files[0]})

    ## Recommended ColabFold GPU environment variables
    #export TF_FORCE_UNIFIED_MEMORY="1"
    #export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"
    #export XLA_PYTHON_CLIENT_ALLOCATOR="platform"
    #export TF_FORCE_GPU_ALLOW_GROWTH="true"

    ls -R
    colabfold_batch \
        --random-seed 42 \
        --num-models 5 \
        --num-recycle 3 \
        --model-type auto \
        --rank auto \
        --num-ensemble 1 \
        --num-seeds 1 \
        "${input_dir}" \
        "${output_dir}" 

    ls -R
    >>>
  
    output {
        Array[File] colabfold_output_files = glob("Results/~{query_name}_screen/predictions/colabfold_outputs/*")
    }

    runtime {
      runtime {
        cpu: 12
        memory: "85 GB"
        disks: "local-disk 2000 SSD"
        docker: 'us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/custom_build_cudabase_mmseqs2bin_colabfold:0.0.7'
        maxRetries: 2 
        zones: "us-central1-c"
        preemptible: 3
        bootDiskSizeGb: 25
        gpuCount: 1
        gpuType: "a2-highgpu-1g"
    }
}
