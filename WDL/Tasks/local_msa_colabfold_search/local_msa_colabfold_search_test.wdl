version 1.0

task LocalMSAColabfoldSearch {
    input {
        String query_name
        Array[File] colabfold_input_files
        String colabfold_db_path
    }

    command <<<
    set -euxo pipefail

    ## GETTING FASTA INPUT DIRECTORY AND LISTING FILES
    echo "Files:"
    for f in ~{sep=" " colabfold_input_files}; do
        echo "$f"
    done

    # Get parent directory of the first file
    input_dir=$(dirname ~{colabfold_input_files[0]})
    echo "Parent directory: ${input_dir}"

    # List all input files
    ls -lh ${input_dir}

    # PREPARING COLABFOLD DATABASE DIRECTORY

    # This is what I got so far, but I am completelly annoyed and I don't see it as a solution...or it can be a last or temporary one.
    # Anyway, how much disk space we need to request at runtime? 2TB SSD? the db is 1.4-1.6TB. The a3m files are also something.
    # TO BE discussed - https://docs.cloud.google.com/storage/docs/cloud-storage-fuse/overview 
    # GC seems to offer solutions for mounting buckets as file systems, but I don't know if it is possible on Terra/Cromwell. 
    # The Fuse may not be compatible because Cloud Storage FUSE is not POSIX compliant
    # Other options we can explore are https://cloud.google.com/filestore?hl=en and "persistent disk"???
    # I cannot test locally anything that involved copying into the VM 1.6T because I don't have space.

    colabfold_db_dir="colabfold_db"
    mkdir -p "${colabfold_db_dir}"
    echo "Downloading ColabFold database to ${colabfold_db_dir} ..."
    gcloud storage ls ~{colabfold_db_path}
    gcloud storage cp -r ~{colabfold_db_path} "${colabfold_db_dir}/"


    # MAKE OUTPUT DIRECTORY
    ls -R
    output_dir="Results/~{query_name}_screen/predictions/a3m_files"
    mkdir -p "${output_dir}"

    # # Recommended ColabFold GPU environment variables
    # export TF_FORCE_UNIFIED_MEMORY="1"
    # export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"
    # export XLA_PYTHON_CLIENT_ALLOCATOR="platform"
    # export TF_FORCE_GPU_ALLOW_GROWTH="true"


    # RUN MSAs WITH COLABFOLD SEARCH
    # colabfold_search "/workspace/localmsa_4input" "/workspace/databases" "localmsa_4output" --gpu 1
    ls -R
    colabfold_search \
        "${input_dir}" \
        "${colabfold_db_dir}" \
        "${output_dir}" \
        --gpu 1

    ls -R
    >>>
  
    output {
        Array[File] colabfold_a3m_files = glob("Results/~{query_name}_screen/predictions/a3m_files/~{query_name}*.a3m")
    }
    runtime {
	    gpuType: "nvidia-tesla-t4"
        gpuCount: 1
        nvidiaDriverVersion: "418.87.00"
        zones: "us-central1-c"
        cpu: 1
        disks: "local-disk 2000 SSD" 
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        memory: "256 GB"
        docker: 'us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/custom_build_cudabase_mmseqs2bin_colabfold:0.0.7'
    }
}
