version 1.0

task PrepareFastaFiles {
    input {
        String query_sequence 
        String query_name 
        File subject_native_sequences_file 
        File subject_scrambled_sequences_file 
    }

    command <<<
    set -euxo pipefail

    screen_dir="Results/~{query_name}_screen"
    mkdir -p ${screen_dir}
    predictions_dir="${screen_dir}/predictions"
    mkdir -p "${predictions_dir}"
    query_protein_path="${screen_dir}/query_~{query_name}.tsv"
    echo -e "query_name\tquery_sequence\n~{query_name}\t~{query_sequence}" > ${query_protein_path}

    #DEBUGGING
    ls -R 
    conda info --envs
    python -c "import protbindscreen" && echo "Available" || echo "Not available"

    python3 -c "from pathlib import Path; from protbindscreen.submission.submission_helper import make_fasta_files; count = make_fasta_files(Path('${query_protein_path}'), Path('~{subject_native_sequences_file}'), Path('~{subject_scrambled_sequences_file}'), Path('${predictions_dir}'))"

    ls -R
    >>>

    output {
        Array[File] colabfold_input_files = glob("Results/~{query_name}_screen/predictions/*.fasta")
    }

    runtime {
        bootDiskSizeGb: 10
        disks: "local-disk 10 HDD" 
        cpu: 1
        memory: "8 GiB" 
        preemptible: 3
        maxRetries: 1
        docker: 'us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/protbindscreen:0.0.7'
    }
}
