version 1.0

task PrepareFastaGroups {
    input {
        String query_sequence #= "BMP6"
        String query_name #= "wdl_test_proteome"
        File subject_native_sequences_file #= "subject_proteome_native_seq.tsv"
        File subject_scrambled_sequences_file #= "subject_proteome_scrambled_seq.tsv"
        File subject_proteome_dictionary_file #= "subject_proteome_dictionary.tsv"
    }

    command <<<
    set -euxo pipefail

    # Create output screen directory and predictions subdirectory
    screen_dir="Results/~{query_name}_screen"
    mkdir -p ${screen_dir}
    predictions_dir="${screen_dir}/predictions"
    mkdir -p "${predictions_dir}"
    query_protein_path="${screen_dir}/query_~{query_name}.tsv"
    echo -e "query_name\tquery_sequence\n~{query_name}\t~{query_sequence}" > ${query_protein_path}

    cp ~{subject_proteome_dictionary_file} ${screen_dir}

    #DEBUGGING
    ls -R 
    conda info --envs
    python -c "import protbindscreen" && echo "Available" || echo "Not available"

    # GENERATE FASTA FILES by combining the query with each of the subject sequences
    # The submission_helper.make_fasta_files should be available from the docker image being used
    # The arguments are:
    # 1. path to the query protein tsv file
    # 2. path to the subject native sequences tsv file
    # 3. path to the subject scrambled sequences tsv file
    # 4. path to the predictions output directory where the fasta files will be written
    # This function also return the number of generated fasta files, which is captured here
    
    fasta_file_count=$(python3 -c "from pathlib import Path; from protbindscreen.submission.submission_helper import make_fasta_files; count = make_fasta_files(Path('${query_protein_path}'), Path('~{subject_native_sequences_file}'), Path('~{subject_scrambled_sequences_file}'), Path('${predictions_dir}')); print(count)")

    echo "$fasta_file_count" > ${predictions_dir}/fasta_file_count.txt

    # GROUP FASTA FILES INTO DIRECTORIES FOR BATCHED PREDICTIONS
    # The arguments are:
    # 1. path to the predictions output directory where the fasta files were just written
    # This function returns:
    # - the number of created groups (directories), 
    # - a placement data file path, 
    # - a count of trailing fasta files (for documentation purposes, but should be 0)
    
    read fasta_groups_count placement_data_path trailing_fasta_files <<< $(python3 -c "from pathlib import Path; from protbindscreen.submission.submission_helper import group_fasta_files; fasta_groups_count, placement_data_path, trailing_fasta_files = group_fasta_files(Path('${predictions_dir}')); print(f'{fasta_groups_count} {placement_data_path} {trailing_fasta_files}')")

    echo "$fasta_groups_count" > ${predictions_dir}/fasta_groups_count.txt
    echo "$trailing_fasta_files" > ${predictions_dir}/fasta_trailing_files_count.txt

    # Without adding the tag in the file name, the grouping is completelly unnecessary. 
    # Grouping might be completelly unnecesarry anyway
    # because task 2 will only see the files copied to it's own inputs in it's own sandbox
    # without knowing anything about the dir structure in task 1

    # for f in Results/~{query_name}_screen/predictions/*/*; do
    #   dir=$(basename "$(dirname "$f")")
    #   base=$(basename "$f")
    #   mv "$f" "$(dirname "$f")/${base%.*}_${dir}.${base##*.}"
    #done

    ls -R
    >>>

  
    output {
        Array[File] colabfold_input_files = glob("Results/~{query_name}_screen/predictions/*/*.fasta")
    }

    runtime {
        cpu: 1
        memory: "8 GiB" 
        disks: "local-disk 10 HDD" 
        bootDiskSizeGb: 10 
        preemptible: 3
        maxRetries: 1
        docker: 'us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/protbindscreen:v0.0.3'
    }
}
