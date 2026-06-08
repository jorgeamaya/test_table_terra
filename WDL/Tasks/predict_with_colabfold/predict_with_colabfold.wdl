version 1.0

task PredictWithColabfold {
    input {
				String bucket_name
				String screen_id
				File fasta_input_group_inventory
    }

    command <<<
        set -euxo pipefail

        # Prepare directory structure

        bucket_name="~{bucket_name}"
        screen_id="~{screen_id}"
        BUCKET_URI="gs://${bucket_name}"
        SCREEN_ROOT="${BUCKET_URI}/screens/${screen_id}_screen"

        mkdir -p local/predictions local/inventories local/logs

        # Set variables

        fasta_input_group_inventory_path="~{fasta_input_group_inventory}"

        fasta_list = here we need to retrieve all the paths from eg below.. 
query_name	full_subject_id	fasta_file	group	group_path	size
PF4	sehs000032n	predictions/Size1KB1/PF4_sehs000032n.fasta	Size1KB1	predictions/Size1KB1	305
PF4	sehs000032s	predictions/Size1KB1/PF4_sehs000032s.fasta	Size1KB1	predictions/Size1KB1	305
PF4	sehs000147n	predictions/Size1KB1/PF4_sehs000147n.fasta	Size1KB1	predictions/Size1KB1	239
PF4	sehs000147s	predictions/Size1KB1/PF4_sehs000147s.fasta	Size1KB1	predictions/Size1KB1	239

just get all the list of fasta files from the fasta file column...make sure you don't get the header as well. 



        # Make local VM directory for predictions
        ls -R
        output_dir="local/predictions"
        mkdir -p "${output_dir}"
    
        # Recommended ColabFold GPU environment variables
        export TF_FORCE_UNIFIED_MEMORY="1"
        export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"
        export XLA_PYTHON_CLIENT_ALLOCATOR="platform"
        export TF_FORCE_GPU_ALLOW_GROWTH="true"

        for f in fasta files do:
try: 
        ls -R
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
            "${output_dir}" 
   if fails write - prediction of fasta files....path failed...
continue..
write this to colabfold_predictions_group..name of the group....log 

        ls -R
    >>>
      
        gcloud storage cp --recursive local/predictions/* "${SCREEN_ROOT}/predictions/"        

    make the invetbnory...write it to a tsv only one column with all should it be pretier? takes more spase...we only need the lsit opf paths stored in a tsv inventory 

then for the output...you outputt he tsv right? yes...

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
            docker: 'ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2'
        }
    }
