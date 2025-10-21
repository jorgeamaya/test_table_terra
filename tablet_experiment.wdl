version 1.0

task MakeNewTable {
  input {
    Array[String] entity_ids
    Array[File] outputs
    String table_name
    String workspace_name
    String workspace_bucket
  }

  command <<<
    set -eux
    out_tsv="~{table_name}.tsv"
    #paste <(printf "%s\n" "${sep='\n' entity_ids}") <(printf "%s\n" "${sep='\n' outputs}") \
    #  | awk 'BEGIN {OFS="\t"; print "entity:record_id", "output"} {print $1, $2}' > ${out_tsv}
    echo -e "entity:record_id\toutput" > ${out_tsv}
    echo -e "sample1\t/path/to/output1" >> ${out_tsv}
    echo -e "sample2\t/path/to/output2" >> ${out_tsv}
    echo -e "sample3\t/path/to/output3" >> ${out_tsv}
    echo -e "sample4\t/path/to/output4" >> ${out_tsv}
    echo -e "sample5\t/path/to/output5" >> ${out_tsv}

    gsutil cp terra_table.tsv ~{workspace_bucket}/terra_table.tsv
    
    # Import the table directly into the workspace
    terra workspace set ~{workspace_name}
    terra data import-entity-table \
      --table sample \
      --file terra_table.tsv

  >>>

  output {
    File table_tsv = "~{table_name}.tsv"
  }

  runtime {
    docker: 'jorgeamaya/protbindscreen:v0.0.1'
    memory: "1 GiB"
    cpu: 1
  }
}

workflow WriteTableExample {
  input {
    Array[String] sample_ids
    Array[File] result_files
    String workspace_name
    String workspace_bucket

  }

  call MakeNewTable {
    input:
      entity_ids = sample_ids,
      outputs = result_files,
      table_name = "predictions",
      workspace_name = workspace_name,
      workspace_bucket = workspace_bucket
  }

  output {
    File terra_table = MakeNewTable.table_tsv
  }
}
