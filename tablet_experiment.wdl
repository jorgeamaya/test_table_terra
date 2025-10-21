version 1.0

task MakeNewTable {
  input {
    Array[String] entity_ids
    Array[File] outputs
    String table_name
  }

  command <<<
    set -eux
    out_tsv="~{table_name}.tsv"
    paste <(printf "%s\n" "${sep='\n' entity_ids}") <(printf "%s\n" "${sep='\n' outputs}") \
      | awk 'BEGIN {OFS="\t"; print "entity:record_id", "output"} {print $1, $2}' > ${out_tsv}

    # Copy directly to Terra's monitored path so it becomes available immediately
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
  }

  call MakeNewTable {
    input:
      entity_ids = sample_ids,
      outputs = result_files,
      table_name = "predictions"
  }

  output {
    File terra_table = MakeNewTable.table_tsv
  }
}
