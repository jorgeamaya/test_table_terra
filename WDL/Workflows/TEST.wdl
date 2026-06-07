version 1.0

import "tasks/TestBucket.wdl" as Tasks

workflow TestEnvironment {

    input {
        String bucket_name
    }

    call Tasks.TestBucket {
        input:
            bucket_name = bucket_name
    }

    output {
        String bucket_name_out = Tasks.TestBucket.bucket_name_out
        String bucket_uri_out = Tasks.TestBucket.bucket_uri_out
    }
}
