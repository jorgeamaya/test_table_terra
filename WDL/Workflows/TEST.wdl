version 1.0

workflow TestTerraEnvironment {
        call TestEnvironment
}

task TestEnvironment {
        command <<<
                set -euox pipefail

                env | sort > env.txt

                echo "${WORKSPACE_BUCKET:-NOT_FOUND}" > workspace_bucket.txt
                echo "${WORKSPACE_NAME:-NOT_FOUND}" > workspace_name.txt
                echo "${WORKSPACE_NAMESPACE:-NOT_FOUND}" > workspace_namespace.txt
        >>>

        output {
                File env_file = "env.txt"
                File workspace_bucket = "workspace_bucket.txt"
                File workspace_name = "workspace_name.txt"
                File workspace_namespace = "workspace_namespace.txt"
        }

        runtime {
                docker: "ubuntu:22.04"
                cpu: 1
                memory: "1 GB"
                disks: "local-disk 10 HDD"
        }
}
