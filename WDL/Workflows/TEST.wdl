version 1.0

workflow TestTerraEnvironment {
        call TestEnvironment
}

task TestEnvironment {
        command <<<
                set -euox pipefail
        
                python3 <<'PY'
        import os
        
        with open("python_bucket.txt", "w") as f:
            f.write(str(os.environ.get("WORKSPACE_BUCKET", "NOT_FOUND")))
        
        with open("python_workspace.txt", "w") as f:
            f.write(str(os.environ.get("WORKSPACE_NAME", "NOT_FOUND")))
        
        with open("python_namespace.txt", "w") as f:
            f.write(str(os.environ.get("WORKSPACE_NAMESPACE", "NOT_FOUND")))
        PY
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
