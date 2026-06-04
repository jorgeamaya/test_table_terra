# Internal README file for devs

`WORKFLOWS` will be named in uppercase.

`tasks` will be named in lowercase. 

# Source Code Privacy Strategy

During the development, testing, and pre-release phases of ProtBindScreen, workflow deployment follows a GitHub → Dockstore → Terra architecture designed to enable workflow execution while maintaining the privacy of the underlying source code and implementation.

## 1. GitHub Repository (`protbindscreen-wdl`) – Public

- Contains only WDL workflow and task definitions required for workflow registration and execution.
- Public availability is required for synchronization with Dockstore and integration with Terra.
- The repository does not contain the ProtBindScreen source code, implementation modules, or proprietary software components.

## 2. Dockstore – Workflow Registry

- Serves as a workflow registry and versioning layer between GitHub and Terra.
- Stores workflow descriptors and metadata but does not contain the underlying ProtBindScreen implementation.
- References WDL workflows hosted in the public GitHub repository.

## 3. Terra Workspaces – Private

- Workflows are executed within private Terra workspaces accessible only to authorized users.
- Input data, outputs, execution history, and workspace resources remain restricted to workspace members.

## 4. Docker Images – Private Google Artifact Registry

- The ProtBindScreen implementation is distributed through Docker images rather than through the public WDL repository.
- Docker images are stored in a private Google Artifact Registry managed by the Salic Lab.
- Access to the images is controlled through Google Cloud IAM permissions.
- Terra execution environments are granted permission to pull the images, while the underlying source code remains inaccessible to unauthorized users.

## Future Release

The privacy measures described above are intended only for the development, testing, reviewer-access, and pre-release phases of the project. Upon official public release of ProtBindScreen, the source code repository, Docker images, workflow definitions, and associated software resources are expected to become publicly available under the project's selected open-source distribution model.

# Review of currently used docker images

Tasks: 

1. prepare_fasta_files.wdl: `us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/protbindscreen:0.0.7`
2. predict_with_colabfold.wdl: `ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2`
3. predict_with_colabfold_test.wdl: `us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/custom_build_cudabase_mmseqs2bin_colabfold:0.0.10`
4. local_msa_colabfold_search_test.wdl: `us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/custom_build_cudabase_mmseqs2bin_colabfold:0.0.10`
5. analyze_screen.wdl: `us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/protbindscreen:0.0.7`

These images (except `ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2`) were all kept on AML's Google Artifact Registry (GAR) `us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo` . I think this is private. 

An additional image is being used for the jupyter lab environment `protbindscreen_viewer_terra-jupyter-base114` this image was stored both under a private and a public GAR repository. I think the Terra workspace was pulling from the public one? not sure need to continue verifications. Below are the two locations:

- us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/protbindscreen_viewer_terra-jupyter-base114
- us-central1-docker.pkg.dev/global-axe-475818-q0/public-pbsv-repo/protbindscreen_viewer_terra-jupyter-base114

I also have some images stored on docker both on public and private repositories.

Private Docker repositories:

- anamihaelalupan/protbindscreen
- anamihaelalupan/cuda12.4.1cudnn_mmseqs2master_colabfold.addfad4

Public Docker repositories:

- anamihaelalupan/public_pbsv_docker_repo

On my local machine, I have the following Docker images:

(base) anamihaelalupan@MBC-SALIC-D-05776:~/myprojects/protbindscreen_dev$ docker images
REPOSITORY                                                                                                              TAG                IMAGE ID       CREATED        SIZE
anamihaelalupan/protbindscreen                                                                                          e0496c5.0.0.6      96ad812ccf06   5 months ago   2.45GB
protbindscreen_e0496c5                                                                                                  0.0.6              96ad812ccf06   5 months ago   2.45GB
protbindscreen_e0496c5                                                                                                  0.0.4              1407b3b04c45   5 months ago   2.45GB
anamihaelalupan/protbindscreen                                                                                          e0496c5.0.0.3      c4148671496b   5 months ago   2.45GB
protbindscreen_e0496c5                                                                                                  0.0.3              c4148671496b   5 months ago   2.45GB
anamihaelalupan/protbindscreen                                                                                          v0.1.0             80283c774124   5 months ago   3.63GB
protbindscreen                                                                                                          0.1.0              80283c774124   5 months ago   3.63GB
anamihaelalupan/protbindscreen                                                                                          e0496c5.0.0.1      d5ddc6df935b   5 months ago   2.67GB
anamihaelalupan/protbindscreen_e0496c5                                                                                  0.0.1              d5ddc6df935b   5 months ago   2.67GB
protbindscreen_e0496c5                                                                                                  0.0.1              d5ddc6df935b   5 months ago   2.67GB
protbindscreen_e0496c5                                                                                                  0.0.2              d5ddc6df935b   5 months ago   2.67GB
anamihaelalupan/protbindscreen                                                                                          e0496c5            30fb5dc4b5f4   5 months ago   2.44GB
protbindscreen                                                                                                          e0496c5            30fb5dc4b5f4   5 months ago   2.44GB
anamihaelalupan/cuda12.4.1cudnn_mmseqs2master_colabfold.addfad4                                                         0.0.10             2e402177b53c   6 months ago   13.5GB
local/custom_build_cudabase_mmseqs2bin_colabfold                                                                        0.0.10             2e402177b53c   6 months ago   13.5GB
s-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/custom_build_cudabase_mmseqs2bin_colabfold    0.0.10             2e402177b53c   6 months ago   13.5GB
us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/custom_build_cudabase_mmseqs2bin_colabfold   0.0.10             2e402177b53c   6 months ago   13.5GB
local/custom_build_cudabase_mmseqs2bin_colabfold                                                                        0.0.9              e89405b1c847   6 months ago   5.68GB
local/custom_build_cudabase_mmseqs2bin_colabfold                                                                        0.0.8              3e74706eec48   6 months ago   8.73GB
local/custom_build_cudabase_mmseqs2bin_colabfold                                                                        0.0.7              bcec99dce90c   6 months ago   13.1GB
us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/custom_build_cudabase_mmseqs2bin_colabfold   0.0.7              bcec99dce90c   6 months ago   13.1GB
local/custom_build_cudabase_mmseqs2bin_colabfold                                                                        0.0.6              e8fae974a5db   6 months ago   12.2GB
local/custom_build_cudabase_mmseqs2bin_colabfold                                                                        0.0.5              9ad30bed05c2   7 months ago   12.2GB
us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/custom_build_cudabase_mmseqs2bin_colabfold   0.0.5              9ad30bed05c2   7 months ago   12.2GB
local/custom_build_cudabase_mmseqs2bin_colabfold                                                                        0.0.4              5f0a0d8277fc   7 months ago   12.2GB
local/custom_build_cudabase_mmseqs2bin_colabfold                                                                        0.0.3              edb9210466a2   7 months ago   12.1GB
us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/custom_build_cudabase_mmseqs2bin_colabfold   0.0.2              bd579e0c7503   7 months ago   12.1GB
us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/custom_build_cudabase_mmseqs2bin_colabfold   0.0.1              23564813459e   7 months ago   12.1GB
us-central1-docker.pkg.dev/global-axe-475818-q0/protbindscreen-docker-repo/protbindscreen                               0.0.7              c10d037a17c6   7 months ago   3.55GB
anamihaelalupan/public_pbsv_docker_repo                                                                                 0.0.6              d79f7c535656   7 months ago   27.3GB
us-central1-docker.pkg.dev/global-axe-475818-q0/public-pbsv-repo/protbindscreen_viewer_terra-jupyter-base114            0.0.4              6a0e701eed4a   7 months ago   27.3GB
anamihaelalupan/protbindscreen                                                                                          v0.0.2             9316dcaa1bca   7 months ago   4.32GB
ubuntu                                                                                                                  <none>             392fa14dddd0   8 months ago   77.9MB
hello-world                                                                                                             latest             1b44b5a3e06a   9 months ago   10.1kB
ghcr.io/sokrypton/colabfold                                                                                             1.5.5-cuda12.2.2   cba799638e84   2 years ago    6.89GB

# Review of workflows

version: 1.2

workflows:
  - name: Experiment_Table
    subclass: WDL
    primaryDescriptorPath: /tablet_experiment.wdl

      - this is obsolete - deleted
    
  - name: Analyze_Table
    subclass: WDL
    primaryDescriptorPath: /WDL/Workflows/wdl_analyze_mode.wdl
  
      ```wdl
      import "../Tasks/analyze_screen/analyze_screen.wdl" as analyze_screen_t
      ```
    
  - name: Submit_Data
    subclass: WDL
    primaryDescriptorPath: /WDL/Workflows/wdl_submit_mode.wdl
    
      ```wdl
      import "../Tasks/prepare_fasta_files/prepare_fasta_files.wdl" as prepare_fasta_files_t
      import "../Tasks/predict_with_colabfold/predict_with_colabfold.wdl" as predict_with_colabfold_t
      ```
  - name: Submit_Data_Test
    subclass: WDL
    primaryDescriptorPath: /WDL/Workflows/wdl_submit_mode_test.wdl

    ```wdl
    import "../Tasks/prepare_fasta_files/prepare_fasta_files.wdl" as prepare_fasta_files_t
    import "../Tasks/local_msa_colabfold_search/local_msa_colabfold_search_test.wdl" as local_msa_colabfold_search_t
    import "../Tasks/predict_with_colabfold/predict_with_colabfold_test.wdl" as predict_with_colabfold_t
    ```

So, I think that current `Analyze_Table` and `Submit_Data` should be the one that were functional. and that `Submit_Data_Test` was the one I was getting failures before I paused dev.

Renaming of the workflows and tasks to better reflect the main protbindscreen implementation. 

Workflows: 

`SCREEN_SUBMISSION_MODE.wdl` - will be the current `Submit_Data` plus any improvements.

`SCREEN_SUBMISSION_MODE_LOCAL_MSA.wdl` - will be the current `Submit_Data_Test` plus any improvements. This is highly experimental and we do not expect to be ready soon.

`SCREEN_ANALYSIS_MODE.wdl` - will be the current `Analyze_Table` plus any improvements.

Tasks, will keep the same names, which are fine as-are: 

1. prepare_fasta_files.wdl - called by `SCREEN_SUBMISSION_MODE.wdl` and `SCREEN_SUBMISSION_MODE_LOCAL_MSA.wdl`
2. predict_with_colabfold.wdl - called by `SCREEN_SUBMISSION_MODE.wdl`
3. analyze_screen.wdl - called by `SCREEN_ANALYSIS_MODE.wdl`
4. predict_with_colabfold_test.wdl - called by `SCREEN_SUBMISSION_MODE_LOCAL_MSA.wdl`
5. local_msa_colabfold_search_test.wdl - called by `SCREEN_SUBMISSION_MODE_LOCAL_MSA.wdl`

The docker images will have other addreses that will reflect salic.lab.group migration but also updated versions of the probindscreen code (to reflect as closely as possible the current state of the mother protbindscreen repository)

I'll get back to each of these remember what were they doing, reconstruct a minimal working set of workfloes for submission and analysis.
The development with GPU-accelerated local MSA remnains a desire of development but will not be attempted before submission.

# Location of Docker images

The docker images will be located in a separate repository called `protbindscreen_docker_images` that will stay private. 

    


