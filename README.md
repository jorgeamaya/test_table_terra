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
