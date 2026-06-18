# Contributing Guide & Workflow

To keep our code clean and prevent the repository from breaking, I have set up some automated "shields."

Please follow this guide strictly.

## Rules
1. **NEVER push directly to `main`.**
    * The `main` branch is locked. If you try to push to it, GitHub will reject you.
2.  **ALWAYS create a new branch.**
    * Create a branch for every lab or feature (e.g., `feature/my_feature`, `fix/boundary-conditions`).
3.  **NO BINARIES or MESHES.**
    * Do not upload `.exe`, `lab-01` (executables), `.o` files, or large mesh files (`.msh`, `.vtk`).

---

## How to Start Developing

### 1. The "Common" Folder Situation
**Read this carefully:** This repository is a fork of the specific university setup.
* **On the University Cluster:** A folder named `common` exists one level up (`../common`).
* **On GitHub/Your PC:** This folder **does not exist**.
* **The Fix:** We have commented out the dependency in `CMakeLists.txt` so the project compiles standalone. **Do not uncomment these lines** unless you are running on the cluster and have the files locally.

### 2. Workflow Steps
1.  **Get the latest code:**
    ```bash
    git checkout main
    git pull origin main
    ```

2.  **Create your branch:**
    ```bash
    git checkout -b feature/my-new-code
    ```

3.  **Write your code:**
    * Place your source files in `src/`.
    * Update `CMakeLists.txt` if you added new `.cpp` files.

4.  **Test Locally (Compilation):**
    If you are on the cluster, load the modules first:
    ```bash
    module load gcc-glibc dealii
    ```
    Then build:
    ```bash
    mkdir build
    cd build
    cmake ..
    make
    ```

5.  **Push your work:**
    ```bash
    git add .
    git commit -m "Implemented the finite element assembly"
    git push origin feature/my-new-code
    ```

---

## Automated "Shields" (CI/CD)

When you open a Pull Request, an **Automated Sanity Check** (`cluster-simulation`) will run immediately.

### What is it doing?
It spins up a Docker container that mimics the university cluster (running **Deal.II v9.6.0**) and tries to compile your code from scratch.

### How to read the results:
* ✅ **Green Check:** Your code compiled successfully. You are safe to merge.
* ❌ **Red Cross:** You broke the build!
    1.  Click **"Details"** on the error.
    2.  Read the compiler error log.
    3.  Fix the code on your computer, commit, and push again. The check will re-run automatically.

---

## Review Process
1.  Once the shield turns **Green**, request a review from the repository owner.
2.  The owner will review the logic.
3.  Once approved, the owner will **Squash and Merge** your code into `main`.
4.  Delete your branch after merging.