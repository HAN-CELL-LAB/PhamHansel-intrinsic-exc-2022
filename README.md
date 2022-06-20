# Intrinsic threshold plasticity: cholinergic activation and role in the neuronal recognition of incomplete input patterns

- Authors: [Tuan Pham](https://github.com/tuanpham96), [Christian Hansel](http://www.hansellab-uchicago.com/)
- Short description: This is the code repository for running the simulations and analyses to study the effects of intrinsic excitability changes on recognition and discrimination in a simple feed-forward network.
- The code repository can be cloned with its submodules by:

  ``` bash
  # the --depth 1 is optional 
  git clone --recurse-submodules --depth 1 \
    https://github.com/HAN-CELL-LAB/PhamHansel-intrinsic-exc-2022
  ```

- The complete data, figures and dependencies (`matlab-ext`) are archived at `figshare`. 

## Repository organization

- `manuscript`: contains `Latex` project to generate manuscript
- `matlab-ext`: contains external submodules.
  - See `matlab-requirements.txt` for specific packages
  - Don't need to run `install-requirements.sh` with if this repository is cloned using `git clone --recurse-submodules`
- Projects:
  - `01-discrim-sim`: discrimination simulation
  - `02-recog-sim`: recognition simulation
  - `03-bioldata-analy`: analyses of biological data
- In side each of the project and manuscript folder, the subfolders `data` and/or `figures` are empty when cloned. To download, see the `figshare` link above. 

## Data

- The follow is run to create a compressed archive to upload onto `figshare`:

  ``` bash
  cd .. # run from parent directory of the local folder of this repository
  tar czfv PhamHansel-intrinsic-exc-2022.tar.xz \
    -X PhamHansel-intrinsic-exc-2022/tar_exclude.txt \
    PhamHansel-intrinsic-exc-2022/
  ```

  - Excludes all the files from `tar_exclude.txt`, so as to only include necessary scripts, dependencies, figures, documents and data.
  - The `matlab-ext` subfolders are also included in case the repositories are no longer available.
- The `tar.xz` is uploaded to `figshare` (see above).

## Figure sources

- **Figure 1**: see [`01-discrim-sim`](01-discrim-sim/README.md). Note: `1a,b` and parts of `1c` were illustrated in `inkscape`.
- **Figure 2, 3, 4**: see [`02-recog-sim`](02-recog-sim/README.md). Note: upper parts of `2a,b`, parts of `3a` and `4c` were illustrated in `inkscape`.
- **Figure 5**: see [`03-bioldata-anly`](03-bioldata-anly/README.md).

