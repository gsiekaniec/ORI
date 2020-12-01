## Build package from github

```bash
git clone https://github.com/gsiekaniec/ORI
conda env create -p ./ori-env --file ./ORI/conda/environment.yml
conda activate ./ori-env
conda install conda-build
conda-build ./ORI/conda/ori
```

## Install

```bash
ORI_PACKAGE_PATH=$(conda-build ./ORI/conda/ori --output)
conda install ${ORI_PACKAGE_PATH}
```