## Install from anaconda cloud

```bash
conda create -p ori_env
conda activate ori_env
conda install -c gsiekaniec -c conda-forge ori
```

## Build package from github

```bash
git clone https://github.com/gsiekaniec/ORI
conda create -p ./ori-env 
conda activate ./ori-env
conda install conda-build -c conda-forge
conda-build ./ORI/conda/ori -c conda-forge
```

## Install

```bash
ORI_PACKAGE_PATH=$(conda-build ./ORI/conda/ori --output)
conda install --update-deps ${ORI_PACKAGE_PATH}
```
