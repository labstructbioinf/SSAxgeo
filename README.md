# SSAxgeo

This software performs protein Secondary Structure Assignment using differential geometry and knot theory descriptors.
---

## How to install?

1. Clone the repository:
```{bash}
git clone --recurse-submodules https://github.com/labstructbioinf/SSAxgeo.git
cd SSAxgeo
```

2. Create conda environment:
```{bash}
conda create -n ssaxgeo python=3.10 -y
conda activate ssaxgeo
conda install -c conda-forge "pyarrow>=16,<21" ldc dub -y
```

3. Install SSAxgeo
```{bash}
pip install -e .
```

4. Compile diffgeo
```{bash}
cd diffgeo
dub build --build=release --compiler=ldc2
cd ..
```

### Test the installation

```{bash}
singularity exec ssaxgeo.sif ssaxgeo [pdb_filepath]
```

