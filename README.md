# SSAxgeo

This software performs protein Secondary Structure Assignment using differential geometry and knot theory descriptors.
---

## How to install?

1. Clone the repository:
```{bash}
git clone --recurse-submodules https://github.com/labstructbioinf/SSAxgeo.git
cd SSAxgeo
```

2. Create and activate conda environment:
```{bash}
conda create -n ssaxgeo python=3.10 "pyarrow>=16,<21" ldc dub -c conda-forge -y
conda activate ssaxgeo
```

3. Install SSAxgeo:
```{bash}
# Be sure to use pip from the activated conda environment!
pip install -e .
```

4. Compile diffgeo:
```{bash}
cd diffgeo
dub build --build=release --compiler=ldc2
cd ..
```

### Test the installation

```{bash}
cd test
ssaxgeo -can_dir ../canonical 3tsi_chain_A.pdb
# For more options see
ssaxgeo -h
```

