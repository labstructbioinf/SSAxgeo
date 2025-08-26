# SSAxgeo

SSAxgeo assigns protein secondary structure elements using differential geometry and knot theory. It provides containerised command-line tools for analysing Protein Data Bank entries.

## Features

- Differential geometry–based descriptors.
- Command-line utilities for sampling PDB entries and computing descriptors.
- Singularity container for reproducible execution.
- Workflows to reproduce the analyses from the associated publication.

## Prerequisites

- Git
- Python ≥3.8
- Singularity 3.x
- localpdb (for reproducing the paper analyses)

```bash
python --version
```

## Installation

```bash
git clone --recurse-submodules https://github.com/labstructbioinf/SSAxgeo.git
cd SSAxgeo
sudo singularity build ssaxgeo.sif SingularityFile
```

## Basic Usage

```bash
singularity exec ssaxgeo.sif ssaxgeo /path/to/structure.pdb
```

## Reproduce Paper Analyses

1. **Prepare a local copy of the PDB**

   ```bash
   localpdb_setup -db_path /path/to/mypdb/ -plugins DSSP PDBClustering PDBChain --fetch_cif --fetch_pdb
   ```

   Downloads PDB files and sets up required plugins at `/path/to/mypdb/`.

2. **Sample a clustered PDB**

   ```bash
   ssaxgeo_getSampleOfClstrPDB /path/to/mypdb/ -out_dir /path/to/mydir/ -redundancy 30 -res_lim 2.0 -ncpus 4 -seed 0
   ```

   Generates a clustered set of structures with the specified redundancy in `/path/to/mydir/`.

3. **Compute differential geometry descriptors**

   ```bash
   ssaxgeo_computePDBxgeo --mylocalpdb_path /path/to/mypdb/ --sampled_clstrd_path /path/to/sampled_clust-30.csv --xgeo_output_dir /path/to/mypdb/xgeo_chains/ --ncpus 8 --out_csv /path/to/sampled_clust-30_updated.csv
   ```

   Produces xgeo descriptor files and an updated sample table.

4. **Cluster residues and generate fragments**

   ```bash
   ssaxgeo_clusterResidues /path/to/sampled_clust-30_updated.csv clust-30 -ncpus 8
   ```

   Normalises xgeo values and outputs residue clusters representing structural fragments.

5. **Select canonical regions**

   Use `notebooks/SetCanonicalRegions.ipynb` to filter fragments for geometrical helices and derive canonical sets.

## Roadmap

### v1.0

- [x] bring diffgeo to be part of ssaxgeo
- [x] migrate code for canonical regions detection to rely on localpdb
- [x] add/update and adapt scripts to reproduce paper results more easily (via CLI)
- [ ] adapt old code scripts to rely on new structure
  - [x] ssaxgeo
  - [x] computePDBxgeo
  - [x] getSampleOfClstrPDB
  - [x] prepare_BCX_df
  - [x] getCanonicalClusters
    - [x] convert getCanonical to a notebook
- [x] improve sanitize pdb function to comply with melodia
- [ ] cosmetics changes
- [ ] update documentation
- [ ] test end-to-end

### v1.1

- [ ] add citation
- [ ] add pymol viz support
- [ ] add xgeo Dlang code support

## Citation

If you use SSAxgeo in your research, please cite the associated publication.

## Contributing

Contributions are welcome. Please open an issue to discuss significant changes before submitting a pull request.

## License

This project is licensed under the MIT License.
