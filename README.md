# HelixContinuum

Here all the methods/scripts/notebooks related to Helix Continuum project are stores.

---

## How to? (quick and dirty)

### Install

SSAxgeo relies on [Melodia](https://github.com/rwmontalvao/Melodia) to compute the differential geometry representation.
After Melodia is available, install SSAxgeo:

```{bash}

git clone https://github.com/labstructbioinf/SSAxgeo.git
cd SSAxgeo/
pip install -e ./
```

### Run ssaxgeo

```{bash}
ssaxgeo [pdb_filepath]
```

Precomputed differential geometry can be provided to ssaxgeo via:

```{bash}
ssaxgeo my.pdb -xgeo_flpath my_pdb_xgeo.csv
```

----
## REPRODUCE PAPER ANALYSES


### 0 - Get a local copy of the PDB

To reproduce the analyses presented on the paper, be sure localpdb is available on your environment.
Then, setup your local pdb copy:

```{bash}
localpdb_setup -db_path /path/to/mypdb/ -plugins DSSP PDBClustering PDBChain --fetch_cif --fetch_pdb  
```
This process most likely will take a long time.

### 1 - Get a sampling of a clustered PDB

Once the local pdb copy is in place, compute a clustered pdb with a given sequence redundancy.
For instance, with the command bellow the user can obtain entries clustered by 30% of redundance and entry with at least 2 angstron resolutions.

```{bash}
getSampleOfClstrPDB /path/to/mypdb/ -out_dir /path/to/mydir/ -redundancy 30 -res_lim 2.0 -ncpus 4 -seed 0 
```

```
usage: getSampleOfClstrPDB [-h] [-redundancy REDUNDANCY] [-out_dir OUT_DIR] [-res_lim RES_LIM] [-ncpus NCPUS] [-seed SEED] mylocalpdb

This script loads data from localpdb, select a given clustered PDB, select randomly one exemplar of each cluster and save results as csv files.

positional arguments:
  mylocalpdb            Path to a local PDB copy (must be obtained by localpdb package)

options:
  -h, --help            show this help message and exit
  -redundancy REDUNDANCY
                        redundancy by sequence identity [100, 95, 90, 70, 50 and 30]
  -out_dir OUT_DIR      Output directory (default=working dir)
  -res_lim RES_LIM      resolution limit of structures to be considered (default=2.0)
  -ncpus NCPUS          number of cpus to use (default = 1)
  -seed SEED            seed for random number generator (default = None
```
### 2 - compute differential geometry descriptors

For each entry on the clustered pdb, we need to compute our differential geometry descriptors:

```{bash}
computePDBxgeo --mylocalpdb_path /path/to/mypdb/ --sampled_clstrd_path /path/to/sampled_clust-30.csv --xgeo_output_dir /path/to/mypdb/xgeo_chains/ --ncpus 8 --out_csv /path/to/myclstrdPDB_updated.csv
```

```
usage: computePDBxgeo [-h] --mylocalpdb_path MYLOCALPDB_PATH --sampled_clstrd_path SAMPLED_CLSTRD_PATH [--xgeo_output_dir XGEO_OUTPUT_DIR] [--ncpus NCPUS]
                      [--out_csv OUT_CSV]

Compute xgeo data for a given set of chain provided.

options:
  -h, --help            show this help message and exit
  --mylocalpdb_path MYLOCALPDB_PATH
                        path to a localpdb database
  --sampled_clstrd_path SAMPLED_CLSTRD_PATH
                        path to a sampled clustered csv (produced by getSampleOfCLstrPDB)
  --xgeo_output_dir XGEO_OUTPUT_DIR
                        path of a dir to store xgeo csv files (default = xgeo_output_dir+"/xgeo_chains/"
  --ncpus NCPUS         Number of cpus to be used (default=1)
  --out_csv OUT_CSV     Description of out_csv
```
### 3 - identify geometrical helices of protein structures



### 4 - Select canonical regions
----
To obtain a dataframe with all the descriptors and run the SS_assignment just run the _prepare_BCX_df.py_
```{bash}
>$ python prepare_BCX_df.py /path/to/lPDB_metadata.csv bc-50 -ncpus 10
```

    usage: prepare_BCX_df.py [-h] [-ncpus NCPUS] [-save_dir SAVE_DIR]
                              [-canonical_dir CANONICAL_DIR]
                              lpdb_csv bc_group_col

    This script load precompute data from a lPDB copy, get a sample from a given
    clustered PDB, compute normalized and smoothed values for each entry, standardize
    coordinates and run the SSEx assignment on all residues.

    positional arguments:
      lpdb_csv              path to a localPDB metadata csv
      bc_group_col          column name with cluster number for each
                            entry (ex. 'bc-90' for BC90)

    optional arguments:
      -h, --help            show this help message and exit
      -ncpus NCPUS          number of cpus to use (default=1)
      -save_dir SAVE_DIR    set directory for output files (default=Working dir)
      -canonical_dir CANONICAL_DIR
                            set directory of canonical regions to be loaded
                            (default=canonical folder at this script directory)


----
### 2 - Run SSAx assignment on a single entry
To run the Secondary Structure Assignment based on Differential Geometry (xgeo) on a given pdb file, you can use the _'RunSSax.py'_ script.


```{bash}
>$ ssaxgeo pdbfile.pdb xgeofile.csv
```

    usage: ssaxgeo [-h] [-can_dir CAN_DIR]
                      [-out_dir OUT_DIR] [-min_dist MIN_DIST]
                      [-pp2_max PP2_MAX] [-prefix PREFIX]
                      pdb_flpath xgeo_flpath

    This script run the Secondary Structure Assignment based on differential
    geometry (xgeo) descriptors on a specified entry and save results as csv
    files.

    positional arguments:
      pdb_flpath          protein chain coordinate file path
      xgeo_flpath         protein chain xgeo file path

    optional arguments:
      -h, --help          show this help message and exit
      -can_dir CAN_DIR    SSE canonical dataframes directory (default=script dir)
      -out_dir OUT_DIR    Output directory (default=working dir)
      -min_dist MIN_DIST  Minimum distance from canonical Pi/Alpha/3(10) to consider
                          as a label (default=0.2)
      -pp2_max PP2_MAX    Maximum distance from canonical PP2 to consider as PP2
                          (default=0.07)
      -prefix PREFIX      Prefix for ouput files (default=None)
---
### 3 - Generate canonical groups for SSE
To compute the canonical regions used for SSAx, one can use the _GetCanonicalClusters.py_.
```{bash}
>$ python GetCanonicalClusters.py /path/to/grp_state_BC_file.p
```

    usage: GetCanonicalClusters.py [-h] [-save_dir SAVE_DIR]
                           [-kmeans_nclstrs KMEANS_NCLSTRS]
                           [-pp2_cluster_n PP2_CLUSTER_N]
                           [-min_cluster_size MIN_CLUSTER_SIZE]
                           [-pp2_max PP2_MAX] [-Pc_MIN PC_MIN]
                           [-alpha_cluster_n ALPHA_CLUSTER_N]
                           [-pi_cluster_n PI_CLUSTER_N]
                           [-three_cluster_n THREE_CLUSTER_N] [-show_plot]
                           [-plot_debug_data]
                           grp_flpath

    This script run the bootstrap solution to identify to the canonical regions
    based on a representative set of protein structures.

    positional arguments:
      grp_flpath            pickle file path for a PDBx.group object

    optional arguments:
      -h, --help            show this help message and exit
      -save_dir SAVE_DIR    directory to save output files (default= working dir)
      -kmeans_nclstrs KMEANS_NCLSTRS
                        The number of the cluster to be used at k-means
                        (default=6).
      -pp2_cluster_n PP2_CLUSTER_N
                        The number of the cluster which corresponds to
                        PP2(default=2). [WARNING: check the plots!]
      -min_cluster_size MIN_CLUSTER_SIZE
                        Minimum number of points that HDBSCAN can should
                        consider for a cluster (default=15)
      -pp2_max PP2_MAX      Maximum distance from canonical PP2 to consider as
                            PP2 (default=0.07)
      -Pc_MIN PC_MIN        Minimum membership probability for a point to be considered
                            as part of the 'core' (default=0.4)

      -alpha_cluster_n ALPHA_CLUSTER_N
                        cluster number assigned as alpha canonical region
      -pi_cluster_n PI_CLUSTER_N
                        cluster number assigned as pi canonical region
      -three_cluster_n THREE_CLUSTER_N
                        cluster number assigned as 3(10) canonical region
      -show_plot            show generated plots at runing (default=False)
      -plot_debug_data      generate plots for intermediate steps for debugging
                            (default=False)

### Algorithm description

---
## TODO: (ver 1.0)

- [x] bring diffgeo to be part of ssaxgeo
- [x] migrate code for canonical regions detection to rely on localpdb
- [x] add/update and adapt scripts to reproduce paper results more easily (via CLI)
- [ ] adapt old code scripts to rely on new structure 
  - [x] ssaxgeo
  - [x] computePDBxgeo
  - [x] getSampleOfClstrPDB
  - [x] prepare_BCX_df
  - [ ] getCanonicalClusters
- [x] improve sanitize pdb function to comply with melodia
- [ ] add ssaxgeo on cli via container (not in path)
- [ ] cosmetics changes
- [ ] update documentation

## TODO: (ver 1.1)
- [ ] add citation
- [ ] add pymol viz support
- [ ] add xgeo Dlang code suport
