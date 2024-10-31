# GMT Integrated Modeling for Seismic Analyses

This repository provides the code to simulate the effect of RLE seismic events on the GMT. The simulations employ the dos-actors framework based on the Rust programming language.

![](data/model.dot.png)

The RLE ground accelerations available from the [1K RLE SSSHA release](https://github.com/GMTO/seismic/releases) must be copied to the [gmt-im-seismic/data/](https://github.com/GMTO-Integrated-Modeling/gmt-im-seismic/tree/main/data) folder to run the simulations. Those files have ground acceleration time series upsampled to 1KHz (see CR-05354 for details on the upsampling method). Simulations also require a zip file containing the modal model and the static solution of the telescope FEM. Structural model files are available on the AWS GMT-IM s3 drive. The environmental variable `FEM_REPO` specifies the location of the zip file.

The simulation results are recorded in a parquet file in the data folder. The results include displacements of the pier nodes, mount main axes (AZ, EL, and GIR), and M1 rigid-body motions and hardpoint displacements.

A Python notebook (`processing.ipynb`) can be used to visualize the data and perform preliminary evaluation. The MatLab script (`assess_sssha_dt.m`) produces a collection of plots to analyze the simulation results thoroughly. The script (`rle_latKs_sensitivity.m`) uses simulation results of different telescope structural models to compare the spectral acceleration (SA) of the bottom of the pier. The following plot compares the SA of the node at the bottom of the pier resulting from models with different seismic isolator lateral stiffness.

![](pierSHA_latKs_sens.png)

The considered models and the corresponding lateral stiffness value are reported in the table below. 

| FEM ID | SIS lateral stiffness [N/m] |
|:---:|:------------------:|
| 20240408_1535 | 0.96e9 |
| 20241003_1800 | 1.36e9 |
| 20230817_1808 | 3.00e9 |
| 20241021_1535 | 9.00e9 |



