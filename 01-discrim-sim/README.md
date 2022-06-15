# Discriminatability in output layer with minimal input activation

Before running any scripts, please move the them into the same folder of this `README.md` file, i.e. the scripts should be in `01-discrim-sim`, not in `01-discrim-sim/scripts`.

| Description  | Script | Input | Output | Manuscript | Notes |
| --- | --- | --- | --- | --- | --- |
| Demonstration of distances | `A0_plot_dist_demo.m` | None | `figures/demo-dist.pdf` | parts of `Fig 1b` | |
| Variations of `N_Y` and `k` with distance optimization | `A1_sim_optim_with_dist.m` | None | `data/varyNYk/sim_results.mat` <br/> `data/varyNYk/summary_analyses.mat` | data source of `Fig 1c,d` | The second file is the mean of the first (hence smaller) |
| Plot mean pairwise distances with minimal activation | `A2_plot_optim_with_dist.m` | `data/varyNYk/summary_analyses.mat` (above) | `figures/dist-plot.pdf` | source of `Fig 1c,d` | |
| Variations of `N_Y` and plot entropies and number of unique in the output layer | `B1_optimnet_dtheta_varyNY.m` | None | `data/vary10XNY/*.mat` <br/> `figures/demo-optim.pdf`  <br/> `figures/discrim-ent-unq-10X-varyNY.pdf` <br/> `figures/discrim-pairs-anly-10X-varyNY.pdf` | Not in MS | Supplemental |
| Variations of `N_Y` and `N_X` (more variations and trials than above) on [NSG](https://www.nsgportal.org/index.html) | `C1_simnsg_optimnet_dtheta_varyNXY.m` | None | `data/varyNXY/*.mat` <br/> `figures/discrim-ent-varyNXY.pdf` <br/> `figures/discrim-unq-varyNXY.pdf` | Not in MS | Supplemental <br/> **NEED** the `../matlab-ext` to be in same folder as script when submitted to NSG |
