# Biological data analysis

## Disclaimer

Please excuse some of the ways I organized the data folder and how I named them. Things changed over time and decisions were made that I did not have the time to re-organize and re-run to make sure everything is tip-top and clean. Just email or raise an issue if there needs to be any clarification.

## V1 data

This was obtained from communication with the authors of the following paper, hence no data will be committed to the `data/v1-data` folder of the data repository. Please contact the authors of the following paper directly for the data.

> Li B, Routh BN, Johnston D, Seidemann E, Priebe NJ. Voltage-Gated intrinsic conductances shape the Input-Output relationship of cortical neurons in behaving primate V1, Neuron, 107(1), (2020), 185–196.e4.

## S1 data

These data were obtained in the Hansel lab, from the following paper by Gill & Hansel 2020, as well as additional recordings (prefixes TP: Tuan Pham, SEB: Silas Busch, CH: Christian Hansel). All data were committed in the data repository (refer to the main `README.md` file) but only selected data were analyzed due to the quality of recordings and cell healths.

> Gill, Daniel F., and Christian Hansel. 2020. Muscarinic Modulation of SK2-Type K+ Channels Promotes Intrinsic Plasticity in L2/3 Pyramidal Neurons of the Mouse Primary Somatosensory Cortex. eNeuro 7 (2).

The IP data are in `data/gill-data` (for the previous paper) and `data/supp-data` (additional recordings of oxo-m and oxo-m + somatic depolarization conditions). The 4 recordings to estimate the threshold upper bounds (using "refractory thresholds" as mentioned in the manuscript) are in `data/extra-data` (combined from `data/{pham-data,hansel-data}`).

## Figures sources

| Description  | Script | Input | Output | Manuscript | Notes |
| --- | --- | --- | --- | --- | --- |
| V1 visualization | `A2_visualize_v1_data.m` | `data/v1-data/v1-invivo_Li2020.mat` (not available in repo) | `figures/v1-data/v1-violin.pdf` | `Fig 5a` | |
| S1 visualization | `E2_visualize_s1_ip_data.m` | `data/ip-data` | `figures/ip-data` <br/> `s1-violin.pdf` <br/> `s1-change.pdf` <br/> `s1-chi2-pie.pdf` | `Fig 5` <br/> `5b` <br/> `5d` <br/> `5e` | |
| Recognition simulation with S1 range | `E4_visualize_model_and_biol_data.m` | `data/model-data/recog_data.mat` <br/> `data/ip-data/ip-select-data-analysis.mat` <br/> `data/extra-s1/processed/hansel-pham-extra-s1-vthres.mat` | `figures/model-data/model-with-biol-range.pdf` | `Fig 5c` | |
