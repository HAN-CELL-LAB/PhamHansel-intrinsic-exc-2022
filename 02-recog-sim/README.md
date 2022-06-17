# Recognition toy example and simulation

Before running any scripts, please move the them into the same folder of this `README.md` file, i.e. the scripts should be in `02-recog-sim`, not in `02-recog-sim/scripts`.

| Description  | Script | Input | Output | Manuscript | Notes |
| --- | --- | --- | --- | --- | --- |
| Toy example for recognition | `A_toy_demo.m` | None | `figures/toy-recog` <br/> `example-uniform_combo-at-X_uniform.pdf` <br/> `example-equal_combo-at-X_equal.pdf` <br/> `summary.pdf` | `Fig 2` <br/> `2a` <br/> `2b` <br/> `2c,d` | |
| Components for demo in `Fig 3` | `B_plot_demo.m` | None | `figures/demo` <br/> `words.pdf` <br/> `ip_heaviside.pdf` <br/> `WXY_alphaW-noverlap.pdf` <br/>`inpout_percentcomplete.pdf` | `Fig 3` <br/> Latex words  <br/> part of `3a` <br/> `3b` <br/> `3c` | |
| Recognition model and simulation | `C_run_recog_sim.m` | None | `data/recog_data.mat` | data of `Fig 4, 5c` | this file is also copied to `03-bioldata-anly/data/model-data` |
| Visualization of the model performances | `D_plot_recog_perf.m` | `data/recog_data.mat` | `figures/main-recog` <br/> `recog-perf.pdf` <br/> `optim-dtheta.pdf` | `Fig 4` <br/> `4a` <br/> `4d` | |
| Visualization of ROC curves | `E_plot_recog_roc.m` | `data/recog_data.mat` | `figures/main-recog` <br/> `roc-curves.pdf` <br/> `roc-auc.pdf` | `Fig 4` <br/> `4b` | `roc-auc.pdf` is just supplemental |
| Visualization of ratios of performance rates | `F_recog_ratio_rate.m` | `data/recog_data.mat` | `figures/main-recog` <br/> `supp_recognition_nPRratio-vs-recognition_TPR.pdf` <br/>  `supp_recognition_PLLR-vs-recognition_TPR.pdf` | Not in MS | supplemental |
