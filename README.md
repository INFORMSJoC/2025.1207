[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Statistical Inference in Conditional Value-at-Risk Optimization

MATLAB programs and results accompanying the paper by Qixin Wang, Hao Cao,
and Jian-Qiang Hu (JOC-2025-03-OA-1207). See [AUTHORS](AUTHORS)
and the [MIT license](LICENSE).

All code in this repository was developed by the authors specifically for this
project. No source code from other projects has been reused.

## Cite

Please cite the paper and its accompanying code:
[paper DOI](https://doi.org/10.1287/ijoc.2025.1207) and
[code DOI](https://doi.org/10.1287/ijoc.2025.1207.cd).

```bibtex
@misc{WangCaoHu2026CVaRCode,
  author = {Qixin Wang and Hao Cao and Jian-Qiang Hu},
  title = {{Statistical Inference in Conditional Value-at-Risk Optimization}},
  publisher = {INFORMS Journal on Computing},
  year = {2026},
  doi = {10.1287/ijoc.2025.1207.cd},
  url = {https://github.com/INFORMSJoC/2025.1207},
  note = {Available for download at https://github.com/INFORMSJoC/2025.1207}
}
```

## Files and paper results

`src/` contains 9 experiment scripts and 7 correlated-noise helper functions.
The portfolio and lognormal scripts use a consistent `1/sqrt(2)` CI scale.
`results/tables/` contains the manuscript's 12 TeX table fragments.

| Paper item | MATLAB file in `src/` | TeX file in `results/tables/` |
| --- | --- | --- |
| Figure 1, Tables 1–2: portfolio | `run_portfolio_experiment.m` | `table_1_portfolio_estimation.tex`, `table_2_portfolio_ci_coverage.tex` |
| Section 7.1: Hu comparison | `run_portfolio_hu_comparison.m` | — |
| Table 3: SAA | `run_portfolio_saa_comparison.m` | `table_3_saa_relative_errors.tex` |
| Figure EC.1, Table EC.1: multimodal | `run_multimodal_experiment.m` | `table_ec1_multimodal_estimation.tex` |
| Figures EC.2–3, Tables EC.2–3: correlated noise | `run_correlated_noise_experiments.m` | `table_ec2_correlated_noise_autocorrelations.tex`, `table_ec3_correlated_noise_estimation.tex` |
| Figure EC.4, Tables EC.4/EC.6: lognormal | `run_lognormal_experiment.m` | `table_ec4_lognormal_estimation.tex`, `table_ec6_lognormal_ci_coverage.tex` |
| Figure EC.5, Tables EC.5/EC.7: regenerative queue | `run_regenerative_queue_experiment.m` | `table_ec5_regenerative_queue_estimation.tex`, `table_ec7_regenerative_queue_ci_coverage.tex` |
| Figure EC.6, Table EC.8: original queue | `run_original_queue_experiment.m` | `table_ec8_original_queue_estimation.tex` |
| Table EC.9: batching comparison | `run_lognormal_batching_comparison.m` | `table_ec9_batching_ci_coverage.tex` |

`results/figures/` contains all 14 figure files used in the paper and its
electronic companion:

| Paper figure | Files in `results/figures/` |
| --- | --- |
| Figure 1: portfolio | `figure_1a_portfolio_trajectories.pdf`, `figure_1b_portfolio_convergence_rates.png` |
| Figure EC.1: multimodal | `figure_ec1a_multimodal_trajectories.pdf`, `figure_ec1b_multimodal_convergence_rates.png` |
| Figure EC.2: correlated noise trajectories | `figure_ec2a_correlated_exponential_trajectories.pdf`, `figure_ec2b_correlated_normal_trajectories.pdf` |
| Figure EC.3: correlated noise rates | `figure_ec3a_correlated_exponential_convergence_rates.png`, `figure_ec3b_correlated_normal_convergence_rates.png` |
| Figure EC.4: lognormal | `figure_ec4a_lognormal_trajectories.pdf`, `figure_ec4b_lognormal_convergence_rates.png` |
| Figure EC.5: regenerative queue | `figure_ec5a_regenerative_queue_trajectories.pdf`, `figure_ec5b_regenerative_queue_convergence_rates.png` |
| Figure EC.6: original queue | `figure_ec6a_original_queue_trajectories.pdf`, `figure_ec6b_original_queue_convergence_rates.png` |

The plotting notebook is [plot_convergence_rates.ipynb](<scripts/plot_convergence_rates.ipynb>).
The RMSE series for the notebook's six cases are also available in
[convergence_rate_data.csv](results/rates/convergence_rate_data.csv), with columns
`case`, `paper_figure`, `quantity`, `n` (budget), and `rmse` (unlogged values).
The notebook redraws Figure 1(b), EC.1(b), EC.3(a,b), EC.5(b), and EC.6(b)
from stored RMSE vectors. Figure EC.4(b) is included as the original manuscript image.

## Use

MATLAB requires Statistics and Machine Learning Toolbox; scripts using symbolic
references also require Symbolic Math Toolbox. Copy the desired script to a working
directory, select the experiment's budgets, repetition count, distribution and
commented output blocks, then run that copy. For correlated noise, copy
the entrypoint and its seven helper files together from `src/` to a working directory,
configure `correlated_paper_settings.m`, and run the entrypoint. The helper files are
`correlated_paper_settings.m`, `run_correlated_paper_tables.m`, `simulate_correlated_case.m`,
`finalize_correlated_tables.m`, `compute_correlated_rmse.m`, `correlated_cvar_reference.m`,
and `export_correlated_matlab_excel.m`. Parallel Computing Toolbox is optional.

Open the notebook with Python/Jupyter, using matplotlib, SciPy and seaborn.
Set `CASE` in the first code cell and run both code cells. Available cases are
`portfolio`, `nonconvex`, `correlated_exp`, `correlated_normal`, `queue_regen`,
and `queue`; each selects the corresponding paper rate plot.
The default `queue` case corresponds to Figure EC.6(b). Rates are fitted using
natural logarithms of the budgets and RMSEs. Set `SAVE_FIGURE = True` to save the
selected plot under `results/generated_rates/`.

