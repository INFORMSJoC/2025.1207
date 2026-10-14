[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Statistical Inference in Conditional Value-at-Risk Optimization

MATLAB programs and results accompanying the paper by Qixin Wang, Hao Cao,
and Jian-Qiang Hu (JOC-2025-03-OA-1207.R3). See [AUTHORS](AUTHORS)
and the [MIT license](LICENSE).

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

`src/` contains 12 experiment scripts and 7 correlated-noise helper functions.
The portfolio and lognormal scripts use a consistent `1/sqrt(2)` CI scale.
`results/tables/` contains the manuscript's 12 TeX table fragments.

| Paper item | MATLAB file in `src/` | TeX file in `results/tables/` |
| --- | --- | --- |
| Figure 1, Tables 1–2: portfolio | `CVaR_SA_duowei_jinrong.m` | `portfolio.tex`, `Coverage_portfolio.tex` |
| Section 7.1: Hu comparison | `CVaR_SA_duowei_jinrong_jiaqiao_duibi.m` | — |
| Table 3: SAA | `CVaR_SA_duowei_jinrong_SAA_duibi.m` | `SAA.tex` |
| Figure EC.1, Table EC.1: multimodal | `CVaR_SA_nonconvex.m` | `nonconvex.tex` |
| Figures EC.2–3, Tables EC.2–3: correlated noise | `CVaR_SA_duowei_ex5_correlated_noise.m` | `k_0.tex`, `Correlated_noise.tex` |
| Figure EC.4, Tables EC.4/EC.6: lognormal | `CVaR_SA_tail_lognormal.m` | `lognormal1.tex`, `Coverage_lognormal.tex` |
| Figure EC.5, Tables EC.5/EC.7: regenerative queue | `CVaR_SA_queue_correlated_case.m` | `queue_regen.tex`, `Coverage_queueing_regen.tex` |
| Figure EC.6, Table EC.8: original queue | `CVaR_SA_duowei_queue.m` | `lognormal2.tex` |
| Table EC.9: batching comparison | `CVaR_SA_tail_lognormal_batch_duibi.m` | `Coverage_batch.tex` |

The retained quadratic-model scripts are `CVaR_SA_duowei_ex5.m`,
`CVaR_SA_duowei_ex5_jiaqiao_duibi.m`, and `CVaR_SA_huatuzhuanyong.m`.

The original plotting notebook is [CVaR_SA_huatu.ipynb](scripts/CVaR_SA_huatu.ipynb).
Its numeric arrays are also available in [convergence_rate_data.csv](results/rates/convergence_rate_data.csv).
The active theta series in cells 1–2 corresponds to Figure EC.6(b).
The CSV identifies each list by `cell`, `source_line` and `array_id`;
`variable` and `state` distinguish assigned series from saved alternatives.
`input_value` is the unlogged input, `stored_scale` records the source transformation,
and `paired_n` is populated for the assigned RMSE series in cell 1.

## Use

MATLAB requires Statistics and Machine Learning Toolbox; scripts using symbolic
references also require Symbolic Math Toolbox. Copy the desired script to a working
directory, select the experiment's budgets, repetition count, distribution and
commented output blocks, then run that copy. For correlated noise, copy
the entrypoint and its seven helper files together from `src/` to a working directory,
configure `correlated_paper_settings.m`, and run the entrypoint. The helper files are
`correlated_paper_settings.m`, `run_correlated_paper_tables.m`, `simulate_correlated_case.m`,
`finalize_correlated_tables.m`, `repair_correlated_rmse.m`, `correlated_cvar_reference.m`,
and `export_correlated_matlab_excel.m`. Parallel Computing Toolbox is optional.

For the active rate plot, open the notebook with Python/Jupyter and run cells 1–2
using matplotlib, SciPy and seaborn. To plot another saved series, assign its array
and matching iteration budgets in a working copy. Rates are fitted using natural
logarithms of the budgets and RMSEs.

The TeX files contain reported table values and use the manuscript's LaTeX preamble
and table macros. View them as text or include them in the manuscript to typeset them.

Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu. All files are licensed under MIT; see [LICENSE](LICENSE).
The CSV license notice is in [convergence_rate_data.csv.license](results/rates/convergence_rate_data.csv.license).
