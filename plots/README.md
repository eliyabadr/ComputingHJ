# Plots: Python vs C++ for the paper runs

Each figure folder holds two single runs side by side: `python/` (the original cluster job) and `cpp/` (the C++ port). Images are the original files, renamed only. No run is completed with images from another run.

`phase_NN.png` is the value function after refinement phase NN (phase 00 = initial value iteration). `iteration_0001_refinement_00.png` is an extra Python plot saved after the first VI sweep.

| Folder | Paper | Python job | Python phases | C++ phases | Same parameters | Python vs C++ |
|---|---|---|---|---|---|---|
| `fig2_dubins_RA_nodiscount` | Fig. 2 | 359053 | 0–13 (job ended in phase 14) | 0–14 | yes | phases 0–13: identical cell counts, values within 1.3e-15 |
| `fig3_evasion_RA_nodiscount` | Fig. 3 | 359054 | 0–9 (job ended in phase 10) | 0–10 | yes | phases 0–9: identical cell counts, values within 1.8e-15 |
| `fig4a_dubins_RA_nodiscount_120` | Fig. 4(a) | 356169 | 0–1 (cancelled in phase 2) | 0–3 | no (ε, δ_max, VI cap) | phase 0: identical, within 8.9e-16 |
| `fig4b_dubins_RA_discount_120` | Fig. 4(b) | 357815 | 0–2 (stopped in phase 3) | 0–6 | no (ε, VI cap) | phases 0–2: identical, within 8.9e-16 |
| `fig5_dubins_avoid_nodiscount_80` | Fig. 5 | 356456 | 0–6 (time limit in phase 7) | 0–6 | no (γ, ε, δ_max) | phases 0–6: identical, within 8.9e-16 |

"Identical" = leaf, boundary, classified-cell and VI-sweep counts per phase, from the run logs.

Plot layouts: for Fig. 2 and Fig. 3 both implementations draw the same value-function figure (C++ rendered from SVG). For Fig. 4(a), 4(b) and 5 the C++ run saved theta-slice rasters instead, so the layout differs from the Python figure.

## Original file names

### fig2_dubins_RA_nodiscount
```
cpp/phase_00.png <- value_function_phase_0_complete.png
cpp/phase_01.png <- value_function_phase_1_complete.png
cpp/phase_02.png <- value_function_phase_2_complete.png
cpp/phase_03.png <- value_function_phase_3_complete.png
cpp/phase_04.png <- value_function_phase_4_complete.png
cpp/phase_05.png <- value_function_phase_5_complete.png
cpp/phase_06.png <- value_function_phase_6_complete.png
cpp/phase_07.png <- value_function_phase_7_complete.png
cpp/phase_08.png <- value_function_phase_8_complete.png
cpp/phase_09.png <- value_function_phase_9_complete.png
cpp/phase_10.png <- value_function_phase_10_complete.png
cpp/phase_11.png <- value_function_phase_11_complete.png
cpp/phase_12.png <- value_function_phase_12_complete.png
cpp/phase_13.png <- value_function_phase_13_complete.png
cpp/phase_14.png <- value_function_phase_14_complete.png
python/iteration_0001_refinement_00.png <- iteration_0001_refinement_00.png
python/phase_00.png <- value_function_phase_0_complete.png
python/phase_01.png <- value_function_phase_1_complete.png
python/phase_02.png <- value_function_phase_2_complete.png
python/phase_03.png <- value_function_phase_3_complete.png
python/phase_04.png <- value_function_phase_4_complete.png
python/phase_05.png <- value_function_phase_5_complete.png
python/phase_06.png <- value_function_phase_6_complete.png
python/phase_07.png <- value_function_phase_7_complete.png
python/phase_08.png <- value_function_phase_8_complete.png
python/phase_09.png <- value_function_phase_9_complete.png
python/phase_10.png <- value_function_phase_10_complete.png
python/phase_11.png <- value_function_phase_11_complete.png
python/phase_12.png <- value_function_phase_12_complete.png
python/phase_13.png <- value_function_phase_13_complete.png
```

### fig3_evasion_RA_nodiscount
```
cpp/phase_00.png <- value_function_phase_0_complete.png
cpp/phase_01.png <- value_function_phase_1_complete.png
cpp/phase_02.png <- value_function_phase_2_complete.png
cpp/phase_03.png <- value_function_phase_3_complete.png
cpp/phase_04.png <- value_function_phase_4_complete.png
cpp/phase_05.png <- value_function_phase_5_complete.png
cpp/phase_06.png <- value_function_phase_6_complete.png
cpp/phase_07.png <- value_function_phase_7_complete.png
cpp/phase_08.png <- value_function_phase_8_complete.png
cpp/phase_09.png <- value_function_phase_9_complete.png
cpp/phase_10.png <- value_function_phase_10_complete.png
python/iteration_0001_refinement_00.png <- iteration_0001_refinement_00.png
python/phase_00.png <- value_function_phase_0_complete.png
python/phase_01.png <- value_function_phase_1_complete.png
python/phase_02.png <- value_function_phase_2_complete.png
python/phase_03.png <- value_function_phase_3_complete.png
python/phase_04.png <- value_function_phase_4_complete.png
python/phase_05.png <- value_function_phase_5_complete.png
python/phase_06.png <- value_function_phase_6_complete.png
python/phase_07.png <- value_function_phase_7_complete.png
python/phase_08.png <- value_function_phase_8_complete.png
python/phase_09.png <- value_function_phase_9_complete.png
```

### fig4a_dubins_RA_nodiscount_120
```
cpp/phase_00.png <- slices_phase_0.png
cpp/phase_01.png <- slices_phase_1.png
cpp/phase_02.png <- slices_phase_2.png
cpp/phase_03.png <- slices_phase_3.png
python/iteration_0001_refinement_00.png <- iteration_0001_refinement_00.png
python/phase_00.png <- value_function_phase_0_complete.png
python/phase_01.png <- value_function_phase_1_complete.png
```

### fig4b_dubins_RA_discount_120
```
cpp/phase_00.png <- slices_phase_0.png
cpp/phase_01.png <- slices_phase_1.png
cpp/phase_02.png <- slices_phase_2.png
cpp/phase_03.png <- slices_phase_3.png
cpp/phase_04.png <- slices_phase_4.png
cpp/phase_05.png <- slices_phase_5.png
cpp/phase_06.png <- slices_phase_6.png
python/iteration_0001_refinement_00.png <- iteration_0001_refinement_00.png
python/phase_00.png <- value_function_phase_0_complete.png
python/phase_01.png <- value_function_phase_1_complete.png
python/phase_02.png <- value_function_phase_2_complete.png
```

### fig5_dubins_avoid_nodiscount_80
```
cpp/phase_00.png <- slices_phase_0.png
cpp/phase_01.png <- slices_phase_1.png
cpp/phase_02.png <- slices_phase_2.png
cpp/phase_03.png <- slices_phase_3.png
cpp/phase_04.png <- slices_phase_4.png
cpp/phase_05.png <- slices_phase_5.png
cpp/phase_06.png <- slices_phase_6.png
python/iteration_0001_refinement_00.png <- iteration_0001_refinement_00.png
python/phase_00.png <- value_function_phase_0_complete.png
python/phase_01.png <- value_function_phase_1_complete.png
python/phase_02.png <- value_function_phase_2_complete.png
python/phase_03.png <- value_function_phase_3_complete.png
python/phase_04.png <- value_function_phase_4_complete.png
python/phase_05.png <- value_function_phase_5_complete.png
python/phase_06.png <- value_function_phase_6_complete.png
```
