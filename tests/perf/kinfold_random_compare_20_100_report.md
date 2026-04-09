# Kinfold Random Benchmark Report

Date: 2026-04-06

## Workload

- Random RNA sequences with seed `0`
- One sequence per length from `20` to `100 nt`
- `81` total sequences
- `1000` Kinfold trajectories per sequence
- Old binary: `/Users/howardsalis/ViennaRNA/old/ViennaRNA/src/Kinfold/Kinfold`
- New binary: `/Users/howardsalis/ViennaRNA/src/Kinfold/Kinfold`
- Compare artifact: `/tmp/kinfold_random_compare_20_100.json`
- CSV artifact: `/tmp/kinfold_random_compare_20_100.csv`

## Summary

| Metric | Old | New | Result |
| --- | ---: | ---: | ---: |
| Total runtime (ns) | `3,916,106,175,790` | `1,065,516,944,788` | `3.675x` speed-up |
| Mean runtime per sequence (ns) | `48,346,989,824.57` | `13,154,530,182.57` | `72.8%` lower runtime |
| Mean folding time across cases | `500.3038106419753` | `500.30382232098765` | effectively identical |
| Most stable energy match rate | - | - | `81/81` |
| Most stable structure match rate | - | - | `81/81` |

The optimized `Kinfold` is substantially faster on this workload while preserving the lowest-energy structure and energy found for every tested sequence.

## Runtime By Length Range

| Length range | Cases | Old total wall ns | New total wall ns | Speed-up |
| --- | ---: | ---: | ---: | ---: |
| `20-39` | `20` | `21,732,399,666` | `18,636,107,666` | `1.166x` |
| `40-59` | `20` | `136,869,928,000` | `66,083,135,457` | `2.071x` |
| `60-79` | `20` | `791,803,588,125` | `242,067,900,125` | `3.271x` |
| `80-100` | `21` | `2,965,700,259,999` | `738,729,801,540` | `4.015x` |

The speed-up increases strongly with sequence length.

## Notable Cases

### Largest Speed-Ups

| Case | Length | Speed-up | Stable energy | Stable structure |
| --- | ---: | ---: | ---: | --- |
| `L098_S000` | `98` | `4.816x` | `-30.3` | `.(((((.(((.(((((..((((.....))))..)))(.((.(((....))).)).)))))))))))...((((....(((.....)))....))))..` |
| `L100_S000` | `100` | `4.606x` | `-33.5` | `((((......(((((((((.((((.((((((((.((((.(((((..(.....)..))))).))))..)))))))).)))).).....)))))))))))).` |
| `L096_S000` | `96` | `4.552x` | `-23.0` | `((...((((..(((.((((...........))))(((.(((((((.((.(((((....))))))).)))).))).))))))..))))...))....` |
| `L090_S000` | `90` | `4.409x` | `-27.9` | `(((.((((.((((((...(((.(((((.((((((....((((...))))...)).)))).)).))).)))..))))))))))..)))...` |
| `L086_S000` | `86` | `4.277x` | `-28.3` | `....(((.(((...(((((.(((((((...(((.((.((((((....)))))).)).)))))).))))....))))).))).))).` |

### Slowest Relative Cases

| Case | Length | Speed-up | Stable energy | Stable structure |
| --- | ---: | ---: | ---: | --- |
| `L021_S000` | `21` | `0.997x` | `-0.8` | `.(((.........))).....` |
| `L025_S000` | `25` | `1.006x` | `-2.6` | `...(((.(((....))).)))....` |
| `L028_S000` | `28` | `1.007x` | `0.0` | `............................` |
| `L022_S000` | `22` | `1.011x` | `-5.6` | `.....((.(((....)))))..` |
| `L023_S000` | `23` | `1.014x` | `-0.6` | `.((((....))))..........` |

Short sequences are roughly flat; the largest gains occur on longer RNAs.

## Representative Lengths

| Length | Speed-up | Stable energy | Stable structure |
| --- | ---: | ---: | --- |
| `20` | `1.579x` | `-3.2` | `...(((((...)))))....` |
| `40` | `1.428x` | `-10.2` | `(((((((.(((((....).)))).....))))))).....` |
| `60` | `2.710x` | `-9.1` | `.............((((..((....)).))))...((((...(((....))).))))...` |
| `80` | `3.955x` | `-18.8` | `......((((((((((.((........)).))))).)))))(((((....(((.......)))..)))))..........` |
| `100` | `4.606x` | `-33.5` | `((((......(((((((((.((((.((((((((.((((.(((((..(.....)..))))).))))..)))))))).)))).).....)))))))))))).` |

## Exactness Notes

- Most stable structure matched exactly in all `81` cases.
- Most stable energy matched exactly in all `81` cases.
- Mean folding time matched exactly in `77/81` cases.
- The remaining `4` cases differed only by tiny floating-point-scale amounts:
  - `L087_S000`: `+8.2e-05`
  - `L095_S000`: `+9.44e-04`
  - `L096_S000`: `+1.0e-06`
  - `L097_S000`: `-8.1e-05`

These mean-fold-time differences are negligible relative to the `~500` time scale of the runs.
