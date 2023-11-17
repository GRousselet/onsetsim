# Using cluster-based permutation tests to estimate MEG/EEG onsets: how bad is it?

Code and data for simulation project on MEEG onset estimation.
[[Preprint]](https://www.biorxiv.org/content/10.1101/2023.11.13.566864v1)

## R code 
| Notebook | Content | Figures |
| ----- | ----- | ----- |
|`examples.Rmd`| Illustrate noise, signal and statistical methods.| Fig 1 & 2 |
|`onsetsim_eeg.Rmd`| Simulate one participant using EEG-like noise with the approach from Yeung et al. (2004).| Fig 3 & 4 |
|`onsetsim_eeg_group.Rmd`| Simulate 20 participants with random onsets.| Fig 5 |
|`onsetsim_1overf.Rmd`| Simulate one participant with 1/f noise.| NA |

## Matlab code
| Notebook | Content |
| ----- | ----- |
|`onsetsim_1overf.m`| Illustrate 1/f noise & methods and perform simulations|
|`onsetsim_eeg.m`| Illustrate EEG-like noise (Yeung et al., 2004) & methods and perform simulations|
|`onsetsim_eeg_group.m`| Perform group simulation with 20 participants|

## Dependencies

Extra code dependencies are in the `code` folder for R and in the `functions` folder for Matlab.

## Simulation results

All simulation results are in the `data` folder. Figures and results in the article can be reproduced without running the simulations.