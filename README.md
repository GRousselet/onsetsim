# Using cluster-based permutation tests to estimate MEG/EEG onsets: how bad is it?

Code and data for simulation project on MEEG onset estimation.
[[Preprint]](https://www.biorxiv.org/content/10.1101/2023.11.13.566864v1)

## R code 

| Notebook | html version | Content | Figures |
| ----- | ----- | ----- | ----- |
|`fdr_demo.Rmd`| [fdr_demo](docs/fdr_demo.md) | Explain FDR using BH95 and BY01. Demonstrate weak control of the FWER. | 1 |
|`weakstrongfwer_demo.Rmd`| [weakstrongfwer_demo](docs/weakstrongfwer_demo.md) | Extend `fdr_demo.Rmd` to CS, CD and MAX algorithms. Consider if centring is necessary and sufficient to achieve strong FWER control. | 2 |
|`examples.Rmd`| [examples](docs/examples.md) | Illustrate noise, signal and statistical methods.| 3-4 |
|`onsetsim_eeg.Rmd`| [onsetsim_eeg](docs/onsetsim_eeg.md) | Simulate one participant using EEG-like noise with the approach from Yeung et al. (2004).| 5-6 |
|`onsetsim_eeg_half.Rmd`| [onsetsim_eeg_half](docs/onsetsim_eeg_half.md) | Same as `onsetsim_eeg` but using a Gaussian signal half in duration.| NA |
|`onsetsim_eeg_group.Rmd`| [onsetsim_eeg_group](docs/onsetsim_eeg_group.md) | Simulate 20 participants with random onsets.| 7-8 |
|`onsetsim_eeg_group_hb.Rmd`| [onsetsim_eeg_group_hb](docs/onsetsim_eeg_group_hb.md) | Example of hierarchical bootstrap onset inference using the CP method | 9-10 |
|`onsetsim_eeg_group_comp.Rmd`| [onsetsim_eeg_group_comp](docs/onsetsim_eeg_group_comp.md) | Compare two independent distributions of onsets estimated using the change point method.| 11 |
|`onsetsim_1overf.Rmd`| [onsetsim_1overf](docs/onsetsim_1overf.md) | Simulate one participant with 1/f noise. | NA |

## Abbreviations

BH95 = Benjamini-Hochberg (1995).  
BY01 = Benjamini-Yekutieli (2001). 
CS = cluster-sum.    
CD = cluster-depth.    
MAX = maximum statistics.   
CP = change point.  

## Matlab code

A subset of simulations were reproduced in Matlab. The code doesn't include methods BY01 and cluster-depth.  

| Notebook | Content |
| ----- | ----- |
|`onsetsim_1overf.m`| Illustrate 1/f noise & methods and perform simulations |
|`onsetsim_eeg.m`| Illustrate EEG-like noise (Yeung et al., 2004) & methods and perform simulations |
|`onsetsim_eeg_group.m`| Perform group simulation with 20 participants |

## Julia code

See replication of some of my results and link to cluster-depth implementation in Julia in [Benedikt Ehinger's blog](https://benediktehinger.de/blog/science/on-the-onsets-of-clusters-a-replication-of-rousselet-2023/).

## Dependencies

Extra code dependencies are in the `code` folder for R and in the `functions` folder for Matlab.

## Simulation results

All simulation results are in the `data` folder.  
Figures and results in the article can be reproduced without running the simulations.

