# Bayesian re-analysis of the PREDESCI trial

The purpose of this project was to provide a Bayesian interpretation of the PREDESCI trial.  This was supplemented with probabilistic simulations to understand the implicatons of changed hazard estimates.  The following files are included to allow the analyses to be repeated in full.

## 00_predesci_data
This file allows analysis of data extracted from figure 2 of the index publication.  Two endpoints are considered: 'binary', cause specific decompensation / liver-related mortality; and 'outcome', decompensation / liver-related before death from another cause to allow competing risk estimates.  In the estimated data, the outcomes are coded, 0 = censored, 1 = decompensation / liver-related death, 2 = death before decompensation