# Bayesian re-analysis of the PREDESCI trial

The purpose of this project was to provide a Bayesian interpretation of the PREDESCI trial.  This was supplemented with probabilistic simulations to understand the implications of changed hazard estimates.  The following files are included to allow the analyses to be repeated in full.

## 00_predesci_data
This file allows analysis of data extracted from figure 2 and supporting text of the index publication.  Two endpoints are considered: 'binary', cause specific decompensation / liver-related mortality; and 'outcome', decompensation / liver-related before death from another cause to allow competing risk estimates.  In the estimated data, the outcomes are coded, 0 = censored, 1 = decompensation / liver-related death, 2 = transplant, and 3 = death before decompensation.

Output from this script is the cumulative incidence plot labelled with cause specific and competing risk summary statistics.

## 01_prior_distributions
The priors defined were as follows.  A moderate neutral ("skeptical") prior, a moderate optimistic prior, and a weak pessimistic prior.  This script plots these and combines them into a single figure.

## 02_bayes_reanalysis
This script runs the Bayesian reanalysis.  The first step is to fit the Bayesian survival model.  The posterior probabilities are then calculated for each of the three priors described above.  These models take some time to run, approximately 15 minutes overall.  The results, posterior probabilities and summaries of these, are exported for downstream analyses.

## 03_posterior_distributions
Outputs of the Bayes reanalyses for the three priors plotted.  Tabulated outcomes available from the summary tables, e.g."outstanding".

## 04_single_trial_simulation
This script runs a single trial including 2,000 patients with compensated cirrhosis at risk of decompensation.  The illustrated example has the annual probability of decompensation set at 0.05.  This can be adjusted as needed.  Simulations of this size take approximately 30 seconds to run.  There is, I'm sure, a more efficient way to run these simulations.  A plot of a single trial, annotated with the decompensation-free lifespan gain is produced.


