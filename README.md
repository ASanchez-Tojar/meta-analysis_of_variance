[![DOI](https://zenodo.org/badge/195809400.svg)](https://zenodo.org/badge/latestdoi/195809400)

# Meta-analysing variances alongside means

This repository contains the R scripts and data used in the following study:

Alfredo Sánchez-Tójar, Nicholas P. Moran, Rose E. O'Dea, Klaus Reinhold, Shinichi Nakagawa. *In press*. **Illustrating the importance of meta-analysing variances alongside means in ecology and evolution**. Journal of Evolutionary Biology. (preprint: https://doi.org/10.32942/osf.io/yhfvk)

The **electronic supplementary material** is available online at: http://asanchez-tojar.github.io/meta-analysis_of_variance/SupplementaryMaterial

More information and materials available at the [OSF project](https://osf.io/yjua8/). DOI for this repository available via Zenodo (see above). For any further information, please contact: [Alfredo Sánchez-Tójar](https://scholar.google.co.uk/citations?hl=en&user=Sh-Rjq8AAAAJ&view_op=list_works&sortby=pubdate), email: alfredo.tojar@gmail.com

## Scripts:

`001_data_exploration.R`: this code was used for our initial exploration of the data, when we were trying to understand the reported and shared data, the calculations, and the strengths and weaknesses of the data. This script did not generate any results for the manuscript. R session info is available [here](https://github.com/ASanchez-Tojar/meta-analysis_of_variance/blob/master/data/data_exploration_R_session.txt).

`002_data_re-extraction.R`: this code was used for organizing the re-extraction of the data. The outputs of this script were later used to double-check and re-extract the data, which was performed by two observers. R session info is available [here](https://github.com/ASanchez-Tojar/meta-analysis_of_variance/blob/master/data_re-extraction/data_re-extraction_R_session.txt).

`003_data_preparation.R`: this code was used to put together the clean data from our different rounds of doublechecking, and to obtain information about effect size direction, group-level data, ratio-scale data, etc... The main output of this script is a database that can be used to calculate the effect size statistics. Be aware that this script contains a large amount of information extracted from many studies and used to make decisions about, for example, the direction of the effect size. It may not look pretty but it essentially provides a way of reproducing our decisions regarding effect size direction, and others. R session info is available [here](https://github.com/ASanchez-Tojar/meta-analysis_of_variance/blob/master/data_re-extraction/clean_data/data_preparation_R_session.txt).

`004_shared_control_implementation.R`: this code was used for dealing with share-control non-independence. R session info is available [here](https://github.com/ASanchez-Tojar/meta-analysis_of_variance/blob/master/data_re-extraction/clean_data/shared_control_R_session.txt).

`005_effect_size_calculation.R`: this code was used to calculate the effect size statistics of interest for our multilevel meta-analyses and meta-regressions. R session info is available [here](https://github.com/ASanchez-Tojar/meta-analysis_of_variance/blob/master/data_re-extraction/clean_data/effect_size_calculation_R_session.txt).

`006_phylogeny_reconstruction.R`: this code was used to build the phylogeny and estimate the phylogenetic relatedness among the species included in our analysis. It also corrects a couple of typos regarding species annotation. R session info is available [here](https://github.com/ASanchez-Tojar/meta-analysis_of_variance/blob/master/data_re-extraction/clean_data/phylogeny_reconstruction_R_session.txt).

`008_brms_meta-analysis_original.R`: this code was used to run all multilevel meta-analyses and meta-regressions. Models could not be added to the GitHub repository due to size limitations, but are available at the [Open Science Framework project](https://osf.io/zy7k2/). R session info is available [here](https://github.com/ASanchez-Tojar/meta-analysis_of_variance/blob/master/models/brms/brms_meta-analysis_R_session.txt).

`009_results_figures_and_tables.R`: this code was used to obtain the results of our analyses in the form of tables and figures. R session info is available [here](https://github.com/ASanchez-Tojar/meta-analysis_of_variance/blob/master/plots/results_figures_and_tables_R_session.txt).

`010_additional_calculations_figures_and_data.R`: this code was used to make a final calculations, including the review showing the use of meta-analysis of variance in ecology and evolution. R session info is available [here](https://github.com/ASanchez-Tojar/meta-analysis_of_variance/blob/master/literature_review/additional_R_session.txt).

`supplementary_information_meta-analysis_of_variance.Rmd`: this .Rmd file corresponds to the supplementary material, which output can be found here: http://asanchez-tojar.github.io/meta-analysis_of_variance/SupplementaryMaterial

### Notes:

26th of August, 2019: all R scripts were run on this day simply to make sure that all commands worked, and to find and fix any potential malfunction. That said, the outputs (e.g. datasets, plots) were not changed, and thus remained as originally. *PS: the missing 007 script was an attempt to compute variance-covariance matrices for the sampling errors that we had to discard for this project due to computational limitations.* 

See details of the licence of this repository in [LICENSE.txt](https://github.com/ASanchez-Tojar/meta-analysis_of_variance/blob/master/LICENSE.txt)
