# Reproductibility codes for miRglmm paper

The codes in this repository, coupled with data located at XXXX, can be used to reproduce analyses and figures presented in (insert reference here).

# real_data_simulations.R

This code contains the simulation function and code to produce the datasets of monocyte data with aritifical effects induced. 

# prog_run_sim_array.R

This code will run all presented models on the specified array of simulated data. Run in conjunction with prog_run_sim_array.sbatch to array across multiple cores.

# prog_process_sim_results.R

This codes processes the model fits outputs produced by prog_run_sim_array.R and retrieves parameters of interest (beta, LRT p, etc). Calls functions found in process_results.R.

# extract_monocyte_data_subset.R, extract_bladdertestes_data_subset.R and study 89/extract_study89_data_subset.R

These codes extract sample types of interest from the miRNAome SummarizedExperiment object. 

# filter_mirna_monocyte.R, filter_mirna_bladder_testes.R, and study 89/filter_mirna_study89.R

These codes filters the extracted data to get a working dataset of samples. This also filters miRNA with insufficient expression to model (log(median CPM) <= 5).

# prog_run.R and study 89/prog_run.R

This code runs all presented models for a real biological dataset (not arrayed like in the simulation code).

# prog_process.R and study 89/prog_process.R

This code will extract parameters of interest from the model fit objects produced by prog_run.R. Calls functions found in process_results.R.

# process_results.R

This code contains functions which extract parameters of interest from the model fit objects produced by prog_run.R

# miRglmnb.R

This is a function for running the Negative Binomial GLM model. Called by prog_run.R and prog_run_sim_array.R. 

# figure codes directory

Codes to reproduce all figures from the manuscript. Filenames correspond to figure numbers. 

# ERCC directory

Codes to reproduce the ERCC analysis. filter_miRNA.R, prog_run.R and prog_process.R function as above but with different input files specified. 

prog_process_miRglmmreducedonly.R will only process results from miRglmm reduced models (that do not contain the random slope effect for sequence/isomiR). 

prog_run_filteragg.R and prog_process_filteragg.R will run the analysis that filters the isomiRs prior to aggregation.

prog_MSE.R will produce MSEs and coverage probabilities in a table. 

