The DTA files "composite_gallup_boot.dta", "cps_for_union_boot.dta", "state_year_union_ineq.dta" are the raw data from the paper:
Farber, H. S., D. Herbst, I. Kuziemko, and S. Naidu (2021). Unions and inequality over the twentieth century: New evidence from survey data. The Quarterly Journal of Economics 136 (3), 1325–1385.

The folder "original_code" contains the original code from the paper above with STATA code.

I used R to present the results.

To draw the plot of Top-10% Share vs Union density across years, the ACF plot,  and the Q-statistics plot, use 1121_c_plot.R.

To extract the covariate matrix Z from the raw data, and the corresponding AR1 coefficients in X and eps, use "1108_c.R".
In this file, stop at the line of 
"save(list = c("dat_target", "eps_perp_ori_wide", 
              "loading1_wide", "loading2_wide",
              "Q_res_c", "Q_Xperp_c", "Q_p_res_c", "Q_p_Xperp_c",
              "res_c_acf_mat", "res_c_pacf_mat", "Xperp_c_acf_mat", "Xperp_c_pacf_mat",
              "Xperp_ori_wide", "Z", 
              "Z.eta", "Z.gamm",
              "bet", "eta", "gamm",
              "northeast_states", "south_states", "statefip_ipt", "west_states", "year_ipt"),
     file = "1108_PermSimu_Setting_c.Rdata")"
The rest of the code is performing some permutation based simulation which I didn't use.


The "1108_PermSimu_Setting_c.Rdata" contains the covariate Z used in my simulation.

The "1121_c.R" runs the simulation used in the paper, i.e., constructing AR1 process of each state for both X and eps, and add aggregate shock on eps.
The "1121_c_parallel.R" runs the same simulation but use do.Parallel.
The simulation results are stored in "1121_c_AR1X&eps_epsShock_changeShock_L1.Rdata", "1121_c_AR1X&eps_epsShock_changeShock_L1_sameseed.Rdata",...
If the Rdata file name contains "sameseed", it means that within each value of shock strength in eps, the seed does not change.

The non-coverage and coverage probabilities plots are in "1121_c_plot.R".

Then I used HPC to calculate the results with 0621_HAC_satterthwaite_prewhite_cluster.R and 0621_HAC_array.sbatch. (The 0621.R file has some local code to try before HPC, but not complete for the df of PW method. The _PW parts should also have its corresponding df.) It still uses the "1108_PermSimu_Setting_c.Rdata" but add several finite sample adjustments. The Yale cluster simulation results are stored in results_0621 folder.

I finally use the HAC_HC1_Satterthwaite in the summary csv's of the result_0621 folder. The logic of Satterthwaite df can be found here: https://chatgpt.com/share/6a386b10-1098-83ea-aacb-c955246cc71a





