# Lasso_Optimal_SSD
A repository of R code for Lasso Optimal SSD (supersaturated screening designs) evaluation.


First, the Design_Catalog folder has two subfolders corresponding to the different scenarios in the paper. 
In each of these subfolders there are two desings "d1" and "d2" that represent the PEDs from Sigh and Stufken et. al.
In the n14_p20_k5 folder, the other design represents the next best design that was not PED selected by HILS. 
In the n9_p10_k3 folder, there are two other designs representing the HILS optimal designs under either the integrated lambda measure, or the optimized lambda measure. These designs were not PED. 

The master file for most of the functions utilized in this repo is Lasso_optimal_SSD_function_library.R.

The CS_probabillity_heatmaps.R file produces the contour plots in Section 3 and can also be used to give the contour plots in Section 5 with slight parameter changes. 


The HILS_n9_p10.R file gives the code to run the HILS algorithm for Scenario 1 in Section 5 using the integrated lambda measure. 
The HILS_n9_p10_opt_lam_reeval.R file gives the code to rerank the designs produced by the previous file using the optimized lamnda measure. 

The HILS_n14_p20_opt_lam.R file gives the code to run the HILS algorithm for Scenario 2 in Section 5 using the optimized lambda measure. 
The HILS_n14_p20_int_reeval.R file gives the code to rerank the designs produced by the previous file using the integrated lamnda measure. 

Lastly, the two files with "comparison" in the title, take the designs in the Design_Catalog and plot them to produce the plots shown in Section 5. 
