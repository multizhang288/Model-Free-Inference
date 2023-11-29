Folder Study 2 contains the simulation codes and results of Study 2, S.4.2, S.4.3 and S.4.4
 -Function-Study2.R: Core functions used

 -Distribution: Section 5.2.1

   -Data_Analysis.R: Codes of Table 3 (also covered in Table_Figures.Rmd) 
   -Simulation_Dist.R: Simulation demo
   -Intermediate Data: 25 .csv files like "Impact_Model_1_n_400_h_7_s_10_d_1_delta_1_choice1_Chi_choice2_Bspline_active_random_end.csv"
      -Model 1 with d =1, n = 400, h = 7, s = 10, delta  = 1 (signal strength), Chi: X~Chi dist, Bspline: Bspline Transformation
   -Results: Table 3
##############################
 -Sparsity: Section 5.2.2

   -Data_Analysis.R: Codes of Table 4 (also covered in Table_Figures.Rmd) 
   -Simulation_Sparsity.R: Simulation demo
   -Intermediate Data: 20 .csv files like "Impact_Model_1_n_400_h_7_s_10_d_1_delta_1_choice1_Chi_choice2_Bspline_active_random_end.csv"
      -Model 1 with d =1, n = 400, h = 7, s = 10, delta  = 1 (signal strength), Chi: X~Chi dist, Bspline: Bspline Transformation
   -Results: Table 4

#############################
 -Transformation: Section S.4.2

   -Data_Analysis.R: Codes of Table S.2 (also covered in Table_Figures.Rmd) 
   -Simulation_Transformation.R: Simulation demo
   -Intermediate Data: 20 .csv files like "Impact_Model_1_n_400_h_7_s_10_d_1_delta_1_choice1_Chi_choice2_Bspline_active_random_end.csv"
      -Model 1 with d =1, n = 400, h = 7, s = 10, delta  = 1 (signal strength), Chi: X~Chi dist, Bspline: Bspline Transformation
   -Results: Table S.2

############################
     -Extreme_d: Section S.4.3
     -Simulation_Extreme_d.R: Generate Table S.3
   -Results: Table S.3
##########################
 -Semi-synthetic: Section S.4.4

   -RatEyeExpression.txt: Real data of rat eyes gene measurements
   -Functions_Semi.R: Core funcitions
   -Semisynthetic.R: Simualtion script 
   -Results: Table S.4
##########################
 -Adversarial Example: Section S.4.6

   -Functions_Adversarial_MGT.R: Core funcitions with MGT
   -Simulation_Adversarial_MGT.R: Simualtion script with MGT
   -Functions_Adversarial_Nomgt.R: Core funcitions without MGT
   -Simulation_Adversarial_Nomgt.R: Simualtion script without MGT

      -Data_Analysis.R: Generate Figure S.3
        -Intermediate data MGT: The simulation output with MGT, 13 .csv files like "Impact_Adversarial__n_400_h_7_d_1_delta_1_df_0_Trans_Bspline_end"
        -Intermediate data Nomgt: The simulation output without MGT, 13 .csv files like "Impact_Adversarial__n_400_h_7_d_1_delta_1_df_0_Trans_Bspline_end"
        -Results: Figure S.3