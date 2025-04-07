
Section 4. 
Compound MSE(DP)-optimal designs search for the blocked case-study experiment

---------------------------------------------------------------------

1. R script titled "Figure_6_case_study.R":

The figure can be generated from the pre-calculated tables of designs' criteria values (default option), or the search can be run first by calling the "<...>_main.R" scripts -- this option is specified within the Figure script, the lines that runs the scripts can be just un-commented.

The figure is saved in the "output/" folder.

---------------------------------------------------------------------

2. If it is desirable to change the number of random starts -- this can be done in the "<..>_main.R" scripts, by changing the variable "Nstarts".

"MSEB_CP_main.R":   search for the designs with fixed center points, tau2 = 1.
"MSEB_CP_Q_main.R": search for the designs with fixed center points, tau2 = 1/q.

"MSEB_main.R":   search for the designs without fixed center points, tau2 = 1.
"MSEB_Q_main.R": search for the designs without fixed center points, tau2 = 1/q.


---------------------------------------------------------------------


3. Each "<...>_main.R" script creates an output file with the tables of designs' criteria values ("<..>_values.csv" files) and R workspace, saved in the corresponding folders within the "output/" folder.


---------------------------------------------------------------------


4. The "\functions" folder contains script files with the functions used for search and design evaluation:


  "Common.functions.R": functions that are common for all criteria
  "Common.functions_CP.R": functions that are common for all criteria, with fixing center points
  
  "MSEB_functions.R": criteria and search functions for the blocked design search
  "MSEB_CP_functions.R": criteria and search functions for the blocked design search, with fixed center points
  
  "Criteria_evaluation.R": functions for design evaluation


---------------------------------------------------------------------

	   