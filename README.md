# EcoEvo_Competition_Adaptation
The folder contains 2 files:
- "Amicone_EcoEvo_Code.R"
- "Amicone_EcoEvo_Analysis.R"

The first is the main script to run the model described in the paper (Amicone and Gordo 2021, Evolution).
The parameters are explained in the script (lines 9-30) and can be changed by the users.
Notice that the code is such that multiple mutation rates can be given (line 20).
The code will then run each of the U conditions while keeping the other parameters constant.

As it is, the code models drift + selection as explained in the manuscript.
To run the neutral model (drift alone) or the pure selection model (no drift), just uncomment lines 148 or 149, respectively.

The full set of data will be output as .RData file named as chosen in line 28.
Data are stored as nested list. for example:
- allStrains[[1]] contains the #genotypes of all the simulated populations with the first given U, over time.
allStrains[[1]][2,13] is the number of genotypes of population 2 at generation 13.
- allAlphas[[2]][[4]] contains all the phenotypes that emerged in population 4 simulated with the second given U.
- allAlphas[[2]][[4]][[20]] contains only the alpha traits of genotype 20. 
- allHistory[[3]][[10]][[20]] contains the mutation history of genotype 20 of population 10 evolved under the third given U condition.

The computational times might be very long for conditions where N*U>10.

Once this code run, the data can be analysed and summarized by the second script: "Amicone_EcoEvo_Analysis.R".
This script takes as input the output of the previous script and outputs summary data:
- Alpha_Dynamics.RData, containing the trait sum over time (e.g. Fig. 2E of Amicone & Gordo 2021)
- M_Dynamics.RData, containing the number of genotypes over time (e.g. Fig. 3A of the same paper)
- U*_Sammary_Data.csv, one for each U condition, containing a table with the final statistics at the end of the simulations (Pi, Tajima's D, #clusters, etc.)
