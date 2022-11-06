## The STRATEGIC study

This repository includes data and codes for reproducing the results of the manuscript:

**"Strategies to inTerrupt RAbies Transmission for the Elimination Goal by 2030 In China (STRATEGIC): a modelling study"**

[Link]: https://github.com/pkuepi/STRATEGIC

In the following sections, we will describe the purposes of each major directory and the scripts in that folder. **The scripts and data here are used for reproduction of the results in the paper only.** Scripts were tested on R-4.0.5.

The analysis code is modified from **the WHO rabies modelling consortium study**     

[Link]: https://github.com/katiehampson1978/rabies_PEP_access   

**Reference:** *The Lancet Infectious Diseases* 2019; 19(1): 102-11.

### 1.Analysis
The `analysis` .Rmd file is set to run the multivariate analysis and univariate analysis, and the results are saved at the `output` folder. If you want to run the codes, you need to: 

#### 1) Change the path

Please change the `path` to the location of the `STRATEGIC` folder you downloaded and run the first chunk to load the R function.

#### 2) Prepare R packages

You may need to install the following R packages if you have not done so yet:

- triangle
- dplyr

Then you can run the second chunk to load the packages.

#### 3) Load the data

Please run the third chunk to load the data needed for the analysis.

#### 4) Probabilistic sensitivity analyses (PSA)

Please run the fourth chunk to complete the simulation of PSA (multivariate analysis).

#### 5) One-way sensitivity analyses 

Please run the fifth chunk to complete the simulation of univariate analysis.

### 2.Tables and figures generation
The `plot` .Rmd file is set to summarize the results of `analysis`. Tables and figures are saved at the `figures and tables` folder. If you want to run the codes, you need to: 

#### 1) Change the path

Please change the `path` to the location of the `STRATEGIC` folder you downloaded and run the first chunk to load the R function.

#### 2) Prepare R packages

You may need to install the following R packages if you have not done so yet:

- dplyr
- tools
- ggplot2
- scales

And then you can run the second chunk to load the packages.

#### 3) Generate tables 

Please run the third to fifth chunk to generate the tables of the predicted human rabies deaths(Table 1 & 2), the incremental cost-effectiveness ratio (ICER) (Supplementary Table S5).

#### 4) Generate figures

Please run the sixth to ninth chunk to generate the figures of the predicted human rabies deaths(Figure 2) and the cost effectiveness analysis(Figure 3), the one-way sensitivity analyses (Supplementary Figure S3), the probabilistic sensitivity analyses (Supplementary Figure S4).

## If you are not interested in the details of this code, you can skip the following descriptions.

### 3.Descriptions of folders

#### Folder `R`
This folder contains major R functions used for our model fitting, parameter estimation and producing tables and figures. 

- `1.rabies probability calculation.R`: calculate the probability of the rabies infection of the person bitten by rabid dogs.
- `2.medical care probability calculation.R`: calculate the probability to seek, receive and complete the medical care of the bitten person.
- `3.DALYs calculation.R`: calculate the DALYs caused by human rabies.
- `4.time of strategy.R`: prepare the simulation period of different strategies.
- `5.univariate analysis.R`: simulate the one-way sensitivity analyses.
- `6.multivariate analysis.R`: simulate the probabilistic sensitivity analyses.
- `7.results summary.R`: summary the results of different scenarios by areas and years.
- `8.plot_function.R`: prepare the results for tables & figures generation.

#### Folder `data`

This folder contains the main data used in the study. 

- `data.csv`: This file contains human population, human-dog ratios, dog bite probabilities, medical care probabilities, and costs parameters.  
- `dogs_pop_traj.csv`: This file contains dog population parameters. 
- `baseline_incidence.csv`: This file contains rabies incidence of dogs.  
- `bio_data.csv`: This file contains the prevented effect parameters of post-exposure prophylaxis (PEP) treatment.  
- `biological_param_p_infect.csv`: This file contains probabilities of human rabies bitten on the different spots.  
- `DALY.csv`: This file contains parameters used to estimate DALYs .   
- `lifetable.csv`: This file contains life expectancy of different ages, used to calculate DALYs parameters.
- `rabies_traj.csv`: This file contains trends of rabies incidence caused by dog vaccination.  
- `vaccine_use.csv`: This file contains parameters of completeness of the PEP regimen.  


#### Folder `output `

This folder stores the results of `analysis.Rmd`. 

#### Folder `figures and tables `

This folder stores the files for main tables and figures. 

- `Tab1_model_summary_scenario.csv`: This file presented Table 1 in the manuscript.  
- `Tab2_model_year_scenario_all.csv`: This file presented Table 2 in the manuscript. 
- `Fig2_deaths_vaccinated_vials.pdf`: This file presented Figure 2 in the manuscript.  
- `Fig3_ICER_province.pdf`: This file presented Figure 3 in the manuscript.  
- `ICER_scenario.csv`,`multivariate.pdf`,`univariate.pdf` presented results in the supplementary files.
