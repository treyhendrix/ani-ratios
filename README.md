#  No evidence for adaptive sex ratio adjustment in a cooperatively breeding bird with helpful helpers
 This is the code for the following manuscript currently under peer review: 
 
 >  Hendrix, T.C. and C. Riehl (In review) No evidence for adaptive sex ratio adjustment in a cooperatively breeding bird with helpful helpers.
 
 
## Project Abstract
The project investigates adaptive sex ratio adjustment in the Greater Ani (*Crotophaga major*) via model selection. Abstract from the manuscript:


> The repayment hypothesis for offspring sex allocation predicts that breeders in cooperatively breeding groups should overproduce the helping sex, particularly when they lack helpers. In birds that produce sexually size-dimorphic nestlings, sex allocation is further predicted to vary across hatching order to maximize brood survival. This could occur by biasing late-hatched nestlings towards either the energetically cheaper sex (the intra-brood sharing-out hypothesis) or the more competitive sex (the intra-brood competitive equilibrium hypothesis). Here we test these hypotheses using data from 553 nestlings in 109 broods of Greater Anis (*Crotophaga major*), a cooperatively breeding bird that breeds in groups composed of multiple reproductive pairs and non-breeding helpers. Helpers are predominantly male and increase the reproductive output of breeders; late-hatched nestlings are more vulnerable to starvation. Therefore, the repayment hypothesis predicts that groups without helpers should produce male-biased broods, whereas the intra-brood sharing-out and competitive equilibrium hypotheses predict that either females (the energetically cheaper sex) or males (the more competitive sex), respectively, should be overproduced at the end of the laying sequence. Contrary to these predictions, population-wide sex ratios did not differ significantly from 50:50, and we found no evidence for facultative sex ratio adjustment within broods by helper presence or by hatch order. These results support a growing consensus that facultative sex allocation is less widespread in birds than once thought, even in cooperatively breeding species with sex-biased helping behavior. 


## How to use this repository
This project is written exclusively in the R language. The code has been compiled so that one R script can reproduce the manuscript's statistical results, tables, figures, and supplemental materials. This R script is available in several formats that are identical in their content but have different use cases (described below). 

### Repository Structure
* `Data/`: upon publication, this folder will contain the .csv file necessary to reproduce analyses. It also includes a README file (`README_GRANH1_Data.md`) that details the contents of this .csv file. 
* `ani-ratios.Rproj`: an R project to facilitate running scripts. 
* `Output/`: A folder containing all the figures and tables presented in the manuscript and its supplemental materials. These are generated/updated by the R scripts in the `Scripts/` folder. 
* `Scripts/`: Contains three formats of the R code for the manuscript:
	+ `GRANH1_R_Code.md`: This is a Markdown document designed to be viewed on the GitHub website. It may be helpful to those interested in viewing the R code and its output (e.g., tables and figures) but not editing or running the code. The `GRANH1_R_Code_figures_for_markdown/` folder contains the .png files necessary to render the Markdown document.
	+ `GRANH1_R_Code.Rmd`: an R Markdown document that produces both the `GRANH1_R_Code.md` and `GRANH1_R_Code.R` documents when "knit." This document may interest those wanting to edit the code or run it in small code "chunks."
	+ `GRANH1_R_Code.R`: an alternative version of the above R Markdown document generated using the `purl` function from the `knitr` R package within the `GRANH1_R_Code.Rmd` document. It may be helpful to users who want to rerun analyses or edit code but prefer working with .R files.     

## Figshare 
Upon publication, the data, supplemental materials, and R code found in this GitHub repository will also be deposited in the Figshare repository as a collection with the following DOI: [https://doi.org/10.6084/m9.figshare.c.6350639](https://doi.org/10.6084/m9.figshare.c.6350639)
