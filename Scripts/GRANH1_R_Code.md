R Code for No evidence for adaptive sex ratio adjustment in a
cooperatively breeding bird with helpful helpers
================
Trey Hendrix
Updated: 2023-03-06

- <a href="#required-packages" id="toc-required-packages">Required
  Packages</a>
- <a href="#data" id="toc-data">Data</a>
- <a href="#preparing-the-data" id="toc-preparing-the-data">Preparing the
  data</a>
- <a href="#sample-sizes-reported-in-manuscript"
  id="toc-sample-sizes-reported-in-manuscript">Sample sizes reported in
  manuscript</a>
  - <a href="#abstract-and-methods-sections"
    id="toc-abstract-and-methods-sections">Abstract and methods sections</a>
  - <a href="#results-section" id="toc-results-section">Results section</a>
- <a href="#population-level-sex-ratio"
  id="toc-population-level-sex-ratio">Population-level sex ratio</a>
  - <a href="#figure-1" id="toc-figure-1">Figure 1</a>
- <a href="#facultative-sex-ratio-adjustment"
  id="toc-facultative-sex-ratio-adjustment">Facultative sex ratio
  adjustment</a>
  - <a href="#generating-list-of-candidate-models"
    id="toc-generating-list-of-candidate-models">Generating list of
    candidate models</a>
  - <a href="#model-selection" id="toc-model-selection">Model Selection</a>
  - <a href="#table-1" id="toc-table-1">Table 1</a>
  - <a href="#table-s1" id="toc-table-s1">Table S1</a>
  - <a href="#figure-2" id="toc-figure-2">Figure 2</a>
  - <a href="#table-s2" id="toc-table-s2">Table S2</a>
- <a href="#assessing-model-residuals-using-dharma"
  id="toc-assessing-model-residuals-using-dharma">Assessing model
  residuals using DHARMa</a>
  - <a href="#figure-s1" id="toc-figure-s1">Figure S1</a>

``` r
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  include = FALSE, 
  fig.path = "GRANH1_R_Code_figures_for_markdown/"
)
```

This document contains the R code used to produce the analyses, figures,
and supplemental material for the manuscript titled “No evidence for
adaptive sex ratio adjustment in a cooperatively breeding bird with
helpful helpers.” This manuscript is currently under peer review.

# Required Packages

This document was created using:

# Data

The data set analyzed in this study will be made publicly available in
the figshare repository upon publication. The DOI
[doi.org/10.6084/m9.figshare.21760979](http://doi.org/10.6084/m9.figshare.21760979)
has been reserved, and it will become active if the paper is accepted
for publication. We apologize for any inconvenience. The data file
necessary to reproduce these analyses is named “GRANH1_anon_221123.csv”
and, although it is not currently publicly available, its contents are
described in a “README” document in the “Data” folder of this
repository. The structure of the data frame is:

The structure of these data set is:

# Preparing the data

# Sample sizes reported in manuscript

## Abstract and methods sections

## Results section

# Population-level sex ratio

To determine if the nestlings in our study population overall deviated
from the expected 50:50 sex ratio, we consider two tests that differ in
the assumptions they make about the independence of a nestling’s sex.
These statistical tests (as well as the proportion of our study
population that is male) are reported in our results section under the
subheading “Population-level sex ratio.”

## Figure 1

Sex ratio of nestlings during each year of our study. Error bars show
95% confidence intervals. The dotted line indicates a 50:50 sex ratio,
and the number of nestlings sexed in each year is shown above the error
bars for each year.

# Facultative sex ratio adjustment

Here, we are interested in understanding which variables of interest, if
any, can predict the sex of individual nestlings. We use logistic
mixed-effects models with a logit link using the glmer function in the
lme4 package. All models included nestling sex (coded as “0” for females
and “1” for males) as the response variable as well as nest location and
year (nested within nest location) as random effects to control for
repeated measures. Fixed effects included hatch order, hatching
synchrony, brood size, clutch size, the presence of helpers, the number
of breeding pairs, and the interaction between synchrony and hatch
order. Please see the manuscript text for definitions of these
variables.

In our model selection we use a “best subsets” method. Candidate models
are compared to every possible model including a subset of the terms. We
consider 80 candidate models (listed below). We do not consider
candidate models that included the synchrony by hatch order interaction
term unless they also contain synchrony and hatch order separately as
fixed effects. Candidate models are evaluated based on Akaike’s
Information Criterion adjusted for finite sample sizes (AICc). Models
that differed from the top model by less than 2 AICc units are
considered to potentially have explanatory value unless they differed
from the top model by the addition of a single parameter.

## Generating list of candidate models

This code produces a list of 79 candidate models (80 if you consider a
null model with no fixed effects a candidate model):

## Model Selection

The code below creates a logistic mixed-effects models for each of our
79 candidate models generated above, compares the model to the null
model, and then ranks/sorts the models by their AICc scores.

## Table 1

Best subset models predicting the sex of individual nestlings. All
models are logistic mixed-effect models with nest location and year
included as a random intercept. Each model is compared to a null model
run on the same dataset including only random effects using a likelihood
ratio test. Models with ∆AICc \< 2 are shown as well as global and null
models. All continuous variables (i.e., clutch size, brood size,
synchrony) were scaled and centered to reduce multicollinearity.

## Table S1

Variance inflation factors for the fixed effects included in the global
model predicting the sex of individual nestlings. These estimates were
generated using the “vif” function from the “car” package (Fox and
Weisberg 2019).

## Figure 2

Sex ratio of nestlings with respect to (a) whether or not a helper was
observed at their nest and (b) position in hatch order. In both, error
bars show 95% confidence intervals, the dotted line indicates a 50:50
sex ratio, and the number of nestlings sexed in each year is shown above
each error bar.

## Table S2

Please not that the formatting of this table has not been reproduced
below. Instead, I present the information as 4 seperate tables. The
caption for this table presented in the Supplementary Information is:
“Parameters and random effect estimates for (a) the best-fit and (b)
global logistic mixed effects models predicting the sex of individual
ani nestlings. Significant parameters (p\<0.05) have been bolded. The
reference levels for the categorical fixed effects are: first hatching
(for hatch order) and two-pair group (for breeding pairs).”

# Assessing model residuals using DHARMa

Here, we use the “DHARMa” package to visualize the global model’s
residuals.

## Figure S1

Scaled residual plots of the global model predicting nestling sex.
