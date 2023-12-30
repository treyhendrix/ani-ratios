R Code for “No evidence for adaptive sex ratio adjustment in a
cooperatively breeding bird with helpful helpers”
================
Trey Hendrix

Updated: 2023-12-30

- [Required Packages](#required-packages)
- [Data](#data)
- [Preparing the data](#preparing-the-data)
- [Sample sizes reported in article](#sample-sizes-reported-in-article)
  - [Abstract and methods sections](#abstract-and-methods-sections)
  - [Results section](#results-section)
- [Population-level sex ratio](#population-level-sex-ratio)
  - [Figure 1](#figure-1)
- [Facultative sex ratio adjustment](#facultative-sex-ratio-adjustment)
  - [Generating list of candidate
    models](#generating-list-of-candidate-models)
  - [Model Selection](#model-selection)
  - [Table 1](#table-1)
  - [Table S1](#table-s1)
  - [Figure S1](#figure-s1)
  - [Table S2](#table-s2)
  - [Figure 2](#figure-2)
  - [Table S2](#table-s2-1)
- [Assessing model residuals using
  DHARMa](#assessing-model-residuals-using-dharma)

This document contains the R code used to produce the analyses, figures,
and supplemental material for the article published in *Behavioral
Ecology and Sociobiology* titled “No evidence for adaptive sex ratio
adjustment in a cooperatively breeding bird with helpful helpers” (DOI:
[10.1007/s00265-023-03355-1](https://doi.org/10.1007/s00265-023-03355-1)).
The complete list of authors is Trey C. Hendrix and Christina Riehl.

# Required Packages

``` r
library(readr)
library(dplyr) 
library(lme4)
library(lmerTest)
library(broom.mixed)
library(MuMIn)
library(purrr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(DHARMa)
library(car)
library(forcats)
library(knitr)
library(callr)
library(rdryad)
library(simr)
```

This document was created using:

``` r
sessionInfo()
```

    ## R version 4.3.2 (2023-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 11 x64 (build 22631)
    ## 
    ## Matrix products: default
    ## 
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] simr_1.0.7          rdryad_1.0.0        callr_3.7.3        
    ##  [4] knitr_1.45          forcats_1.0.0       car_3.1-2          
    ##  [7] carData_3.0-5       DHARMa_0.4.6        tidyr_1.3.0        
    ## [10] gridExtra_2.3       ggplot2_3.4.4       stringr_1.5.1      
    ## [13] purrr_1.0.2         MuMIn_1.47.5        broom.mixed_0.2.9.4
    ## [16] lmerTest_3.1-3      lme4_1.1-35.1       Matrix_1.6-4       
    ## [19] dplyr_1.1.4         readr_2.1.4        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.0    RLRsim_3.1-8        fastmap_1.1.1      
    ##  [4] digest_0.6.33       mime_0.12           lifecycle_1.0.4    
    ##  [7] processx_3.8.2      magrittr_2.0.3      compiler_4.3.2     
    ## [10] rlang_1.1.2         tools_4.3.2         plotrix_3.8-4      
    ## [13] utf8_1.2.4          yaml_2.3.7          curl_5.1.0         
    ## [16] plyr_1.8.9          abind_1.4-5         httpcode_0.3.0     
    ## [19] withr_2.5.2         numDeriv_2016.8-1.1 grid_4.3.2         
    ## [22] stats4_4.3.2        fansi_1.0.5         colorspace_2.1-0   
    ## [25] future_1.33.0       globals_0.16.2      scales_1.3.0       
    ## [28] iterators_1.0.14    MASS_7.3-60         crul_1.4.0         
    ## [31] cli_3.6.1           rmarkdown_2.25      generics_0.1.3     
    ## [34] rstudioapi_0.15.0   binom_1.1-1.1       tzdb_0.4.0         
    ## [37] minqa_1.2.6         splines_4.3.2       parallel_4.3.2     
    ## [40] vctrs_0.6.4         boot_1.3-28.1       jsonlite_1.8.7     
    ## [43] hms_1.1.3           pbkrtest_0.5.2      listenv_0.9.0      
    ## [46] glue_1.6.2          parallelly_1.36.0   nloptr_2.0.3       
    ## [49] codetools_0.2-19    ps_1.7.5            stringi_1.8.2      
    ## [52] gtable_0.3.4        munsell_0.5.0       tibble_3.2.1       
    ## [55] furrr_0.3.1         pillar_1.9.0        rappdirs_0.3.3     
    ## [58] htmltools_0.5.7     R6_2.5.1            hoardr_0.5.3       
    ## [61] evaluate_0.23       lattice_0.21-9      backports_1.4.1    
    ## [64] broom_1.0.5         Rcpp_1.0.11         zip_2.3.0          
    ## [67] nlme_3.1-163        mgcv_1.9-0          xfun_0.41          
    ## [70] pkgconfig_2.0.3

# Data

The data set analyzed in this study is available in two locations:

- The “Data” folder of the present GitHub repository,
  [ani-ratios](https://github.com/treyhendrix/ani-ratios).
- In the Figshare repository, DOI:
  [10.6084/m9.figshare.c.6350639.v1](https://doi.org/10.6084/m9.figshare.c.6350639.v1).

The data file necessary to reproduce these analyses is named
“GRANH1_anon_221123.csv”. Its contents are described in a “README”
document in the “Data” folder of the ani-ratios repo and in the file
description on Figshare. The structure of the data frame is:

``` r
working_directory <- getwd()
project_file_path <- str_remove(working_directory, pattern = stringr::fixed("/Scripts")) # This code might need to be modified on your machine to produce the correct file path

file_path_to_GRANH1_csv <- paste0(project_file_path, "/Data/GRANH1_anon_221123.csv") # Using the copy of the csv on the ani-ratios repo
```

The structure of these data set is:

``` r
df <- read_csv(file_path_to_GRANH1_csv)
glimpse(df, show_col_types = FALSE)
```

    ## Rows: 553
    ## Columns: 11
    ## $ Year            <dbl> 2007, 2007, 2007, 2007, 2007, 2007, 2007, 2007, 2007, …
    ## $ Location        <chr> "F4", "F4", "F4", "F4", "L2", "L2", "L2", "L2", "H2", …
    ## $ Band_number     <chr> "O88", "Z98", "R12", "T01", "X37", "V24", "S61", "S54"…
    ## $ Sex             <dbl> 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, …
    ## $ Sex_ratio       <dbl> 0.2500000, 0.2500000, 0.2500000, 0.2500000, 0.5000000,…
    ## $ Hatch_order     <chr> "First", "Middle", "Middle", "Last", "First", "Middle"…
    ## $ Synchrony       <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, NA, NA, NA, …
    ## $ Brood_size      <dbl> 6, 6, 6, 6, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, NA, NA, NA, …
    ## $ Clutch_size     <dbl> 8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 6, 6, 6, 6, NA, NA, NA, …
    ## $ Helper_presence <dbl> 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ Breeding_pairs  <chr> "Two", "Two", "Two", "Two", "Two", "Two", "Two", "Two"…

# Preparing the data

``` r
# Adding a year_location variable to df to denote unique broods
df <- df %>%
  mutate(Year_location = paste(Year, Location, sep = "_"))

# Standardizing continuous fixed effects to reduce multicollinearity
df <- df %>%
  mutate(across(c(Brood_size, Clutch_size, Synchrony), 
                .fns = list(scaled = ~scale(., center = TRUE, scale = TRUE)[,1])))

# Setting factor levels to aid in interpretation
df <- df %>% 
  mutate(Hatch_order = factor(Hatch_order, levels = c("First", "Middle", "Last")), 
         Breeding_pairs = factor(Breeding_pairs, levels = c("Two", "Three")), 
         Hatch_order = relevel(Hatch_order, ref = "First"), 
         Breeding_pairs = relevel(Breeding_pairs, ref = "Two"))

# Fixed effects to be used in models
fixed_effects <- c("Hatch_order", "Brood_size_scaled", 
                   "Clutch_size_scaled","Helper_presence", 
                   "Synchrony_scaled",  "Breeding_pairs", 
                   "Hatch_order:Synchrony_scaled")

# Saving a complete case version of the data frame for use in models
dfcc <- df[complete.cases(select(df, fixed_effects[1:6])),]
```

# Sample sizes reported in article

## Abstract and methods sections

``` r
# All nestlings with sex information 
n_nestlings <- df %>% pull(Band_number) %>% length()
n_broods <- df %>% pull(Year_location) %>% unique() %>% length()
n_locations <- df %>% pull(Location) %>% unique() %>% length() 

paste(n_nestlings, "nestlings of known sex from", n_broods, "broods and", n_locations, "locations.")
```

    ## [1] "553 nestlings of known sex from 109 broods and 67 locations."

``` r
# Sample sizes for complete case analysis
n_cc_nestlings <- dfcc %>% pull(Band_number) %>% length()
n_cc_broods <- dfcc %>% pull(Year_location) %>% unique() %>% length()
n_cc_locations <- dfcc %>% pull(Location) %>% unique() %>% length() 

paste("For complete case analysis:", n_cc_nestlings, "nestlings,", n_cc_broods, "broods, and", n_cc_locations, "locations.")
```

    ## [1] "For complete case analysis: 472 nestlings, 87 broods, and 56 locations."

## Results section

``` r
# On average, how many nestlings were sampled in years where a female bias was observed vs. when no bias was observed?
df %>%
  group_by(Year) %>%
  summarize(Mean_sex_ratio = mean(Sex_ratio), 
            n = n(),
            se = sd(Sex_ratio)/sqrt(n)) %>%
  ungroup() %>% 
  mutate(Upper_CI =  Mean_sex_ratio + 1.96*se, 
         Lower_CI =  Mean_sex_ratio - 1.96*se, 
         Bias = case_when(Lower_CI > 0.5 ~ "Male", 
                          Upper_CI < 0.5 ~ "Female", 
                          TRUE ~ "None")) %>% 
  group_by(Bias) %>% 
  summarize(Mean_n = mean(n)) %>% 
  ungroup() %>% 
  mutate(Mean_n = round(Mean_n)) %>% 
  kable()
```

| Bias   | Mean_n |
|:-------|-------:|
| Female |     28 |
| None   |     77 |

# Population-level sex ratio

To determine if the nestlings in our study population overall deviated
from the expected 50:50 sex ratio, we consider two tests that differ in
the assumptions they make about the independence of a nestling’s sex.
These statistical tests (as well as the proportion of our study
population that is male) are reported in our results section under the
subheading “Population-level sex ratio.”

``` r
# Proportion of our study population that is male
mean(df$Sex) # percentage male  
```

    ## [1] 0.4810127

``` r
sum(df$Sex) # total males
```

    ## [1] 266

``` r
length(df$Sex) # total nestlings
```

    ## [1] 553

``` r
# Binomial test: assumes independence of nestlings
binom.test(x = sum(df$Sex), n = length(df$Sex), p = 0.5, alternative = "two.sided", conf.level = 0.95) 
```

    ## 
    ##  Exact binomial test
    ## 
    ## data:  sum(df$Sex) and length(df$Sex)
    ## number of successes = 266, number of trials = 553, p-value = 0.3951
    ## alternative hypothesis: true probability of success is not equal to 0.5
    ## 95 percent confidence interval:
    ##  0.4386631 0.5235661
    ## sample estimates:
    ## probability of success 
    ##              0.4810127

``` r
# Wilcoxon signed rank test: A paired test that does not assume that nestlings are independent
wilcoxon_df <- df %>%
  group_by(Year_location) %>%
  summarize(Prop_male = mean(Sex == 1), Prop_female = mean(Sex == 0)) %>%
  ungroup()

wilcox.test(x = wilcoxon_df$Prop_male, y = wilcoxon_df$Prop_female, alternative = "two.sided", paired = TRUE, conf.level = 0.95) 
```

    ## 
    ##  Wilcoxon signed rank test with continuity correction
    ## 
    ## data:  wilcoxon_df$Prop_male and wilcoxon_df$Prop_female
    ## V = 1882, p-value = 0.3163
    ## alternative hypothesis: true location shift is not equal to 0

To account for the effects of brood and year, we also assess the
population-level sex ratio by making a generalized linear mixed-effects
model with a logit link (lme4 package). The response variable is
nestling sex (coded as “0” for females and “1” for males). Nest location
and year (nested within nest location) are included as random effects.
We do not include any fixed effects and calculate the odds ratio of the
intercept, which may deviate from 1 if a population-level sex bias is
present.

``` r
# Only use data from completely sampled broods 
complete_brood_ids <- df %>% 
  count(Year_location, Brood_size) %>% 
  mutate(Complete_brood = Brood_size == n) %>% 
  filter(Complete_brood) %>% 
  pull(Year_location)

dfcb <- df %>%
  filter(Year_location %in% complete_brood_ids)

# Function for creating a single model
logistic_model_maker <- function(dv, iv, re, df = dfcc){
  formula <- paste(dv, "~", iv, "+", re, sep = " ") %>% as.formula()
  safe_glmer <- safely(glmer)
  safe_output <- safe_glmer(formula, data = df, family = "binomial", control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  model <- safe_output$result
}


population_sex_bias_model <- logistic_model_maker(dv = "Sex", iv = "1", re = "(1 | Year/Location)", df = dfcb)
summary(population_sex_bias_model)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: Sex ~ 1 + (1 | Year/Location)
    ##    Data: ..2
    ## Control: ..4
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    380.2    391.0   -187.1    374.2      267 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.0149 -1.0149  0.9853  0.9853  0.9853 
    ## 
    ## Random effects:
    ##  Groups        Name        Variance Std.Dev.
    ##  Location:Year (Intercept) 0        0       
    ##  Year          (Intercept) 0        0       
    ## Number of obs: 270, groups:  Location:Year, 49; Year, 11
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)  0.02963    0.12173   0.243    0.808
    ## optimizer (bobyqa) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

``` r
odds_ratio_of_intercept <- broom.mixed::tidy(population_sex_bias_model) %>% 
  filter(term == "(Intercept)") %>% 
  mutate(upper = estimate + 1.96*std.error, 
         lower = estimate - 1.96*std.error) %>% 
  mutate(across(c(estimate, upper, lower), ~ round(exp(.x), 2))) %>% 
  select(term, estimate, lower, upper)

odds_ratio_of_intercept
```

    ## # A tibble: 1 × 4
    ##   term        estimate lower upper
    ##   <chr>          <dbl> <dbl> <dbl>
    ## 1 (Intercept)     1.03  0.81  1.31

## Figure 1

Sex ratio of nestlings during each year of our study. Error bars show
95% confidence intervals. The dotted line indicates a 50:50 sex ratio,
and the number of nestlings sexed in each year is shown above the error
bars for each year.

``` r
# ggplot theme that will be reused for plots 
theme_simple <- theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(), 
                      axis.line = element_line(colour = "black"), 
                      legend.key = element_rect(fill = NA), 
                      text = element_text(size = 15)) 


figure_1 <- df %>%
  group_by(Year) %>%
  summarize(Mean_sex_ratio = mean(Sex_ratio), 
            n = n(),
            se = sd(Sex_ratio)/sqrt(n)) %>%
  ggplot(aes(x = Year, y = Mean_sex_ratio)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = Mean_sex_ratio - 1.96*se, ymax = Mean_sex_ratio + 1.96*se)) + 
  theme_simple + 
  ylab("Proportion of males in nest") + 
  xlab("Year") + 
  geom_hline(yintercept = 0.5, lty = 2, color = "black") + 
  scale_x_continuous(breaks = seq(2007, 2019, 1), labels = paste0("'", str_pad(7:19, width = 2, side = "left", pad = "0"))) + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) + 
  geom_text(aes(label = paste(n), y = Mean_sex_ratio + 1.96*se + 0.10), size = 4, color = "black")

figure_1
```

![](GRANH1_R_Code_figures_for_markdown/Figure%201-1.png)<!-- -->

``` r
output_file_path <- paste0(project_file_path, "/Output")

ggsave(paste0(output_file_path, "/Figure_1.pdf"), figure_1, width = 7, height = 3)
```

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
order. Please see the paper’s text for definitions of these variables.

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

``` r
# Function for generating possible combinations of fixed effects 
fixed_effect_combiner <- function(fixed_effects){
  
  fixed_effect_combinations <- lapply(as.list(1:length(fixed_effects)), function(x){combn(fixed_effects, m = x)})
  
  column_getter <- function(l, n){
    lapply(1:ncol(l[[n]]), function(x){l[[n]][,x]})
  }
  
  fixed_effect_list <- lapply(1:length(fixed_effect_combinations), function(x){column_getter(l = fixed_effect_combinations, n = x)}) %>%
    purrr::flatten() %>% 
    lapply(function(x){paste(x, collapse = " + ")}) %>% 
    unlist()
  
  return(fixed_effect_list)
}

# Function for removing models that include interactions made from effects not in the model 
interaction_indexer <- function(fixed_effects_for_one_model){
  
  plus_remover <- function(x){
    str_split(x, pattern = "\\+") %>% 
      unlist() %>% 
      str_trim(side = "both")
  }
  
  all_fixed_effects <- plus_remover(fixed_effects_for_one_model)
  interactions <- all_fixed_effects[str_detect(all_fixed_effects, ":")]
  noninteractions <- all_fixed_effects[!str_detect(all_fixed_effects, ":")]
  
  if(length(interactions) == 0){
    outcome <- TRUE
  } else {
    interactions_split <- lapply(as.list(interactions), function(x){str_split(x, ":")}) %>% 
      unlist() %>%
      unique()
    
    inclusion_test <- mean(interactions_split %in% noninteractions)
    
    outcome <- ifelse(inclusion_test == 1, TRUE, FALSE)
  }
  
  return(outcome)
  
}

# Combining fixed effects
index_of_combinations_to_keep <- lapply(fixed_effect_combiner(fixed_effects), interaction_indexer) %>% unlist() 

candidate_models <- fixed_effect_combiner(fixed_effects)[index_of_combinations_to_keep]
```

## Model Selection

The code below creates a logistic mixed-effects models for each of our
79 candidate models generated above, compares the model to the null
model, and then ranks/sorts the models by their AICc scores.

``` r
# Additional modeling functions to be reused

# Cleaning up terms from the model to make them easier to read
term_cleaner <- function(string){
  new <- str_replace_all(string, stringr::fixed("Hatch_orderMiddle"), "Middle nestling")
  new <- str_replace_all(new, stringr::fixed("Hatch_orderLast"), "Last nestling")
  new <- str_replace_all(new, stringr::fixed("Hatch_orderMiddle:Synchrony_scaled"), "Middle nestling:Synchrony")
  new <- str_replace_all(new, stringr::fixed("Hatch_orderLast:Synchrony_scaled"), "Last nestling:Synchrony")
  new <- str_replace_all(new, stringr::fixed("Helper_presence"), "Helper")
  new <- str_replace_all(new, stringr::fixed("Brood_size_scaled"), "Brood size")
  new <- str_replace_all(new, stringr::fixed("Clutch_size_scaled"), "Clutch size")
  new <- str_replace_all(new, stringr::fixed("Synchrony_scaled"), "Synchrony")
  new <- str_replace_all(new, stringr::fixed("Breeding_pairsTwo"), "Two-pair group")
  new <- str_replace_all(new, stringr::fixed("Breeding_pairsThree"), "Three-pair group")
  new <- str_replace_all(new, stringr::fixed("_"), " ") # Remove underscores
  return(new)
}


# Function for visualizing a single model
logistic_visualizer <- function(model){
  
  df <- model %>%
    tidy() %>%
    filter(effect == "fixed") %>%
    select(- effect) %>%
    filter(!(term %in% c("(Intercept)", "sd_(Intercept).Location:Year", "sd_(Intercept).Year"))) %>%
    mutate(term = term_cleaner(term)) %>%
    mutate(term = fct_reorder(term, estimate)) %>%
    mutate(lower = estimate - 1.96*std.error,
           upper = estimate + 1.96*std.error)
  
  df %>%
    ggplot(aes(x = estimate, y = term)) + 
    geom_point() + 
    geom_linerange(aes(xmin = lower, xmax = upper)) + 
    labs(x = "Model Estimate", y = "Fixed Effect") + 
    geom_vline(xintercept = 0, lty = 2, color = "red") +
    theme_simple
}

# Function for running a list of models, ranking them, and producing a summary table 
logistic_model_comparer <- function(dv, iv, re, df = df) {
  
  # Complete case analysis
  fixed_effects_vector <- iv %>% str_split(pattern = stringr::fixed(" + ")) %>% unlist()
  fixed_effects_vector <- fixed_effects_vector[str_detect(fixed_effects_vector, pattern = stringr::fixed(":"), negate = TRUE)]
  df_index <- df[,fixed_effects_vector]
  df <- df[complete.cases(df_index),]
  
  # Running model
  formula <- paste(dv, "~", iv, "+", re, sep = " ") %>% as.formula()
  safe_glmer <- safely(glmer)
  safe_output <- safe_glmer(formula, data = df, family = "binomial", control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  model <- safe_output$result
  
  # Null model 
  null_formula <- paste(dv, "~", "1",  "+", re, sep = " ") %>% as.formula()
  null_model <- safe_glmer(null_formula, data = df, family = "binomial", control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))$result
  
  # Making dataframe 
  model_output <- ifelse(length(safe_output$result) > 0, AICc(safe_output$result), NA_real_)
  error_output <- ifelse(length(safe_output$error) > 0, as.character(safe_output$error), "No error")
  result_df <- tibble(Formula = iv,
                      AICc = model_output,
                      Error = error_output)
  
  # ANOVA
  if(error_output == "No error") {
    anova_df <- anova(model, null_model) %>% tidy() %>% filter(term == "model") %>% select(-c("term", "df"))
    
    final_df <- result_df %>% 
      bind_cols(anova_df)
    
    
  } else {
    final_df <- result_df
  }
  
  return(final_df)
}


add_evidence_ratio_logistic <- function(dv, odf, re, df = dfcc){
  
  null_formula <- paste(dv, "~", "1",  "+", re, sep = " ") %>% as.formula()
  null_model <- glmer(null_formula, data = df, family = "binomial", control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  
  df <- data.frame(Formula = "Null", 
                   AICc = AICc(null_model),
                   logLik = logLik(null_model)) %>% 
    bind_rows(odf) %>%
    arrange(AICc) %>% 
    mutate(delta_AICc = abs(AICc - min(AICc)),
           rel_likelihood = exp(-1/2 * delta_AICc),
           w = rel_likelihood/(sum(rel_likelihood))) %>%
    arrange(desc(w)) %>% 
    mutate(w_rank = 1:nrow(.))
  
  null_w <- df %>%
    filter(Formula == "Null") %>%
    pull(w)
  
  null_logLik <- df %>%
    filter(Formula == "Null") %>%
    pull(logLik)
  
  final_df <- df %>%
    mutate(pseudo_R2 = 1 - logLik/null_logLik)
  
  
  final_df$evidence_ratio <- NA
  
  for(i in 1:(nrow(final_df)-1)){
    final_df$evidence_ratio[i] <- final_df$w[i]/final_df$w[i+1]
  }
  
  final_df <- final_df %>%
    arrange(AICc) %>%
    mutate(evidence_ratio = round(evidence_ratio, digits = 3)) %>%
    select(-w_rank)
  
  return(final_df)
}
```

``` r
# Global model 
global_model <- logistic_model_maker(dv = "Sex", iv = candidate_models[length(candidate_models)], re = "(1|Year/Location)", df = dfcc)
logistic_visualizer(global_model) 
```

![](GRANH1_R_Code_figures_for_markdown/Running%20Models-1.png)<!-- -->

``` r
summary(global_model) 
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: Sex ~ Hatch_order + Brood_size_scaled + Clutch_size_scaled +  
    ##     Helper_presence + Synchrony_scaled + Breeding_pairs + Hatch_order:Synchrony_scaled +  
    ##     (1 | Year/Location)
    ##    Data: ..2
    ## Control: ..4
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    664.9    714.8   -320.4    640.9      460 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.4312 -0.8985 -0.7179  0.9731  1.4307 
    ## 
    ## Random effects:
    ##  Groups        Name        Variance Std.Dev.
    ##  Location:Year (Intercept) 0.1122   0.3349  
    ##  Year          (Intercept) 0.0000   0.0000  
    ## Number of obs: 472, groups:  Location:Year, 87; Year, 11
    ## 
    ## Fixed effects:
    ##                                     Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)                         0.007552   0.188072   0.040   0.9680  
    ## Hatch_orderMiddle                  -0.371700   0.221110  -1.681   0.0928 .
    ## Hatch_orderLast                    -0.481147   0.269199  -1.787   0.0739 .
    ## Brood_size_scaled                   0.048410   0.163568   0.296   0.7673  
    ## Clutch_size_scaled                  0.110012   0.180163   0.611   0.5414  
    ## Helper_presence                     0.195424   0.219117   0.892   0.3725  
    ## Synchrony_scaled                   -0.208214   0.157569  -1.321   0.1864  
    ## Breeding_pairsThree                 0.173583   0.310100   0.560   0.5756  
    ## Hatch_orderMiddle:Synchrony_scaled  0.287487   0.219774   1.308   0.1908  
    ## Hatch_orderLast:Synchrony_scaled    0.113054   0.270425   0.418   0.6759  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Htch_M Htch_L Brd_s_ Cltc__ Hlpr_p Synch_ Brdn_T H_M:S_
    ## Htch_rdrMdd -0.544                                                        
    ## Htch_rdrLst -0.430  0.379                                                 
    ## Brd_sz_scld -0.012 -0.134  0.101                                          
    ## Cltch_sz_sc  0.218  0.071 -0.066 -0.633                                   
    ## Helpr_prsnc -0.351 -0.039 -0.054  0.014 -0.019                            
    ## Synchrny_sc  0.079 -0.001 -0.042 -0.200  0.110  0.071                     
    ## Brdng_prsTh -0.441  0.001 -0.002  0.058 -0.514 -0.003 -0.166              
    ## Htch_rdM:S_ -0.046 -0.088  0.026  0.033 -0.089  0.015 -0.620  0.041       
    ## Htch_rdL:S_ -0.016  0.021  0.164 -0.011 -0.034  0.004 -0.487 -0.004  0.368
    ## optimizer (bobyqa) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

``` r
# Model selection (summary table)
models_df <- lapply(as.list(candidate_models), 
                    function(x){logistic_model_comparer(dv = "Sex", iv = x, re = "(1|Year/Location)", df = dfcc)}) %>%
  bind_rows() %>%
  arrange(AICc) %>%
  mutate(delta_AICc = abs(AICc - min(AICc))) %>%
  select(Formula, delta_AICc, AICc, everything())

models_df <- add_evidence_ratio_logistic(dv = "Sex", odf = models_df, re = "(1|Year/Location)", df = df)
```

## Table 1

Best subset models predicting the sex of individual nestlings. All
models are logistic mixed-effect models with nest location and year
included as a random intercept. Each model is compared to a null model
run on the same dataset including only random effects using a likelihood
ratio test. Models with ∆AICc \< 2 are shown as well as global and null
models. All continuous variables (i.e., clutch size, brood size,
synchrony) were scaled and centered to reduce multicollinearity.

``` r
table_1 <- models_df %>%
  mutate(Formula = ifelse(Formula == candidate_models[length(candidate_models)], "Global", Formula)) %>% 
  filter(delta_AICc < 2 | Formula %in% c("Global", "Null")) %>% 
  mutate(Formula = term_cleaner(Formula))

# Removing uninformative parameters (Arnold 2010)
table_1 <- table_1 %>% 
  filter(!(Formula %in% c("Hatch order + Clutch size + Helper", "Hatch order + Clutch size + Synchrony", "Hatch order + Clutch size + Breeding pairs"))) 

# Tidying up table 1
table_1 <- table_1 %>%
  select(Fixed_effects = Formula, delta_AICc, AICc, chi_sq = statistic, P_value = p.value, w, pseudo_R2) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3)))

table_1 %>% 
  kable()
```

| Fixed_effects                | delta_AICc |    AICc | chi_sq | P_value |     w | pseudo_R2 |
|:-----------------------------|-----------:|--------:|-------:|--------:|------:|----------:|
| Hatch order + Clutch size    |      0.000 | 656.314 |  7.027 |   0.071 | 0.058 |     0.157 |
| Clutch size                  |      0.133 | 656.447 |  2.799 |   0.094 | 0.054 |     0.152 |
| Hatch order                  |      0.703 | 657.017 |  4.272 |   0.118 | 0.041 |     0.154 |
| Hatch order + Breeding pairs |      0.952 | 657.266 |  6.075 |   0.108 | 0.036 |     0.156 |
| Breeding pairs               |      1.062 | 657.376 |  1.870 |   0.172 | 0.034 |     0.150 |
| Hatch order + Brood size     |      1.266 | 657.580 |  5.761 |   0.124 | 0.031 |     0.156 |
| Brood size                   |      1.425 | 657.739 |  1.507 |   0.220 | 0.028 |     0.150 |
| Clutch size + Helper         |      1.470 | 657.784 |  3.505 |   0.173 | 0.028 |     0.153 |
| Clutch size + Synchrony      |      1.652 | 657.966 |  3.323 |   0.190 | 0.025 |     0.152 |
| Hatch order + Helper         |      1.967 | 658.281 |  5.060 |   0.167 | 0.022 |     0.155 |
| Global                       |      9.242 | 665.556 | 10.284 |   0.328 | 0.001 |     0.161 |
| Null                         |    113.993 | 770.307 |     NA |      NA | 0.000 |     0.000 |

``` r
write_csv(table_1, paste0(output_file_path, "/Table_1.csv"))
```

``` r
# Top performing model (not significant)
top_model <- logistic_model_maker(dv = "Sex", iv = "Hatch_order + Clutch_size_scaled", re = "(1|Year/Location)", df = dfcc) 
logistic_visualizer(top_model)
```

![](GRANH1_R_Code_figures_for_markdown/Top%20model-1.png)<!-- -->

``` r
summary(top_model) 
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: Sex ~ Hatch_order + Clutch_size_scaled + (1 | Year/Location)
    ##    Data: ..2
    ## Control: ..4
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    656.1    681.1   -322.1    644.1      466 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.3923 -0.8982 -0.7378  0.9678  1.3938 
    ## 
    ## Random effects:
    ##  Groups        Name        Variance Std.Dev.
    ##  Location:Year (Intercept) 0.1517   0.3895  
    ##  Year          (Intercept) 0.0000   0.0000  
    ## Number of obs: 472, groups:  Location:Year, 87; Year, 11
    ## 
    ## Fixed effects:
    ##                    Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)          0.1226     0.1568   0.781   0.4346  
    ## Hatch_orderMiddle   -0.3517     0.2178  -1.615   0.1063  
    ## Hatch_orderLast     -0.4742     0.2644  -1.793   0.0729 .
    ## Clutch_size_scaled   0.1802     0.1069   1.686   0.0918 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Htch_M Htch_L
    ## Htch_rdrMdd -0.678              
    ## Htch_rdrLst -0.554  0.405       
    ## Cltch_sz_sc  0.004 -0.050  0.017
    ## optimizer (bobyqa) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

## Table S1

Variance inflation factors for the fixed effects included in the global
model predicting the sex of individual nestlings. These estimates were
generated using the “vif” function from the “car” package (Fox and
Weisberg 2019).

``` r
table_s1 <- vif(global_model)  %>%
  as_tibble() %>% 
  mutate(Predictor = row.names(vif(global_model)), 
         Predictor = term_cleaner(Predictor), 
         VIF = round(GVIF, 3)) %>% 
  select(Predictor, VIF)

table_s1 %>% 
  kable()
```

| Predictor             |   VIF |
|:----------------------|------:|
| Hatch order           | 1.121 |
| Brood size            | 2.462 |
| Clutch size           | 3.014 |
| Helper                | 1.019 |
| Synchrony             | 2.358 |
| Breeding pairs        | 1.840 |
| Hatch order:Synchrony | 2.258 |

``` r
write_csv(table_s1, paste0(output_file_path, "/Table_S1.csv"))
```

## Figure S1

Scaled residual plots of the global model predicting nestling sex.

``` r
simulateResiduals(fittedModel = global_model, n = 1000, plot = T)
```

![](GRANH1_R_Code_figures_for_markdown/Figure%20s1-1.png)<!-- -->

    ## Object of Class DHARMa with simulated residuals based on 1000 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
    ##  
    ## Scaled residual values: 0.6475772 0.05230089 0.3210632 0.2074911 0.5818266 0.6430925 0.2300968 0.1512932 0.588742 0.7298962 0.5108008 0.3400541 0.6101324 0.9822496 0.7006727 0.03414776 0.2935815 0.2436966 0.2279934 0.5502256 ...

``` r
# Saving this figure separately as Figure s1
figure_s1 <- simulateResiduals(fittedModel = global_model, n = 1000)
pdf(paste0(output_file_path, "/Figure_S1.pdf"))
plot(figure_s1)
dev.off()
```

    ## png 
    ##   2

## Table S2

Table S2 To estimate the minimal detectable effect size for our modeling
approach, we use the simr package. Please note that, in an effort to
avoid the pitfalls of a post-hoc power analysis, we do not use our
empirical data to estimate the effect size parameter (i.e., we are not
estimating power using the effect size observed in our study).

The data used are taken from [Khwaja et
al. 2017](https://doi.org/10.1086/693532) and can be accessed
[here](https://datadryad.org/stash/dataset/doi:10.5061/dryad.9bk88).

``` r
Helper_model <- lmer(Sex ~ Helper_presence + (1 | Year/Location), dfcc)

# Define reasonable effect sizes from the literature 
khwaja_file_path <- dryad_download("10.5061/dryad.9bk88")[[1]][2]
khwaja_df <- read_tsv(khwaja_file_path)

khwaja_helped_df <- khwaja_df %>% 
  filter(predictor == "helped") %>% 
  arrange(effect_size) %>% 
  mutate(power = NA) %>% 
  group_by(authors, effect_size) %>% 
  mutate(Significant_predictor = case_when(is.na(low_95) | is.na(up_95) ~ NA_character_, 
                                           between(0, low_95, up_95) ~ "No", 
                                           TRUE ~ "Yes")) %>% 
  ungroup()

# Define conventional estimates of effect sizes as per Lipsey and Wilson 2001
lipsey_df <- tibble(authors = "Lipsey and Wilson", year = 2001, species = c("Small effect size", 
                                                                            "Medium effect size", 
                                                                            "Large effect size"), 
                    effect_size = c(0.2, 0.5, 0.8))


power_df <- bind_rows(khwaja_helped_df, lipsey_df) %>% 
  arrange(effect_size)


# Power analysis using effect size df
# Set seed for reproducibility 
set.seed(67)
for (i in 1:nrow(power_df)) {
  # Replace the original effect size (study data) with ones from literature
  fixef(Helper_model)['Helper_presence'] =  power_df$effect_size[i]
  # Power simulation 
  pwr.summary <- summary(powerSim(Helper_model, test = simr::fixed('Helper_presence', "t"),
                                  nsim = 100, progress = FALSE))
  
  power_df$power[i] <- as.numeric(pwr.summary)[3]
}

table_s2 <- power_df %>% 
  mutate(Study = paste0(authors, " (", year, ")")) %>% 
  select(Study, Species = species, Significant_predictor, Effect_size = effect_size, Power_estimate = power) %>% 
  mutate(Effect_size = round(Effect_size, 3))

table_s2 %>% 
  kable()
```

| Study                    | Species              | Significant_predictor | Effect_size | Power_estimate |
|:-------------------------|:---------------------|:----------------------|------------:|---------------:|
| Koenig et al. (2001)     | Acorn woodpecker     | No                    |       0.031 |           0.11 |
| Nam et al. (2011)        | Long-tailed tit      | No                    |       0.050 |           0.15 |
| Khwaja et al. (2016)     | Rifleman             | No                    |       0.079 |           0.35 |
| Gressler et al. (2014)   | White-banded tanager | No                    |       0.083 |           0.34 |
| Legge et al. (2001)      | Laughing kookaburra  | No                    |       0.172 |           0.86 |
| Lipsey and Wilson (2001) | Small effect size    | NA                    |       0.200 |           0.96 |
| Doutrelant et al. (2004) | Sociable weaver      | Yes                   |       0.338 |           1.00 |
| Komdeur et al. (1997)    | Seychelles warbler   | Yes                   |       0.364 |           1.00 |
| Lipsey and Wilson (2001) | Medium effect size   | NA                    |       0.500 |           1.00 |
| Lipsey and Wilson (2001) | Large effect size    | NA                    |       0.800 |           1.00 |

## Figure 2

Sex ratio of nestlings with respect to (a) whether or not a helper was
observed at their nest and (b) position in hatch order. In both, error
bars show 95% confidence intervals, the dotted line indicates a 50:50
sex ratio, and the number of nestlings sexed in each year is shown above
each error bar.

``` r
helper_plot <- df %>%
  filter(!is.na(Helper_presence)) %>% 
  group_by(Helper_presence) %>%
  summarize(Mean_sex = mean(Sex), 
            n = n(), 
            se = sd(Sex)/sqrt(n)) %>%
  ggplot(aes(x = Helper_presence, y = Mean_sex)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = Mean_sex - 1.96*se, ymax = Mean_sex + 1.96*se), width = 0.75) + 
  geom_text(aes(label = paste(n), y = Mean_sex + 1.96*se + 0.10), color = "black") + 
  theme_simple + 
  ylab("Proportion of males in nest") + 
  xlab("Helper Presence") + 
  scale_x_continuous(breaks = seq(0, 1, 1), labels = c("No helpers", "Helpers")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) + 
  geom_hline(aes(yintercept = 0.5), color = "black", lty = 2) + 
  ggtitle("a.") 

hatch_fraction_plot <- df %>%
  filter(!is.na(Hatch_order)) %>% 
  group_by(Hatch_order) %>%
  summarize(Mean_sex = mean(Sex), 
            n = n(), 
            se = sd(Sex)/sqrt(n)) %>%
  ungroup() %>%
  mutate(Hatch_order = factor(Hatch_order, levels = c("First", "Middle", "Last"))) %>%
  ggplot(aes(x = Hatch_order, y = Mean_sex)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = Mean_sex - 1.96*se, ymax = Mean_sex + 1.96*se), width = 0.75) + 
  geom_text(aes(label = paste(n), y = Mean_sex + 1.96*se + 0.10), color = "black") + 
  theme_simple + 
  ylab("") + 
  xlab("Position in hatch order") + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) + 
  geom_hline(aes(yintercept = 0.5), color = "black", lty = 2) + 
  ggtitle("b.")

figure_2 <- grid.arrange(helper_plot, hatch_fraction_plot, nrow = 1, respect = TRUE)
```

![](GRANH1_R_Code_figures_for_markdown/Figure%202-1.png)<!-- -->

``` r
ggsave(paste0(output_file_path, "/Figure_2.pdf"), figure_2, width = 7, height = 3)
```

## Table S2

Please not that the formatting of this table has not been reproduced
below. Instead, I present the information as 4 separate tables. The
caption for this table presented in the Supplementary Information is:
“Parameters and random effect estimates for (a) the best-fit and (b)
global logistic mixed effects models predicting the sex of individual
ani nestlings. Significant parameters (p\<0.05) have been bolded. The
reference levels for the categorical fixed effects are: first hatching
(for hatch order) and two-pair group (for breeding pairs).”

``` r
# AICc for top model presented in table 
AICc(top_model) %>% round(1)
```

    ## [1] 656.3

``` r
# Best-fit model table S2
top_short_table <- top_model %>% 
  tidy() %>%
  mutate(term = case_when(group == "Location:Year" ~ "sd_(Intercept)|Location:Year", 
                          group == "Year" ~ "sd_(Intercept)|Year",
                          is.na(group) ~ term)) %>%
  select(- c("group", "effect"))

set.seed(42) # for reproducibility 
top_cis <- lme4::confint.merMod(top_model, oldNames = FALSE, method = "profile") 

top_cis_df <- top_cis %>%
  as.data.frame() %>%
  mutate(term = rownames(.)) %>%
  mutate(CI_lower = `2.5 %`, 
         CI_upper = `97.5 %`) %>%
  select(term, CI_lower, CI_upper)

rownames(top_cis_df) <- NULL

table_s2_best_fit <- top_short_table %>%
  left_join(top_cis_df, by = "term") %>%
  mutate_if(is.numeric, function(x){round(x, 3)}) %>% 
  mutate(term = term_cleaner(term))

table_s2_best_fit %>% 
  kable()
```

| term                          | estimate | std.error | statistic | p.value | CI_lower | CI_upper |
|:------------------------------|---------:|----------:|----------:|--------:|---------:|---------:|
| (Intercept)                   |    0.123 |     0.157 |     0.781 |   0.435 |   -0.185 |    0.435 |
| Middle nestling               |   -0.352 |     0.218 |    -1.615 |   0.106 |   -0.785 |    0.072 |
| Last nestling                 |   -0.474 |     0.264 |    -1.793 |   0.073 |   -0.999 |    0.040 |
| Clutch size                   |    0.180 |     0.107 |     1.686 |   0.092 |   -0.034 |    0.395 |
| sd (Intercept)\|Location:Year |    0.390 |        NA |        NA |      NA |    0.000 |    0.739 |
| sd (Intercept)\|Year          |    0.000 |        NA |        NA |      NA |    0.000 |    0.302 |

``` r
write_csv(table_s2_best_fit, paste0(output_file_path, "/Table_S2_best_fit.csv"))
```

``` r
# AICc for global model presented in table  
AICc(global_model) %>% round(1)
```

    ## [1] 665.6

``` r
global_short_table <- global_model %>% 
  tidy() %>%
  mutate(term = case_when(group == "Location:Year" ~ "sd_(Intercept)|Location:Year", 
                          group == "Year" ~ "sd_(Intercept)|Year",
                          is.na(group) ~ term)) %>%
  select(- c("group", "effect"))


set.seed(42) # for reproducibility 
global_cis <- lme4::confint.merMod(global_model, oldNames = FALSE, method = "profile") 

global_cis_df <- global_cis %>%
  as.data.frame() %>%
  mutate(term = rownames(.)) %>%
  mutate(CI_lower = `2.5 %`, 
         CI_upper = `97.5 %`) %>%
  select(term, CI_lower, CI_upper)

rownames(global_cis_df) <- NULL


table_s2_global <- global_short_table %>%
  left_join(global_cis_df, by = "term") %>%
  mutate_if(is.numeric, function(x){round(x, 3)}) %>% 
  mutate(term = term_cleaner(term))

table_s2_global %>% 
  kable()
```

| term                          | estimate | std.error | statistic | p.value | CI_lower | CI_upper |
|:------------------------------|---------:|----------:|----------:|--------:|---------:|---------:|
| (Intercept)                   |    0.008 |     0.188 |     0.040 |   0.968 |   -0.364 |    0.382 |
| Middle nestling               |   -0.372 |     0.221 |    -1.681 |   0.093 |   -0.811 |    0.058 |
| Last nestling                 |   -0.481 |     0.269 |    -1.787 |   0.074 |   -1.016 |    0.043 |
| Brood size                    |    0.048 |     0.164 |     0.296 |   0.767 |   -0.277 |    0.375 |
| Clutch size                   |    0.110 |     0.180 |     0.611 |   0.541 |   -0.249 |    0.469 |
| Helper                        |    0.195 |     0.219 |     0.892 |   0.372 |   -0.242 |    0.635 |
| Synchrony                     |   -0.208 |     0.158 |    -1.321 |   0.186 |   -0.527 |    0.100 |
| Three-pair group              |    0.174 |     0.310 |     0.560 |   0.576 |   -0.440 |    0.801 |
| Middle nestling:Synchrony     |    0.287 |     0.220 |     1.308 |   0.191 |   -0.145 |    0.722 |
| Last nestling:Synchrony       |    0.113 |     0.270 |     0.418 |   0.676 |   -0.426 |    0.644 |
| sd (Intercept)\|Location:Year |    0.335 |        NA |        NA |      NA |    0.000 |    0.698 |
| sd (Intercept)\|Year          |    0.000 |        NA |        NA |      NA |    0.000 |    0.293 |

``` r
write_csv(table_s2_global, paste0(output_file_path, "/Table_S2_global.csv"))
```

# Assessing model residuals using DHARMa

Here, we use the “DHARMa” package to visualize the global model’s
residuals.

``` r
model_interrogator <- function(model, fixed_effects, df){
  
  # Running simulation
  sim <- simulateResiduals(fittedModel = model, plot = T, n = 1000)
  
  # General histogram looking for uniform distribution
  hist(sim)
  
  lapply(as.list(str_split(fixed_effects, pattern = stringr::fixed(" + "), 
                           simplify = TRUE)), 
         function(x){
           predictor <- unlist(df[, x])
           
           if(str_detect(x, stringr::fixed("scaled"))){
             predictor <- round(predictor, 3)
           }
           
           boxplot(sim$scaledResiduals ~ predictor, 
                   xlab = x, 
                   ylab = "Scaled Residuals")
         })
  
  
  # Goodness-of-fit tests
  testResiduals(sim)
  testZeroInflation(sim)
  testTemporalAutocorrelation(sim, df$unique_relative_time)
}
```

``` r
model_interrogator(global_model, "Hatch_order + Brood_size_scaled + Clutch_size_scaled + Helper_presence + Synchrony_scaled + Breeding_pairs", dfcc)
```

![](GRANH1_R_Code_figures_for_markdown/Global%20model%20DHARMa-1.png)<!-- -->![](GRANH1_R_Code_figures_for_markdown/Global%20model%20DHARMa-2.png)<!-- -->![](GRANH1_R_Code_figures_for_markdown/Global%20model%20DHARMa-3.png)<!-- -->![](GRANH1_R_Code_figures_for_markdown/Global%20model%20DHARMa-4.png)<!-- -->![](GRANH1_R_Code_figures_for_markdown/Global%20model%20DHARMa-5.png)<!-- -->![](GRANH1_R_Code_figures_for_markdown/Global%20model%20DHARMa-6.png)<!-- -->![](GRANH1_R_Code_figures_for_markdown/Global%20model%20DHARMa-7.png)<!-- -->![](GRANH1_R_Code_figures_for_markdown/Global%20model%20DHARMa-8.png)<!-- -->![](GRANH1_R_Code_figures_for_markdown/Global%20model%20DHARMa-9.png)<!-- -->![](GRANH1_R_Code_figures_for_markdown/Global%20model%20DHARMa-10.png)<!-- -->

    ## $uniformity
    ## 
    ##  Asymptotic one-sample Kolmogorov-Smirnov test
    ## 
    ## data:  simulationOutput$scaledResiduals
    ## D = 0.035441, p-value = 0.5937
    ## alternative hypothesis: two-sided
    ## 
    ## 
    ## $dispersion
    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 1.0009, p-value = 0.934
    ## alternative hypothesis: two.sided
    ## 
    ## 
    ## $outliers
    ## 
    ##  DHARMa bootstrapped outlier test
    ## 
    ## data:  simulationOutput
    ## outliers at both margin(s) = 0, observations = 472, p-value = 1
    ## alternative hypothesis: two.sided
    ##  percent confidence interval:
    ##  0 0
    ## sample estimates:
    ## outlier frequency (expected: 0 ) 
    ##                                0

![](GRANH1_R_Code_figures_for_markdown/Global%20model%20DHARMa-11.png)<!-- -->![](GRANH1_R_Code_figures_for_markdown/Global%20model%20DHARMa-12.png)<!-- -->

    ## 
    ##  Durbin-Watson test
    ## 
    ## data:  simulationOutput$scaledResiduals ~ 1
    ## DW = 2.0061, p-value = 0.9474
    ## alternative hypothesis: true autocorrelation is not 0
