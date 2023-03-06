# README File for GRANH1 Data
The "GRANH1_anon_221123.csv" file is the data set used in the analyses in the manuscript "No evidence for adaptive sex ratio adjustment in a cooperatively breeding bird with helpful helpers." Each row represents a single Greater Ani (*Crotophaga major*) nestling, and some columns describe nestling-level characteristics while others describe brood-level characteristics (which are shared among several nestlings). In this dataset, no ani breeding group produced more than one brood of nestlings that were sampled in a single year. For this reason, combining the "Year" and "Location" columns (described below) creates a unique brood identifier. 

This data set will be made publicly available in both the `Data/` folder in this GitHub repository and on the Figshare repository upon publication. The DOI [doi.org/10.6084/m9.figshare.21760979](http://doi.org/10.6084/m9.figshare.21760979) has been reserved, and it will become active if the paper is accepted for publication. We apologize for any inconvenience.

## Columns Describing Nestling Characteristics:  
`Band_number` = This is a unique identifier of individual nestlings in our dataset<br>
`Sex` = The sex of a nestling; this column has only two possible values: "1" for males and "0" for females; the data have been filtered such that no nestlings of unknown sex are included<br> 
`Hatch_order` = The position in the hatching order of an ani nesting; this column has three possible values "First," "Middle," and "Last;" "First” and “Last” describe the absolute earliest and latest positions in the hatching order, respectively, and “Middle” describes all intermediate positions as in Riehl (2016)


## Columns Describing Brood Characteristics:
`Year` = The 4-digit year in which a nestling was sampled/in which an observation was made<br>
`Location` = An identifier of the location of a breeding group's territory; this variable is correlated with, but is not identical to, the name of the breeding group (i.e., parental identity) since breeding groups use the same territories across years but the membership of the groups can change<br>
`Sex_ratio` = The proportion of nestlings within a brood that are male; this takes numeric values between 0 (an all-female brood) and 1 (an all-male brood); this column was calculated by grouping the data by "Year" and "Location" (i.e., grouping by brood) and calculating the average of the "Sex" column.<br> 
`Synchrony` = The range of calendar days over which hatching occurred in a nest<br>
`Brood_size` = The number of eggs that hatched within a nest<br>
`Clutch_size` = The number of eggs that were incubated in a nest; this does not include eggs that were ejected from the nest<br>
`Helper_presence` = Binary variable describing whether or not a breeding group included a non-breeding helper; "1" values indicate that a nest had a helper; "0" values indicate that a nest did not have a helper<br>
`Breeding_pairs` =  The number of monogamous pairs of breeding anis that were observed at a nest while it was active and contributed to the communal clutch of eggs (either "Two" or "Three")
 
## Missing data
Missing data is represented throughout this dataset as "NA" 