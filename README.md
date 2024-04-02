# Reproduction scripts for: Semi-parametric benchmark dose analysis with monotone additive models

This repository contains scripts for reproducing the analyses in the submitted paper named above.
First, download and install the `semibmd` `R` package:
```
remotes::install_github("awstringer1/semibmd")
```
The following scripts reproduce the results in the paper:
1. `01-simulations.R`: the simulation results of Section 4.
2. `02-pae-data.R`: the PAE data analysis of Section 5.

To run the scripts, you need to install the `semibmd` package above or download it to your computer.
An informative error will be thrown if this has not been done.
You then need to create a directory to store the results, and then change the corresponding path
in the script. Other than that, they should run as-is, with no other user changes required.
They will attempt to install all other required packages from `CRAN` automatically.

To run `02-pae-data.R`, you need the Prenatal Alcohol Exposure (PAE) data.
This is obtained from the [Harvard DataVerse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/SCLH9W);
the link is also given in the code file itself. All you have to do is convert the `.xlsx` file to a `.csv` and then place the data
on your computer as described in the script.
