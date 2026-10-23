# MissForest: reproduction of the categorical data experiment

Reproduction of the categorical data experiment from Stekhoven & Bühlmann (2012), [*MissForest – non-parametric missing value imputation for mixed-type data*](https://arxiv.org/abs/1105.0828), made for Assignment 2 of Data Wrangling (INFOMDW).

## Background
MissForest imputes missing values with random forests. It starts with a simple guess (the most frequent category for categorical variables). Then, for each variable, it fits a random forest on the observed values and predicts the missing ones. It repeats this until the imputed data stop improving. It needs no tuning parameters, and its out-of-bag (OOB) error gives an estimate of the imputation error without a test set.

## Experiment
For each dataset, 10%, 20% and 30% of the values are removed completely at random and imputed with missForest, with 50 runs per setting. Each run records:
- the **true PFC**: the proportion of imputed values that don't match the removed ones
- missForest's **OOB estimate** of the PFC
- the runtime

The results are compared with Figure 2, Figure 4 and Table 1 of the paper. Only missForest is run, not the comparison methods.

## Data
All three datasets come from the UCI Machine Learning Repository:
- [SPECT Heart](https://archive.ics.uci.edu/dataset/95/spect+heart): 267 × 22
- [Promoter Gene Sequences](https://archive.ics.uci.edu/dataset/67/molecular+biology+promoter+gene+sequences): 106 × 57
- [Lymphography](https://archive.ics.uci.edu/dataset/63/lymphography): 148 × 19

Put the downloaded files in a `data/` folder: `SPECT.train`, `SPECT.test`, `promoters.data`, `lymphography.data`.

## How to run
1. Download the data into `data/`.
2. Open `missforest_categorical.R` in RStudio and set the working directory to the script's folder (**Session → Set Working Directory → To Source File Location**).
3. Click **Source**. Missing packages (`missForest`, `randomForest`) are installed automatically.

Set `QUICK <- TRUE` for a quick test run (3 simulations instead of 50).

**Main parameters:** `ntree = 100`, `mtry = ⌊√p⌋`, `maxiter = 10`, the `randomForest` backend as in the paper, and seed 2012.

**Output:** `results/raw_results.csv`, a summary table in the console, and the reproduced figures.

`results/raw_results.csv` contains one row per simulation run (450 rows in total). Each row gives:
- the dataset and the missingness rate
- the simulation number
- the true PFC and missForest's OOB estimate of it
- the imputation time in seconds
- how many times the missingness pattern had to be redrawn because of a known missForest bug

## Results
The reproduced PFC values and OOB differences closely match the paper. PFC increases with the missingness rate, and the OOB estimate underestimates the true error, more so at higher missingness rates.
