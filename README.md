# Censored Survival Modeling with Regression Trees

Implementation and evaluation of a median regression tree for censored time-to-event data, benchmarked against a parametric Weibull regression model.

## Overview

This project compares two approaches for modeling censored survival data:

- **Median Regression Tree** — captures nonlinear relationships and subgroup effects through residual-based recursive partitioning
- **Weibull Regression** — provides a parametric baseline for survival-time prediction

The analysis evaluates both predictive accuracy and survival ranking performance.

## Dataset

The project uses the **rat tumor survival dataset** available through the R `survival` package.

Variables used in the analysis include:

- `time` — observed survival or censoring time
- `status` — event indicator
- `rx` — treatment group
- `sex` — sex
- `litter` — litter identifier

The data can be accessed directly in R through the `survival` package, so no external dataset is required.

## Methods

- Kaplan–Meier survival estimation
- Censoring-adjusted residuals
- Residual-based split selection
- Recursive partitioning
- Median regression tree
- Weibull accelerated failure time regression
- 70/30 train-test split
- Mean Absolute Error (MAE)
- Concordance Index (C-index)

## Model Evaluation

The regression tree and Weibull model are compared using:

- **MAE** for prediction accuracy
- **C-index** for ranking survival outcomes
- Separate evaluation of all observations and uncensored observations

The analysis also examines how censoring and tied tree predictions affect model performance.

## Technologies

R · Survival Analysis · Statistical Learning · Kaplan–Meier Estimation

## Repository Contents

- `survival_tree.R` — complete R analysis
- `survival_tree_report.pdf` — full project report

## Running the Project

Install the required packages:

```r
install.packages(c("survival", "ggplot2"))
