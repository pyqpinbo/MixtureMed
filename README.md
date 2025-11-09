# MixtureMed <img src="https://img.shields.io/badge/version-0.1.0-blue.svg" alt="version"/> <img src="https://img.shields.io/badge/R-≥4.0.0-success.svg" alt="R version"/>

**MixtureMed** is an R package for *heterogeneous causal mediation analysis*
based on finite mixture models, allowing both parametric and tree-based
mixing proportions. It implements EM algorithms and nonparametric bootstrap
for inference on subgroup-specific natural indirect and direct effects (NIE/NDE).

---

## 🧭 Installation

You can install the development version of **MixtureMed** from GitHub:

```r
# install.packages("remotes")  # if not already installed
remotes::install_github("pyqpinbo/MixtureMed")


# Load the package
library(MixtureMed)


# Set random seed for reproducibility
set.seed(2025)

# Step 1: Simulate heterogeneous mediation data
sim <- data_generate(n = 300, setting = 2)

# Step 2: Fit the tree-based heterogeneous mediation model
fit <- Med_EM_tree(
  K          = 2,
  data_model = sim$data,
  compute_ci = TRUE,
  sims_ci    = 100,
  verbose    = FALSE
)

# Step 3: Display estimated subgroup-specific parameters
fit$coefs

# Step 4: Display estimated confidence intervals
fit$ci


