
# Libs --------------------------------------------------------------------
library(readr)
library(dplyr)
library(lcmm)
library(flexmix)

# Setup -------------------------------------------------------------------
# Function to assign clusters from flexmix per-observation input/output to
# per subject level
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Load --------------------------------------------------------------------
data <- read_csv('data_s2m_long_total.csv') %>% 
  mutate(
    id = as.numeric(gsub('S2D', '1', gsub('S2M', '2', id))),
    GLP1 = glucose,
    time = t
  ) %>% 
  group_by(id) %>%
  filter(!any(is.na(GLP1) | is.na(time))) %>%
  ungroup() %>% 
  data.frame()


# Step 1 ------------------------------------------------------------------
model1 <- hlme(
  GLP1 ~ time + I(time^2) + I(time^3),
  #random = ~ 1,
  subject = 'id',
  data = data,
  ng = 1
)

max_classes <- 6
models <- list()
BICs <- c(model1$BIC)

for (k in 2:6) {
  ng_k <- k
  modelk <- gridsearch(
    rep = 50,
    maxiter = 100,
    minit = model1,
    hlme(
      GLP1 ~ time + I(time^2) + I(time^3),
      mixture = GLP1 ~time + I(time^2) + I(time^3),
      #random = ~ 1,
      subject = 'id',
      ng = ng_k,
      data = data
    )
  )
  models[[k]] <- modelk
  BICs <- c(BICs, modelk$BIC)
}

save.image('hlme.Rdata')

# Plot BIC
plot(
  1:max_classes,
  BICs,
  type = 'b',
  xlab = 'Number of Classes',
  ylab = 'BIC'
)

# Choose model based on BIC
optimal_ng <- which.min(BICs)
model_final <- if (optimal_ng == 1) model1 else models[[optimal_ng]]

# Get probabilities
pprob <- model_final$pprob

# Merge data with posterior probabilities from hlme
data_with_probs <- merge(data, pprob, by = "id", all.x = TRUE)

# Extract columns with probabilities
prob_cols <- grep("^prob", names(data_with_probs))
prob_matrix <- as.matrix(data_with_probs[, prob_cols])

# Run flexmix mixture
flexmix_model <- flexmix(
  GLP1 ~ time + I(time^2) + I(time^3) | id,
  model = list(
    FLXMRglm(GLP1 ~ time + I(time^2) + I(time^3), family = 'gaussian')
  ),
  data = data_with_probs,
  cluster = prob_matrix,  # Start with hlme's class probabilities
  control = list(
    iter.max = 1,         # Stop after the first iteration
    classify = "auto",    # Base classificator
    minprior = 0,         # Prevent clusters from being dropped
    tolerance = 1e10      # High tolerance to prevent further optimization
  )
)

# Comparison matrix -------------------------------------------------------
# Extract class assignments from flexmix per observation
data_with_probs$flexmix_class <- clusters(flexmix_model)

# Aggregate flexmix assignments to per-subject level
flexmix_classes_per_id <- data_with_probs %>%
  group_by(id) %>%
  summarize(class = Mode(flexmix_class))

# Prepare hlme class assignments
hlme_classes_per_id <- data.frame(id = pprob$id, class = pprob$class)

# Merge and compare
class_comparison <- merge(
  hlme_classes_per_id,
  flexmix_classes_per_id,
  by = 'id',
  suffixes = c('_hlme', '_flexmix')
)

# Create comparison table
comparison_table <- table(
  hlme = class_comparison$class_hlme,
  flexmix = class_comparison$class_flexmix
)

print(comparison_table)
