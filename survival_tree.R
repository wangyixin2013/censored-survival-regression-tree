# =============================================================================
# Censored Survival Modeling with Regression Trees
# Original analysis code — formatting cleaned, logic preserved
#
# Notes:
# - No original code lines were removed or reordered.
# - Only pasted Markdown/Unicode artifacts were normalized.
# - Additional comments were added for readability.
# =============================================================================
#This is all R-code needed to for the report, purpose of R-code in each section is labeled.
#Run through the R code by sections to build models and test models.
library(survival)
# -----------------------------------------------------------------------------
# SECTION 1 — Data setup and preprocessing
# -----------------------------------------------------------------------------
# Setting up data
data(cancer, package = "survival")
ls()
colnames(rats)
# Keep variables needed
rats2 <- rats[, c("time", "status", "rx", "sex", "litter")]
# Convert types
rats2$rx  <- factor(rats2$rx, levels = c(0, 1), labels = c("Control", "Drug"))
rats2$sex <- factor(rats2$sex)
# Data without litter
rats_simple <- rats2[, c("time", "status", "rx", "sex","litter")]
# Check
head(rats_simple)
str(rats_simple)
summary(rats_simple)
# -----------------------------------------------------------------------------
# SECTION 2 — Median regression tree construction
# -----------------------------------------------------------------------------
#building the median regression tree model
# Step 1: find median
# create survival object
SurvObj <- Surv(rats_simple$time, rats_simple$status)
# Kaplan-Meier fit
km_fit <- survfit(SurvObj ~ 1)
summary(km_fit)
# find median survival time
theta_hat <- km_fit$time[which.min(abs(km_fit$surv - 0.5))]
theta_hat
# Step 2: compute residuals e_i
# To estimate G(t), the censoring survival function,
# treat censoring as the "event":
#   original status: 1 = tumor, 0 = censored
#   censor_event:    1 = censored, 0 = tumor
censor_event <- 1 - rats_simple$status
# KMRAT estimate of censoring survival function G(t)
G_fit <- survfit(Surv(rats_simple$time, censor_event) ~ 1)
# evaluate KM survival estimate at a chosen time t
get_step_surv <- function(fit, t) {
  s <- summary(fit, times = t, extend = TRUE)$surv
  if (length(s) == 0 || is.na(s)) return(1)
  s
}
# Ghat(theta_hat)
G_theta <- get_step_surv(G_fit, theta_hat)
if (G_theta <= 0) stop("")
# Paper residuals
rats_simple$ei <- as.numeric(rats_simple$time >= theta_hat) / G_theta - 0.5
# Inspect
theta_hat
G_theta
head(rats_simple[, c("time", "status", "rx", "sex", "ei")])
summary(rats_simple$ei)
table(rats_simple$time >= theta_hat)
# Step 3: choose split variable using residual analysis
# residual sign
rats_simple$resid_sign <- ifelse(rats_simple$ei > 0, "pos", "neg")
rats_simple$resid_sign <- factor(rats_simple$resid_sign, levels = c("neg", "pos"))
# rx
tab_rx <- table(rats_simple$resid_sign, rats_simple$rx)
p_rx <- fisher.test(tab_rx)$p.value
# sex
tab_sex <- table(rats_simple$resid_sign, rats_simple$sex)
p_sex <- fisher.test(tab_sex)$p.value
# litter
# group litter into quartiles
q <- quantile(rats$litter, probs = seq(0, 1, 0.25), na.rm = TRUE)
# if duplicate breakpoints happen, keep unique ones
q <- unique(q)
rats$litter_group <- cut(rats$litter, breaks = q, include.lowest = TRUE)
tab_litter <- table(rats_simple$resid_sign, rats$litter_group)
# if expected counts are small, use simulated p-value
p_litter <- chisq.test(tab_litter, simulate.p.value = TRUE)$p.value
# compare all p-values
pvals <- c(rx = p_rx, sex = p_sex, litter = p_litter)
pvals
# choose split variable
split_var <- names(which.min(pvals))
split_var
#we get that the first split variable is sex
# Step 4: split by sex
node_f <- subset(rats_simple, sex == "f")
node_m <- subset(rats_simple, sex == "m")
# female node median
km_f <- survfit(Surv(time, status) ~ 1, data = node_f)
theta_f <- summary(km_f)$table["median"]
if (is.na(theta_f)) {
  theta_f <- max(km_f$time)
} else {
  theta_f <- as.numeric(theta_f)
}
# male node median
km_m <- survfit(Surv(time, status) ~ 1, data = node_m)
theta_m <- summary(km_m)$table["median"]
if (is.na(theta_m)) {
  theta_m <- max(km_m$time)
} else {
  theta_m <- as.numeric(theta_m)
}
# current tree structure:
#         All(split by sex)
#         /   \\
#   famale     male
#we work on female node first
# Female residual
node_f$ei <- ifelse(node_f$time >= theta_f, 0.5, -0.5)
head(node_f[, c("time", "status", "rx", "ei")])
table(node_f$ei)
# 1. residual sign
node_f$resid_sign <- ifelse(node_f$ei > 0, "pos", "neg")
# 2. rx test
tab_rx_f <- table(node_f$resid_sign, node_f$rx)
p_rx_f <- fisher.test(tab_rx_f)$p.value
# 3. litter test
# split into 4 groups
node_f$litter_group <- cut(node_f$litter,
                           breaks = quantile(node_f$litter, probs = seq(0,1,0.25)),
                           include.lowest = TRUE)
tab_litter_f <- table(node_f$resid_sign, node_f$litter_group)
# chisq test
p_litter_f <- chisq.test(tab_litter_f)$p.value
# 4. compare
p_rx_f
p_litter_f
# 5. choose variable
pvals <- c(rx = p_rx_f, litter = p_litter_f)
split_var_f <- names(which.min(pvals))
split_var_f
# we find the the next variable we split is rx
node_f_control <- subset(node_f, rx == "Control")
node_f_drug    <- subset(node_f, rx == "Drug")
nrow(node_f_control)
nrow(node_f_drug)
# current tree structure:
#           All(split by sex)
#         /                 \\
#   famale(split by rx)     male
#     /         \\
# control       drug
# female + control
km_f_control <- survfit(Surv(time, status) ~ 1, data = node_f_control)
theta_f_control <- km_f_control$time[which.min(abs(km_f_control$surv - 0.5))]
# female + drug
theta_f_drug <- summary(km_f_drug)$table["median"]
if (is.na(theta_f_drug)) {
  theta_f_drug <- 104
} else {
  theta_f_drug <- as.numeric(theta_f_drug)
}
theta_f_control
theta_f_drug
# we now go into control group of female to see if we need to split further
# residual
node_f_control$ei <- ifelse(node_f_control$time >= theta_f_control, 0.5, -0.5)
# residual sign
node_f_control$resid_sign <- ifelse(node_f_control$ei > 0, "pos", "neg")
# test litter
node_f_control$litter_group <- cut(
  node_f_control$litter,
  breaks = quantile(node_f_control$litter, probs = seq(0,1,0.25)),
  include.lowest = TRUE
)
tab_litter_fc <- table(node_f_control$resid_sign, node_f_control$litter_group)
p_litter_fc <- chisq.test(tab_litter_fc)$p.value
p_litter_fc
# p-value is greater than 0.05, we don't need to split further
# now we go into drug group of female to see if we need to split further
# residual
node_f_drug$ei <- ifelse(node_f_drug$time >= theta_f_drug, 0.5, -0.5)
# residual sign
node_f_drug$resid_sign <- ifelse(node_f_drug$ei > 0, "pos", "neg")
# test litter
node_f_drug$litter_group <- cut(
  node_f_drug$litter,
  breaks = quantile(node_f_drug$litter, probs = seq(0,1,0.25)),
  include.lowest = TRUE
)
tab_litter_fd <- table(node_f_drug$resid_sign, node_f_drug$litter_group)
p_litter_fd <- chisq.test(tab_litter_fd)$p.value
p_litter_fd
# p-value>0.05, so we don't need to split further
# now we go into the male branch
# compute residuals for male node
node_m$ei <- ifelse(node_m$time >= theta_m, 0.5, -0.5)
# residual sign
node_m$resid_sign <- ifelse(node_m$ei > 0, "pos", "neg")
node_m$resid_sign <- factor(node_m$resid_sign, levels = c("neg", "pos"))
# rx
tab_rx_m <- table(node_m$resid_sign, node_m$rx)
p_rx_m <- fisher.test(tab_rx_m)$p.value
# litter
# group into quartiles
q_m <- quantile(node_m$litter, probs = seq(0, 1, 0.25), na.rm = TRUE)
q_m <- unique(q_m)
node_m$litter_group <- cut(node_m$litter, breaks = q_m, include.lowest = TRUE)
tab_litter_m <- table(node_m$resid_sign, node_m$litter_group)
# use chi-square with simulation
p_litter_m <- chisq.test(tab_litter_m, simulate.p.value = TRUE)$p.value
# compare
p_rx_m
p_litter_m
pvals_m <- c(rx = p_rx_m, litter = p_litter_m)
pvals_m
# choose split variable
split_var_m <- names(which.min(pvals_m))
split_var_m
node_m_control <- subset(node_m, rx == "Control")
node_m_drug    <- subset(node_m, rx == "Drug")
nrow(node_m_control)
nrow(node_m_drug)
# male + control
km_m_control <- survfit(Surv(time, status) ~ 1, data = node_m_control)
# male + drug
km_m_drug <- survfit(Surv(time, status) ~ 1, data = node_m_drug)
get_km_pred <- function(km_fit) {
  idx <- which(km_fit$surv <= 0.5)[1]
  if (!is.na(idx)) {
    return(km_fit$time[idx])
  } else {
    return(max(km_fit$time))
  }
}
theta_m_control <- get_km_pred(km_m_control)
theta_m_drug <- get_km_pred(km_m_drug)
theta_m_control
theta_m_drug
# we go into the control group of male to see if we need to split further
theta_m_control <- as.numeric(quantile(km_m_control, probs = 0.5)$quantile)
if(is.na(theta_m_control)) {
  theta_m_control <- 104
  # median is undefined, so we don't need to split further
} else {
  # residual
  node_m_control$ei <- compute_node_residuals(node_m_control, theta_m_control)
  # residual sign
  node_m_control$resid_sign <- ifelse(node_m_control$ei > 0, "pos", "neg")
  # test litter
  q_mc <- unique(quantile(node_m_control$litter, probs = seq(0,1,0.25), na.rm = TRUE))
  node_m_control$litter_group <- cut(node_m_control$litter, breaks = q_mc, include.lowest = TRUE)
  tab_litter_mc <- table(node_m_control$resid_sign, node_m_control$litter_group)
  p_litter_mc <- chisq.test(tab_litter_mc, simulate.p.value = TRUE)$p.value
  p_litter_mc
  # p-value > 0.05, so we don't need to split further
}
# now we go into drug group of male to see if we need to split further
theta_m_drug <- as.numeric(quantile(km_m_drug, probs = 0.5)$quantile)
if(is.na(theta_m_drug)) {
  theta_m_drug <- 104
  # median is undefined, so we don't need to split further
} else {
  # residual
  node_m_drug$ei <- compute_node_residuals(node_m_drug, theta_m_drug)
  # residual sign
  node_m_drug$resid_sign <- ifelse(node_m_drug$ei > 0, "pos", "neg")
  # test litter
  q_md <- unique(quantile(node_m_drug$litter, probs = seq(0,1,0.25), na.rm = TRUE))
  node_m_drug$litter_group <- cut(node_m_drug$litter, breaks = q_md, include.lowest = TRUE)
  tab_litter_md <- table(node_m_drug$resid_sign, node_m_drug$litter_group)
  p_litter_md <- chisq.test(tab_litter_md, simulate.p.value = TRUE)$p.value
  p_litter_md
}
# our final tree structure
# current tree structure:
#             All(split by sex)
#            /                 \\
#   famale(split by rx)        male
#     /         \           /        \\
# control       drug    control     drug
# 102           103      104          104
# -----------------------------------------------------------------------------
# SECTION 3 — Train/test evaluation
# -----------------------------------------------------------------------------
#testing, spliting data into 70% training data and 30% testing data
set.seed(123)
n <- nrow(rats_simple)
train_idx <- sample(1:n, size = 0.7 * n)
train <- rats_simple[train_idx, ]
test  <- rats_simple[-train_idx, ]
# fit Weibull on training set
weib_model <- survreg(Surv(time, status) ~ rx + sex,
                      data = train, dist = "weibull")
# predict median survival time
train$pred_weibull <- predict(weib_model, newdata = train, type = "quantile", p = 0.5)
test$pred_weibull  <- predict(weib_model, newdata = test,  type = "quantile", p = 0.5)
# tree prediction using your fitted tree rules
predict_tree <- function(df) {
  ifelse(df$sex == "f" & df$rx == "Control", 102,
         ifelse(df$sex == "f" & df$rx == "Drug",    103,
                ifelse(df$sex == "m" & df$rx == "Control", 104, 104)))
}
train$pred_tree <- predict_tree(train)
test$pred_tree  <- predict_tree(test)
# error for all data
train_err_tree <- mean(abs(train$time - train$pred_tree))
train_err_weib <- mean(abs(train$time - train$pred_weibull))
test_err_tree <- mean(abs(test$time - test$pred_tree))
test_err_weib <- mean(abs(test$time - test$pred_weibull))
train_err_tree
train_err_weib
test_err_tree
test_err_weib
#error for uncensored data
train_err_tree_uncensored <- mean(abs(train$time[train$status == 1] -
                             train$pred_tree[train$status == 1]))
test_err_tree_uncensored <- mean(abs(test$time[test$status == 1] -
                            test$pred_tree[test$status == 1]))
train_err_weib_uncensored <- mean(abs(train$time[train$status == 1] -
                             train$pred_weibull[train$status == 1]))
test_err_weib_uncensored <- mean(abs(test$time[test$status == 1] -
                            test$pred_weibull[test$status == 1]))
train_err_tree_uncensored
train_err_weib_uncensored
test_err_tree_uncensored
test_err_weib_uncensored
# -----------------------------------------------------------------------------
# SECTION 4 — Visualization of prediction error
# -----------------------------------------------------------------------------
#present result in graph
library(ggplot2)
# ALL DATA
df_all <- data.frame(
  Dataset = rep(c("Train", "Test"), each = 2),
  Model = rep(c("Tree", "Weibull"), times = 2),
  MAE = c(13.60476, 101.9958,
          13.05556, 97.43306)
)
# UNCENSORED
df_unc <- data.frame(
  Dataset = rep(c("Train", "Test"), each = 2),
  Model = rep(c("Tree", "Weibull"), times = 2),
  MAE = c(26.25806, 48.78597,
          23.63636, 41.25857)
)
# GRAPH 1: ALL DATA
p1 <- ggplot(df_all, aes(x = Dataset, y = MAE, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Tree" = "steelblue", "Weibull" = "tomato")) +
  labs(title = "Prediction Error (All Data)",
       y = "Mean Absolute Error",
       x = "") +
  theme_minimal(base_size = 14)
print(p1)
# GRAPH 2: UNCENSORED
p2 <- ggplot(df_unc, aes(x = Dataset, y = MAE, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Tree" = "steelblue", "Weibull" = "tomato")) +
  labs(title = "Prediction Error (Uncensored Only)",
       y = "Mean Absolute Error",
       x = "") +
  theme_minimal(base_size = 14)
print(p2)
# -----------------------------------------------------------------------------
# SECTION 5 — Concordance-index evaluation
# -----------------------------------------------------------------------------
#C-index test
#test if predicted order is same as observed order
# Weibull
cindex_weib <- concordance(Surv(time, status) ~ pred_weibull, data = test)
# Tree
cindex_tree <- concordance(Surv(time, status) ~ pred_tree, data = test)
cindex_weib$concordance
cindex_tree$concordance
# -----------------------------------------------------------------------------
# SECTION 6 — Weibull model fit on the full dataset
# -----------------------------------------------------------------------------
#weibull model for all data
library(survival)
# whole data
weib_model_all <- survreg(Surv(time, status) ~ rx + sex + litter,
                          data = rats_simple,
                          dist = "weibull")
summary(weib_model_all)
