# Based on some previous literatures, we could use 3 strategy to build the score now
# 1. Covariates: Time-varying BMI, HBP, HBC, family history of DB, physical activity, energy intake, smoking status, SES, marital status
# 2. Based on Frank's suggestion, excluding all of the alcohol components and adj alcohol intake.
rm(list = ls())
gc()
options(bitmapType='cairo')

setwd("/udd/n2gji/micro/dietpt/final/data/")

library(glmnet)
library(survival)
library(dplyr)
library(haven)
library(patchwork)

# install.packages("table1")
# install.packages("pammtools")
# install.packages("pheatmap")
# install.packages("tidyverse")

library(table1)

# remotes::install_version("rms", version = "6.8-2", repos = "https://cloud.r-project.org")
# install.packages("pammtools")

library(pammtools)
library(tidyverse)
library(ggplot2)

# install.packages("pheatmap")

library(pheatmap)
library(RColorBrewer)
library(broom)

# library(chanmetab)
# library(patchwork)
# library(forcats)
# library(Biobase)
# library(tidyr)
# library(gridExtra)
# library(data.table)
# library(cowplot)
# library(purrr)
# library(ggrepel)
# library(stringr)
# library(MASS)
# library(statar)

nhs <- read_sas("nhs.sas7bdat")

fd <- c(
  "procm_serv","rmeat_serv","liver_serv","seaf_serv",
  "sseaf_serv","poul_serv","eggs_serv","butter_serv",
  "marg_serv","milk_serv","yogurt_serv","cheese_serv",
  "daiswt_serv","tea_serv","coff_serv","sugbev_serv",
  "lowbev_serv","citf_serv","fruj_serv","berr_serv",
  "appear_serv","ban_serv","rais_serv","peach_serv",
  "melon_serv","avo_serv","cel_serv","tom_serv",
  "squash_serv","lfveg_serv","cruveg_serv",
  "stachveg_serv","nightveg_serv","legu_serv",
  "mixveg_serv","othveg_serv","tompro_serv",
  "apiveg_serv","fries_serv","wgrain_serv",
  "rgrain_serv","snack_serv","pnuts_serv",
  "otnuts_serv","choco_serv","swt_serv","brdo_serv",
  "swtr_serv","pie_serv","cokh_serv", "cokr_serv",
  "cakh_serv", "cakr_serv","cond_serv",
  "saldre_serv","crsoup_serv","pizza_serv"
)
nhs <- nhs %>%
  group_by(interval) %>%
  mutate(
    across(
      .cols = all_of(fd),
      .fns = ~ as.numeric(scale(log(.x + 1))),
      .names = "{.col}_log_z"
    )
  ) %>%
  ungroup()

food_vars_log_z <- paste0(fd, "_log_z")

nhs <- nhs %>%
  group_by(id) %>%
  mutate(bmibase = bmicon[interval == 1][1]) %>%
  mutate(actbase = actcon[interval == 1][1]) %>%
  mutate(alcobase = alcocon[interval == 1][1]) %>%
  ungroup()
nhs <- nhs %>%
  mutate(
    bmi_cat = cut(
      bmicon,
      breaks = c(-Inf, 21.0, 25.0, 30.0, 32.0, Inf),
      right = FALSE,
      labels = c("u21", "21_25.9", "25_30", "30_32", "o32")
    ),
    
    bmi_u21    = as.numeric(bmi_cat == "u21"),
    bmi_25_30  = as.numeric(bmi_cat == "25_30"),
    bmi_30_32  = as.numeric(bmi_cat == "30_32"),
    bmi_o32    = as.numeric(bmi_cat == "o32")
  )

table(nhs$bmi_cat, useNA = "always")

nhs <- nhs %>%
  group_by(interval) %>%
  mutate(
    act_cat = ntile(actcon, 5),
    act_q2 = as.numeric(act_cat == 2),
    act_q3 = as.numeric(act_cat == 3),
    act_q4 = as.numeric(act_cat == 4),
    act_q5 = as.numeric(act_cat == 5)
  ) %>%
  ungroup()

table(nhs$act_cat, useNA = "always")

covariates_tdc <- c(
  food_vars_log_z,
  "bmibase", "actbase", "bmicon", "actcon",
  "alcobase", "alcocon", "calorscon",
  "smkstatus2", "smkstatus3", "interval",
  "marry", "ses", "dbfh", "basehighbp",
  "basehightc", "highbp", "hightc",
  "alc2", "alc3", "alc4", "alc5", "alc6", "mv"
)

baseline_age_var <- "age86"

df_base1 <- nhs %>%
  group_by(id) %>%
  arrange(id, interval) %>%
  dplyr::summarise(
    entry_time_calendar = first(irt86),
    event_time_calendar = first(overall_tend),
    entry_age = first(.data[[baseline_age_var]]),
    event_status = max(db2, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    follow_up_years = as.numeric(event_time_calendar - entry_time_calendar) / 12,
    event_age = entry_age + follow_up_years
  ) %>%
  filter(event_age > entry_age)

baseline_info_df <- nhs %>%
  filter(interval == 1) %>%
  dplyr::select(
    id,
    entry_time_calendar = irt86,
    entry_age = .data[[baseline_age_var]]
  ) %>%
  distinct(id, .keep_all = TRUE)

df_long_covariates1 <- nhs %>%
  left_join(baseline_info_df, by = "id") %>%
  mutate(
    visit_time_calendar = case_when(
      interval == 1  ~ irt86,
      interval == 2  ~ irt88,
      interval == 3  ~ irt90,
      interval == 4  ~ irt92,
      interval == 5  ~ irt94,
      interval == 6  ~ irt96,
      interval == 7  ~ irt98,
      interval == 8  ~ irt00,
      interval == 9  ~ irt02,
      interval == 10 ~ irt04,
      interval == 11 ~ irt06,
      interval == 12 ~ irt08,
      interval == 13 ~ irt10,
      TRUE ~ NA_real_
    ),
    
    visit_age = entry_age +
      as.numeric(visit_time_calendar - entry_time_calendar) / 12
  ) %>%
  dplyr::select(id, visit_age, all_of(covariates_tdc)) %>%
  filter(!is.na(visit_age))

df_cox_final1 <- tmerge(
  data1 = df_base1,
  data2 = df_base1,
  id = id,
  tstart = entry_age,
  tstop = event_age,
  event = event(event_age, event_status)
)

for (var in covariates_tdc) {
  
  args_list <- list(
    data1 = df_cox_final1,
    data2 = df_long_covariates1,
    id = quote(id)
  )
  
  tdc_call <- call("tdc", quote(visit_age), as.name(var))
  
  args_list[[var]] <- tdc_call
  
  df_cox_final1 <- do.call(tmerge, args_list)
  
  cat("Coming to:", var, "\n")
}

# save(df_cox_final1,file = "final_cox_data.RData")

df_cox_final1 <- df_cox_final1 %>%
  mutate(
    bmi_cat = cut(
      bmicon,
      breaks = c(-Inf, 21.0, 25.0, 30.0, 32.0, Inf),
      right = FALSE,
      labels = c("u21", "21_25.9", "25_30", "30_32", "o32")
    ),
    
    bmi_u21   = as.numeric(bmi_cat == "u21"),
    bmi_25_30 = as.numeric(bmi_cat == "25_30"),
    bmi_30_32 = as.numeric(bmi_cat == "30_32"),
    bmi_o32   = as.numeric(bmi_cat == "o32")
  )

table(df_cox_final1$bmi_cat, useNA = "always")
table(df_cox_final1$bmi_u21, useNA = "always")
table(df_cox_final1$bmi_25_30, useNA = "always")
table(df_cox_final1$bmi_30_32, useNA = "always")
table(df_cox_final1$bmi_o32, useNA = "always")

bmi_dummy_vars <- c(
  "bmi_u21",
  "bmi_25_30",
  "bmi_30_32",
  "bmi_o32"
)

df_cox_final1 <- df_cox_final1 %>%
  group_by(interval) %>%
  mutate(
    act_cat = ntile(actcon, 5),
    act_q2 = as.numeric(act_cat == 2),
    act_q3 = as.numeric(act_cat == 3),
    act_q4 = as.numeric(act_cat == 4),
    act_q5 = as.numeric(act_cat == 5)
  ) %>%
  ungroup()

act_dummy_vars <- c("act_q2", "act_q3", "act_q4", "act_q5")

df_cox_final1$interval <- factor(df_cox_final1$interval, levels = 1:13)

covars_to_adjust <- c(
  act_dummy_vars,
  "bmicon", "calorscon",
  "smkstatus2", "smkstatus3",
  "highbp", "hightc",
  "ses", "marry", "dbfh",
  "alcocon", "mv"
)

X_vars <- c(food_vars_log_z, covars_to_adjust)

X_vars

Y_base <- Surv(
  time = df_cox_final1$tstart,
  time2 = df_cox_final1$tstop,
  event = df_cox_final1$event
)

Y <- stratifySurv(Y_base, strata = df_cox_final1$interval)

X <- as.matrix(df_cox_final1[, X_vars])

n_food <- length(food_vars_log_z)
n_confounder <- length(covars_to_adjust)

penalty_vec <- c(rep(0.5, n_food), rep(0, n_confounder))
alpha_seq <- seq(0.1, 0.9, by = 0.1)
cv_fit_list <- list()
n_folds <- 10
#Avoid data leakage
id_col <- df_cox_final1$id 
unique_ids <- unique(id_col)
n_unique_ids <- length(unique_ids)
set.seed(456)
id_fold_assignment <- sample(rep(1:n_folds, length.out = n_unique_ids))
names(id_fold_assignment) <- as.character(unique_ids)
fold_id <- id_fold_assignment[as.character(id_col)]
fold_id <- unname(fold_id) 
#fold_id <- sample(1:n_folds, size = nrow(X), replace = TRUE)

rm(nhs, df_base1, df_long_covariates1)

for (a in alpha_seq) {
  cat("Running alpha =", a, "...\n")
  fit <- cv.glmnet(
    x = X,
    y = Y,
    family = "cox",
    trace.it = FALSE,
    alpha = a,
    standardize = FALSE,
    type.measure = "C",
    parallel = FALSE,
    penalty.factor = penalty_vec,
    foldid = fold_id
  )
  cv_fit_list[[as.character(a)]] <- fit
}
save(
  cv_fit_list,
  file = "final_alphas_01_to_09_withoutalc_new.RData"
)
