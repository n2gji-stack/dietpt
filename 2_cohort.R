rm(list = ls())
gc()
options(bitmapType='cairo')
setwd("/udd/n2gji/micro/dietpt/final/data/")
library(glmnet)
library(survival)
library(dplyr)
library(haven)
library(patchwork)
#install.packages("table1")
#install.packages("pammtools")
#install.packages("pheatmap")
#install.packages("tidyverse")
library(table1)
#remotes::install_version("rms", version = "6.8-2", repos = "https://cloud.r-project.org")
#install.packages("pammtools")
library(pammtools)
library(tidyverse) 
library(ggplot2)
#install.packages("pheatmap")
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
library(forestploter)
library(grid)
library(cowplot)
library(scales)
library(dplyr)
library(ggplot2)
library(reshape2)
nhs <- read_sas("nhs.sas7bdat")

vars <- c("cokh_serv", "cokr_serv", "brwn_serv", "donut_serv", 
          "cakh_serv", "cakr_serv", "srolh_serv", "srolr_serv", 
          "pieh_serv", "pier_serv")

nhs <- nhs %>%
  group_by(id) %>%
  filter(!any(if_any(all_of(vars), is.na))) %>%
  ungroup()

fd <- c("procm_serv","rmeat_serv","liver_serv","seaf_serv",
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
        "otnuts_serv","choco_serv","swt_serv","brdo_serv","swtr_serv","pie_serv",
        "cokh_serv", "cokr_serv","cakh_serv", "cakr_serv",
        "cond_serv","saldre_serv","crsoup_serv","pizza_serv")

nhs <- nhs %>%
  group_by(interval) %>%
  mutate(across(
    .cols = all_of(fd),
    .fns = ~ as.numeric(scale(log(.x + 1))),
    .names = "{.col}_log_z"
  )) %>%
  ungroup()

food_vars_log_z <- paste0(fd, "_log_z")

nhs <- nhs %>%
  group_by(id) %>%                          
  mutate(bmibase = bmicon[interval == 1][1]) %>%                          
  mutate(actbase = actcon[interval == 1][1]) %>%                          
  mutate(alcobase = alcocon[interval == 1][1]) %>%
  ungroup()

table(nhs$interval)

summary(nhs[food_vars_log_z])

nhs <- nhs %>%
  mutate(
    bmi_cat = cut(bmicon,
                  breaks = c(-Inf, 21.0,  25.0,  30.0, 32.0,  Inf),
                  right = FALSE,
                  labels = c("u21", "21_25.9", "25_30",
                             "30_32", "o32")),
    bmi_u21   = as.numeric(bmi_cat == "u21"),
    bmi_25_30 = as.numeric(bmi_cat == "25_30"),
    bmi_30_32 = as.numeric(bmi_cat == "30_32"),
    bmi_o32   = as.numeric(bmi_cat == "o32")
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

covars_to_adjust <- c("act_q2", "act_q3", "act_q4", "act_q5",
                      "bmicon","calorscon", "smkstatus2", "smkstatus3",
                      "highbp","hightc","ses","marry","dbfh","alcocon","mv")

load("final_alphas_01_to_09_withoutalc_new.RData")

performance_summary <- data.frame()
selected_vars_list <- list()

for (a_char in names(cv_fit_list)) {
  fit <- cv_fit_list[[a_char]]
  idx_1se <- which(fit$lambda == fit$lambda.1se)
  idx_min <- which(fit$lambda == fit$lambda.min)
  c_index_1se <- fit$cvm[idx_1se]
  c_index_min <- fit$cvm[idx_min]
  coef_1se <- coef(fit, s = "lambda.1se")
  active_vars <- rownames(coef_1se)[which(coef_1se != 0)]
  food_selected <- setdiff(active_vars, covars_to_adjust)
  
  if (length(food_selected) > 0) {
    selected_vars_list[[a_char]] <- data.frame(
      Alpha = a_char,
      Variable = food_selected,
      Coefficient = coef_1se[food_selected, 1]
    )
  } else {
    selected_vars_list[[a_char]] <- data.frame(
      Alpha = a_char,
      Variable = "None",
      Coefficient = NA
    )
  }
  
  temp_df <- data.frame(
    Alpha = as.numeric(a_char),
    Lambda_1se = fit$lambda.1se,
    C_index_1se = c_index_1se,
    Num_Food_Selected = length(food_selected),
    Lambda_min = fit$lambda.min,
    C_index_min = c_index_min
  )
  
  performance_summary <- rbind(performance_summary, temp_df)
}

performance_summary <- performance_summary[order(performance_summary$Alpha), ]
rownames(performance_summary) <- NULL

cat("\n=== PERFORMANCE ===\n")
print(performance_summary)

best_alpha_row <- performance_summary[which.max(performance_summary$C_index_1se), ]
cat("\n=== Best model (Based on lambda.1se's highest C-index) ===\n")
cat("Best Alpha:", best_alpha_row$Alpha, "\n")
cat("Highest C-index (1se):", round(best_alpha_row$C_index_1se, 4), "\n")
cat("Num of food groups selected:", best_alpha_row$Num_Food_Selected, "\n")

best_alpha_char <- as.character(best_alpha_row$Alpha)
print(selected_vars_list[[best_alpha_char]])

coef_df <- as.matrix(coef(cv_fit_list$`0.4`, s = "lambda.1se")) %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  rename(Coefficient = 2) %>%   
  filter(Coefficient != 0) %>% 
  mutate(Coefficient = -Coefficient) %>% 
  arrange(Coefficient) %>%
  filter(!(Variable %in% covars_to_adjust)) %>%
  mutate(Variable = factor(Variable, levels = Variable))
coef_df$Variable <- reorder(coef_df$Variable, coef_df$Coefficient)
write.csv(coef_df,"/udd/n2gji/micro/dietpt/final/data/coef_df.csv")
food_group_labels <- c(
  "cokh_serv_log_z"    = "Home-made cookies",
  "cokr_serv_log_z"    = "Ready-made cookies",
  "wine_serv_log_z"    = "Wine",
  "crsoup_serv_log_z"  = "Cream Soups",
  "yogurt_serv_log_z"  = "Yogurt",
  "fries_serv_log_z"   = "French Fries",
  "liver_serv_log_z"   = "Liver",
  "sugbev_serv_log_z"  = "Sugar-Sweetened Beverages",
  "lowbev_serv_log_z"  = "Low-Calorie Beverages",
  "tea_serv_log_z"     = "Tea",
  "rais_serv_log_z"    = "Raisins and Prunes",
  "marg_serv_log_z"    = "Margarine",
  "othfr_serv_log_z"   = "Other Fruits",
  "coff_serv_log_z"    = "Coffee",
  "procm_serv_log_z"   = "Processed Red Meats",
  "fruj_serv_log_z"    = "Fruit Juice",
  "beer_serv_log_z"    = "beer",
  "liq_serv_log_z"     = "Liquor",
  "eggs_serv_log_z"    = "Eggs",
  "toppro_serv_log_z"  = "Tomato products",
  "othveg_serv_log_z"  = "Other Vegetables",
  "appear_serv_log_z"  = "Apples and Pears",
  "dess_serv_log_z"    = "Desserts",
  "rmeat_serv_log_z"   = "Unprocessed Red Meat",
  "h2o_serv_log_z"     = "Water",
  "cruveg_serv_log_z"  = "Cruciferous Vegetables",
  "wgrain_serv_log_z"  = "Whole Grains",
  "saldre_serv_log_z"  = "Salad Dressings",
  "nightveg_serv_log_z"= "Nightshade and Cucurbitaceae Vegetables",
  "lfveg_serv_log_z"   = "Leafy Vegetables",
  "stachveg_serv_log_z"= "Starchy Vegetables",
  "rgrain_serv_log_z"  = "Refined Grains"
)

coef_df$Variable <- reorder(coef_df$Variable, coef_df$Coefficient)

p1 <- ggplot(coef_df, aes(x = Coefficient, y = Variable)) +
  geom_col(
    aes(fill = Coefficient > 0),
    width = 0.7,
    alpha = 0.9
  ) +
  geom_text(
    aes(
      label = sprintf("%.4f", Coefficient), 
      hjust = ifelse(Coefficient > 0, -0.2, 1.2) 
    ),
    size = 4,
    fontface = "bold",
    color = "black"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    linewidth = 0.7,
    color = "black",
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = c("FALSE" = "#FF6B6B", "TRUE" = "#4D8BD9"),
    guide = "none"
  ) +
  scale_x_continuous(expand = expansion(mult = 0.2)) + 
  scale_y_discrete(labels = food_group_labels) +
  labs(
    x = "Coefficient",
    y = "Food Groups"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", color = "grey85"),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
    plot.margin = unit(c(1, 1.5, 1, 1), "cm")
  )

print(p1)
ggsave(filename = "fig1.pdf", 
       plot = p1, 
       width = 16, height = 8) 

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
# plot-related data preparation
nhs <- nhs %>%
  mutate(
    bmi_cat = cut(bmicon,
                  breaks = c(-Inf, 21.0,  25.0,  30.0, 32.0,  Inf),
                  right = FALSE,
                  labels = c("u21", "21_25.9", "25_30",
                             "30_32", "o32")),
    bmi_u21   = as.numeric(bmi_cat == "u21"),
    bmi_25_30 = as.numeric(bmi_cat == "25_30"),
    bmi_30_32 = as.numeric(bmi_cat == "30_32"),
    bmi_o32   = as.numeric(bmi_cat == "o32")
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

missing_cols <- setdiff(coef_df$Variable, names(nhs))
print(missing_cols)

food_matrix <- as.matrix(nhs[, coef_df$Variable])
weights <- coef_df$Coefficient
diet_score_raw <- food_matrix %*% weights
nhs$Diet_Score <- as.vector(diet_score_raw)

nhs <- nhs %>%
  group_by(interval) %>%
  mutate(
    diet_s = as.numeric(scale(Diet_Score))
  ) %>%
  ungroup()

score_q <- c("aheimsvq","hpdirsvq","amedrsvq","dashqrsvq","ediprsvq",
             "edihrsvq","phdrsvq","drsrsvq","drscrsvq")
score_s <- c("Diet_Score","aheimsvcon","hpdirsvcon","amedrsvcon","dashqrsvcon",
             "ediprsvcon","edihrsvcon","drsrsvcon","drscrsvcon")

nhs <- nhs %>%
  group_by(interval) %>%
  mutate(across(
    .cols = all_of(score_s),
    .fns = ~ (.x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE),
    .names = "{.col}_sscore"
  )) %>%
  ungroup()

nhs <- nhs %>%
  group_by(interval) %>%
  mutate(
    diet_5 = factor(ntile(Diet_Score, 5)),
    diet_s = as.numeric(scale(Diet_Score))
  ) %>%
  ungroup()

table(nhs$marry)

nhs <- nhs %>%
  group_by(id) %>%
  mutate(marry = as.integer(any(marry == 1, na.rm = TRUE))) %>%
  ungroup()

T1 <- table1::table1(~agecon+actcon+alcocon+factor(smkstatus)+factor(bmic3cat)+
                       factor(highbp)+factor(hightc)+calorscon+aheimscon+factor(asp)+
                       bmicon+factor(dbfh)+factor(mifh)+factor(mv)+factor(marry)+ses |
                       diet_5,
                     data = nhs[nhs$interval == 1,])

T1
T1 <- data.frame(T1)
write.csv(T1,"nhs_tb1_interval1.csv")

# Cox regression main scores
score_s <- c("diet_s","aheimsvcon_sscore",
             "hpdirsvcon_sscore","amedrsvcon_sscore",
             "dashqrsvcon_sscore","ediprsvcon_sscore",
             "edihrsvcon_sscore","drsrsvcon_sscore")

covariates_list <- list(
  "Model 1" = "strata(agecon)+strata(interval)",
  "Model 2" = "strata(agecon)+factor(smkstatus)+factor(act_cat)+factor(alc)+highbp+hightc+calorscon+marry+dbfh+mv+ses+strata(interval)",
  "Model 3" = "strata(agecon)+factor(smkstatus)+bmicon+factor(alc)+factor(act_cat)+highbp+hightc+calorscon+dbfh+mv+marry+ses+strata(interval)"
)

score_q <- c("diet_5","aheimsvq","hpdirsvq","amedrsvq","dashqrsvq",
             "ediprsvq","edihrsvq","drsrsvq")

results_list_q <- list()

for (model_name in names(covariates_list)) {
  current_covariates <- covariates_list[[model_name]]
  
  cat("\n==================================\n")
  cat("Running categorical model: ", model_name, "\n")
  cat("==================================\n")
  
  for (exposure_variable in score_q) {
    cat("  - Exposure: ", exposure_variable, "\n")
    
    formula_str <- sprintf(
      "Surv(tdb2, db2) ~ as.factor(%s) + %s",
      exposure_variable,
      current_covariates
    )
    
    cox_formula <- as.formula(formula_str)
    
    model_result <- tryCatch({
      cox_model <- coxph(
        formula = cox_formula,
        cluster = id,
        data = nhs
      )
      
      tidy_output <- tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
      
      exposure_rows <- tidy_output %>% 
        filter(grepl(paste0("as.factor\\(", exposure_variable, "\\)"), term)) %>%
        mutate(
          Model = model_name,
          Variable = exposure_variable,
          Level = str_extract(term, "\\d+$") 
        )
      
      exposure_rows
    }, error = function(e) {
      cat("    !!! Model failed:", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (!is.null(model_result)) {
      key <- paste(model_name, exposure_variable, sep = "_")
      results_list_q[[key]] <- model_result
    }
  }
}

final_results_q <- bind_rows(results_list_q)

clean_results_q <- final_results_q %>%
  dplyr::select(Variable, Model, Level, estimate, conf.low, conf.high, p.value) %>%
  rename(
    Quintile = Level,
    HR = estimate,
    CI_Lower = conf.low,
    CI_Upper = conf.high,
    P_Value = p.value
  ) %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, CI_Lower, CI_Upper),
    Quintile_Label = paste0("Q", Quintile, " vs Q1")
  )

print(head(clean_results_q, 10))

ref_rows <- clean_results_q %>%
  dplyr::select(Variable, Model) %>%
  distinct() %>%
  mutate(
    Quintile = "1",
    HR_CI = "1.00 (Ref)"
  )

full_data <- bind_rows(clean_results_q, ref_rows) %>%
  mutate(Quintile = factor(Quintile, levels = c("1", "2", "3", "4", "5")))

final_table_wide <- full_data %>%
  dplyr::select(Variable, Model, Quintile, HR_CI) %>%
  pivot_wider(
    names_from = Quintile, 
    values_from = HR_CI,
    names_prefix = "Q"
  ) %>%
  arrange(Variable, Model)

print(final_table_wide)
write.csv(final_table_wide, file = "cox_nhs_q.csv")

# continuous score models
results_list <- list()

for (model_name in names(covariates_list)) {
  current_covariates <- covariates_list[[model_name]]
  
  cat("\n==================================\n")
  cat("Running continuous model: ", model_name, "\n")
  cat("==================================\n")
  
  for (exposure_variable in score_s) {
    cat("  - Exposure: ", exposure_variable, "\n")
    
    formula_str <- sprintf(
      "Surv(tdb2, db2) ~ %s + %s",
      exposure_variable,
      current_covariates
    )
    
    cox_formula <- as.formula(formula_str)
    
    model_result <- tryCatch({
      cox_model <- coxph(
        formula = cox_formula,
        cluster = id,
        data = nhs
      )
      
      tidy_output <- tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
      
      tidy_output %>% 
        filter(term == exposure_variable) %>% 
        mutate(Model = model_name)
    }, error = function(e) {
      cat("    !!! Model failed:", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (!is.null(model_result)) {
      key <- paste(model_name, exposure_variable, sep = "_")
      results_list[[key]] <- model_result
    }
  }
}

final_results_df <- bind_rows(results_list)
final_results_df

clean_results <- final_results_df %>%
  dplyr::select(Model, term, estimate, conf.low, conf.high, p.value, std.error) %>%
  rename(
    Exposure = term,
    HR = estimate,
    CI_Lower = conf.low,
    CI_Upper = conf.high,
    P_Value = p.value
  )

print(head(clean_results))

wide_table <- clean_results %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, CI_Lower, CI_Upper)
  ) %>%
  dplyr::select(Exposure, Model, HR_CI) %>%
  pivot_wider(names_from = Model, values_from = HR_CI)

print(wide_table)
write.csv(wide_table, file = "cox_nhs_con.csv")

# person-years by quintile
result_df <- data.frame(
  variable = character(),
  level = character(),
  case = integer(),
  PY = double(),
  stringsAsFactors = FALSE
)

for (var in score_q) {
  if (!is.factor(nhs[[var]])) {
    nhs[[var]] <- as.factor(nhs[[var]])
  }
  levels_v <- levels(nhs[[var]])
  for (lvl in levels_v) {
    subset_nhs <- subset(nhs, get(var) == lvl)
    case <- sum(subset_nhs$db2)
    PY <- sum(subset_nhs$tdb2)
    result_df <- rbind(result_df, data.frame(
      variable = var,
      level = lvl,
      case = case,
      PY = PY
    ))
  }
}

result_df
write.csv(result_df, "py_nhs.csv", row.names = FALSE)

# Forest plot: main patterns
data_raw <- subset(clean_results_q,
                   clean_results_q$Model == "Model 3" &
                     clean_results_q$Quintile_Label == "Q5 vs Q1")
plot_data <- data_raw

name_map <- c(
  "diet_5"   = "DDP",
  "aheimsvq" = "AHEI",
  "hpdirsvq" = "HPDI",
  "amedrsvq" = "aMED",
  "ediprsvq" = "rEDIP",
  "edihrsvq" = "rEDIH",
  "drsrsvq"  = "DRRD",
  "dashqrsvq" = "DASH"
)

plot_data$Variable <- as.character(plot_data$Variable)
plot_data$Variable <- ifelse(plot_data$Variable %in% names(name_map), 
                             name_map[plot_data$Variable], 
                             plot_data$Variable)

plot_data$Variable <- factor(plot_data$Variable, levels = plot_data$Variable)
ordered_vars <- plot_data$Variable

color_palette <- scales::hue_pal()(10)

forest_df <- plot_data %>%
  dplyr::select(Variable, HR_CI) %>%
  dplyr::rename(
    `Dietary Patterns` = Variable,
    `HR (95% CI)` = HR_CI
  ) %>%
  mutate(` ` = paste(rep(" ", 20), collapse = " "))

tm <- forest_theme(
  base_size = 10,
  ci_pch = 19,
  ci_lwd = 2,
  refline_col = "grey50",
  core = list(bg_params = list(fill = c("white")))
)

p_forest <- forest(
  data = forest_df,
  est = plot_data$HR,
  lower = plot_data$CI_Lower,
  upper = plot_data$CI_Upper,
  ci_column = 3,
  ref_line = 1,
  xlim = c(0, 1.1),
  ticks_at = c(0.2, 0.4, 0.6, 0.8, 1), 
  theme = tm
)
p_forest

for (i in 1:nrow(plot_data)) {
  p_forest <- edit_plot(
    p_forest,
    row = i,
    col = 3,
    which = "ci",
    gp = gpar(col = color_palette[i], fill = color_palette[i])
  )
}
p_forest

# heatmap of HR differences
plot_data <- plot_data %>%
  mutate(
    logHR = log(HR),
    SE = (log(CI_Upper) - log(CI_Lower)) / 3.92
  )

calculate_p_diff <- function(i, j, data) {
  beta1 <- data$logHR[i]; se1 <- data$SE[i]
  beta2 <- data$logHR[j]; se2 <- data$SE[j]
  se_diff <- sqrt(se1^2 + se2^2)
  z_score <- (beta1 - beta2) / se_diff
  p_val <- 2 * (1 - pnorm(abs(z_score)))
  return(p_val)
}

n <- nrow(plot_data)
p_matrix <- matrix(NA, nrow = n, ncol = n)
diff_matrix <- matrix(NA, nrow = n, ncol = n)

for (i in 1:n) {
  for (j in 1:n) {
    diff_matrix[i, j] <- plot_data$HR[i] - plot_data$HR[j]
    p_matrix[i, j] <- calculate_p_diff(i, j, plot_data)
  }
}

rownames(p_matrix) <- plot_data$Variable
colnames(p_matrix) <- plot_data$Variable
rownames(diff_matrix) <- plot_data$Variable
colnames(diff_matrix) <- plot_data$Variable

heatmap_data <- reshape2::melt(diff_matrix)
colnames(heatmap_data) <- c("Var1", "Var2", "Diff")
p_data_long <- reshape2::melt(p_matrix)
heatmap_data$P_Value_Diff <- p_data_long$value

heatmap_data <- heatmap_data %>%
  mutate(
    stars = case_when(
      P_Value_Diff < 0.001 ~ "***",
      P_Value_Diff < 0.01  ~ "**",
      P_Value_Diff < 0.05  ~ "*",
      TRUE ~ ""
    ),
    label_text = paste0(sprintf("%.2f", Diff), stars),
    label_text = ifelse(Var1 == Var2, "0", label_text)
  )

heatmap_data$Var1 <- factor(heatmap_data$Var1, levels = rev(plot_data$Variable))
heatmap_data$Var2 <- factor(heatmap_data$Var2, levels = plot_data$Variable)

heatmap_data$Var1 <- factor(heatmap_data$Var1, levels = rev(levels(plot_data$Variable)))

p_heatmap <- ggplot(heatmap_data, aes(x = Var2, y = Var1, fill = Diff)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label_text), size = 3, color = "black", lineheight = 0.8) +
  scale_fill_gradient2(
    low = "#F8766D", mid = "white", high = "#7CAE00", midpoint = 0,
    name = "Mean difference in HRs"
  ) +
  scale_x_discrete(position = "top") +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    axis.text.x.top = element_text(
      angle = 45, hjust = 0, vjust = 0,
      margin = margin(b = 6)
    ),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  )

p_heatmap

final_plot <- plot_grid(
  grid::grid.grabExpr(grid::grid.draw(p_forest)), 
  p_heatmap,
  ncol = 2,
  rel_widths = c(1.6, 1),
  labels = "AUTO"
)

print(final_plot)
ggsave("fig2_nhs_forest.pdf", p_forest, width = 5, height = 4)
ggsave("fig2_nhs_heatmap.pdf", p_heatmap, width = 12, height = 8)
