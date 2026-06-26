########In LVS##########
rm(list = ls())
gc()
options(bitmapType='cairo')
setwd("/udd/n2gji/micro/dietpt/final/data/")
library(glmnet)
library(survival)
library(dplyr)
library(haven)
library(patchwork)
library(forestploter)
library(dplyr)
library(tibble)
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
# install.packages("forestploter")
library(forestploter)
library(grid)
library(cowplot)
library(scales)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(RColorBrewer)
setwd("/udd/n2gji/micro/dietpt/final/data")
covars_to_adjust <- c( "act_q2", "act_q3", "act_q4", "act_q5","bmicon","calorscon", "smkstatus2", "smkstatus3","highbp","hightc","ses","marry","dbfh","alcocon","mv") 
coef_df<-read.csv("/udd/n2gji/micro/dietpt/final/data/coef_df.csv")
###############Rep in LVS#################
load("/udd/n2gji/data/lvs.RData")
foodgroup<-read.csv('/udd/n2gji/micro/dietpt/lvs_diet_new.csv')
foodgroup<-foodgroup[,c("id","procm07","rmeat07","liver07","seaf07","sseaf07","poul07","eggs07","butter07","marg07","milk07","yogurt07","cheese07","daiswt07","wine07","beer07_serv","liq07_serv","tea07","coff07","h2o07","sugbev07","lowbev07","citf07","fruj07","berr07","appear07","ban07",
                        "rais07","peach07","othfr07","lfveg07","cruveg07","stachveg07",
                        "nightveg07","legu07","othveg07","tompro07","apiveg07",
                        "fries07","wgrain07","rgrain07","snack07","pnuts07","otnuts07",
                        "choco07","swt07","dess07","cond07","saldre07","crsoup07","pizza07","cokh07","cokr07")]
temp_lvs$white<-ifelse(temp_lvs$race=="White",1,0)
temp_lvs<-merge(temp_lvs,foodgroup,by="id",all.x = T)
food_group <- c("procm07","rmeat07","liver07","seaf07","sseaf07","poul07","eggs07","butter07",
                "marg07","milk07","yogurt07","cheese07","daiswt07","wine07","beer07_serv","liq07_serv","tea07",
                "coff07","h2o07","sugbev07","lowbev07","citf07","fruj07","berr07","appear07","ban07",
                "rais07","peach07","othfr07","lfveg07","cruveg07","stachveg07",
                "nightveg07","legu07","othveg07","tompro07","apiveg07",
                "fries07","wgrain07","rgrain07","snack07","pnuts07","otnuts07",
                "choco07","swt07","dess07","cond07","saldre07","crsoup07","pizza07","cokh07","cokr07")
summary(temp_lvs[food_group])
# procm07           rmeat07          liver07             seaf07          sseaf07            poul07           eggs07        
# Min.   :0.00000   Min.   :0.0000   Min.   :0.000000   Min.   :0.0000   Min.   :0.00000   Min.   :0.0000   Min.   :-0.05571  
# 1st Qu.:0.04209   1st Qu.:0.3005   1st Qu.:0.000000   1st Qu.:0.0000   1st Qu.:0.00000   1st Qu.:0.1890   1st Qu.: 0.10821  
# Median :0.19312   Median :0.5667   Median :0.000000   Median :0.1134   Median :0.00000   Median :0.3808   Median : 0.27803  
# Mean   :0.29475   Mean   :0.6758   Mean   :0.006396   Mean   :0.1665   Mean   :0.06600   Mean   :0.4272   Mean   : 0.35970  
# 3rd Qu.:0.40462   3rd Qu.:0.9265   3rd Qu.:0.000000   3rd Qu.:0.2460   3rd Qu.:0.08239   3rd Qu.:0.5827   3rd Qu.: 0.49992  
# Max.   :2.58111   Max.   :3.6680   Max.   :1.017857   Max.   :1.2109   Max.   :1.30417   Max.   :2.2928   Max.   : 3.60571  
# butter07           marg07           milk07          yogurt07         cheese07         daiswt07           wine07      
# Min.   : 0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.00000   Min.   :0.0000  
# 1st Qu.: 0.2063   1st Qu.:0.0000   1st Qu.:0.0820   1st Qu.:0.0000   1st Qu.:0.3695   1st Qu.:0.00000   1st Qu.:0.0000  
# Median : 0.8803   Median :0.0000   Median :0.3660   Median :0.1036   Median :0.6275   Median :0.08504   Median :0.1702  
# Mean   : 2.0428   Mean   :0.3378   Mean   :0.5666   Mean   :0.1893   Mean   :0.7608   Mean   :0.14568   Mean   :0.4902  
# 3rd Qu.: 2.2823   3rd Qu.:0.4301   3rd Qu.:0.8083   3rd Qu.:0.2859   3rd Qu.:1.0053   3rd Qu.:0.20141   3rd Qu.:0.7550  
# Max.   :37.2429   Max.   :9.1091   Max.   :4.7586   Max.   :2.0815   Max.   :4.5495   Max.   :1.30103   Max.   :5.4545  
# alcob07             tea07            coff07           h2o07           sugbev07         lowbev07           citf07       
# Min.   : 0.00000   Min.   : 0.000   Min.   :0.0000   Min.   : 0.000   Min.   :0.0000   Min.   : 0.0000   Min.   :0.00000  
# 1st Qu.: 0.00000   1st Qu.: 0.000   1st Qu.:0.3055   1st Qu.: 1.277   1st Qu.:0.0000   1st Qu.: 0.0000   1st Qu.:0.00000  
# Median : 0.06943   Median : 0.348   Median :1.5259   Median : 2.509   Median :0.0430   Median : 0.1118   Median :0.00000  
# Mean   : 0.43077   Mean   : 0.806   Mean   :1.5540   Mean   : 2.913   Mean   :0.1708   Mean   : 0.7623   Mean   :0.09401  
# 3rd Qu.: 0.48518   3rd Qu.: 1.180   3rd Qu.:2.3320   3rd Qu.: 4.014   3rd Qu.:0.1816   3rd Qu.: 1.0051   3rd Qu.:0.11450  
# Max.   :15.17546   Max.   :11.019   Max.   :8.1663   Max.   :16.215   Max.   :4.5467   Max.   :17.0000   Max.   :1.85769  
# fruj07            berr07          appear07          ban07            rais07           peach07           othfr07       
# Min.   :0.00000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.00000   Min.   :0.00000   Min.   :0.00000  
# 1st Qu.:0.00127   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.00000   1st Qu.:0.00000   1st Qu.:0.00000  
# Median :0.17225   Median :0.1250   Median :0.0821   Median :0.1589   Median :0.05265   Median :0.00000   Median :0.02163  
# Mean   :0.39186   Mean   :0.2670   Mean   :0.1551   Mean   :0.2453   Mean   :0.21240   Mean   :0.08526   Mean   :0.09724  
# 3rd Qu.:0.63430   3rd Qu.:0.3727   3rd Qu.:0.2238   3rd Qu.:0.3670   3rd Qu.:0.27948   3rd Qu.:0.09595   3rd Qu.:0.12248  
# Max.   :3.61058   Max.   :3.2946   Max.   :1.7598   Max.   :2.1312   Max.   :2.65879   Max.   :2.77806   Max.   :2.07992  
# lfveg07          cruveg07         stachveg07       nightveg07         legu07           othveg07          tompro07       
# Min.   :0.0000   Min.   :0.00000   Min.   :0.0000   Min.   :0.0000   Min.   :0.00000   Min.   :0.00000   Min.   :  0.0000  
# 1st Qu.:0.0803   1st Qu.:0.05848   1st Qu.:0.1066   1st Qu.:0.3035   1st Qu.:0.07878   1st Qu.:0.07748   1st Qu.:  0.0298  
# Median :0.2006   Median :0.17548   Median :0.2137   Median :0.6987   Median :0.19181   Median :0.25058   Median :  0.2900  
# Mean   :0.2807   Mean   :0.26525   Mean   :0.2660   Mean   :0.9466   Mean   :0.27750   Mean   :0.37450   Mean   : 13.9151  
# 3rd Qu.:0.3922   3rd Qu.:0.37289   3rd Qu.:0.3659   3rd Qu.:1.3357   3rd Qu.:0.38108   3rd Qu.:0.51616   3rd Qu.: 11.3023  
# Max.   :2.4361   Max.   :4.76571   Max.   :2.0282   Max.   :8.6862   Max.   :2.19670   Max.   :3.98297   Max.   :403.1529  
# apiveg07          fries07            wgrain07         rgrain07         snack07          pnuts07          otnuts07       
# Min.   :0.00000   Min.   :0.000000   Min.   :  0.00   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.000000  
# 1st Qu.:0.02778   1st Qu.:0.000000   1st Qu.:  1.78   1st Qu.:0.8466   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.006526  
# Median :0.16071   Median :0.005173   Median : 21.94   Median :1.3358   Median :0.1276   Median :0.1690   Median :0.130070  
# Mean   :0.35035   Mean   :0.039224   Mean   : 33.06   Mean   :1.5225   Mean   :0.2356   Mean   :0.3433   Mean   :0.296248  
# 3rd Qu.:0.44659   3rd Qu.:0.055114   3rd Qu.: 54.00   3rd Qu.:2.0296   3rd Qu.:0.3290   3rd Qu.:0.4554   3rd Qu.:0.378861  
# Max.   :6.95839   Max.   :0.742760   Max.   :259.02   Max.   :6.7943   Max.   :2.3742   Max.   :4.4755   Max.   :4.765306  
# choco07            swt07             dess07           cond07            saldre07         crsoup07          pizza07       
# Min.   :0.00000   Min.   :0.00000   Min.   :0.0000   Min.   : 0.00000   Min.   :0.0000   Min.   :0.00000   Min.   :0.00000  
# 1st Qu.:0.00000   1st Qu.:0.08862   1st Qu.:0.2190   1st Qu.: 0.00000   1st Qu.:0.2502   1st Qu.:0.00000   1st Qu.:0.00000  
# Median :0.01948   Median :0.27918   Median :0.6178   Median : 0.00000   Median :0.4450   Median :0.00000   Median :0.00000  
# Mean   :0.12083   Mean   :0.41050   Mean   :0.8813   Mean   : 0.87402   Mean   :0.5647   Mean   :0.04048   Mean   :0.07795  
# 3rd Qu.:0.11672   3rd Qu.:0.56720   3rd Qu.:1.2203   3rd Qu.: 0.04762   3rd Qu.:0.7540   3rd Qu.:0.06056   3rd Qu.:0.12881  
# Max.   :6.72890   Max.   :4.84643   Max.   :9.3555   Max.   :24.84524   Max.   :5.0914   Max.   :0.64705   Max.   :0.86678
zero_vals <- colMeans(temp_lvs[food_group] == 0, na.rm = TRUE)
zero_table <- data.frame(
  Food_Item = names(zero_vals),
  Proportion = as.numeric(zero_vals),
  Percentage = paste0(round(zero_vals * 100, 2), "%")
)
zero_table <- zero_table[order(zero_table$Proportion, decreasing = TRUE), ]
rownames(zero_table) <- NULL
head(zero_table)
# Food_Item Proportion Percentage
# 1   liver07  0.9640468      96.4%
# 2    cond07  0.7391304     73.91%
# 3  crsoup07  0.6655518     66.56%
# 4   peach07  0.6170569     61.71%
# 5   sseaf07  0.5535117     55.35%
# 6    marg07  0.5484950     54.85%
food_group <- c("procm07","rmeat07","liver07","seaf07","sseaf07","poul07","eggs07","butter07",
                "marg07","milk07","yogurt07","cheese07","daiswt07","wine07","beer07_serv","liq07_serv","tea07",
                "coff07","h2o07","sugbev07","lowbev07","citf07","fruj07","berr07","appear07","ban07",
                "rais07","peach07","othfr07","lfveg07","cruveg07","stachveg07",
                "nightveg07","legu07","othveg07","tompro07","apiveg07",
                "fries07","wgrain07","rgrain07","snack07","pnuts07","otnuts07",
                "choco07","swt07","dess07","cond07","saldre07","crsoup07","pizza07","cokh07","cokr07")

fd <- c("procm_serv","rmeat_serv","liver_serv","seaf_serv","sseaf_serv","poul_serv",
        "eggs_serv","butter_serv","marg_serv","milk_serv","yogurt_serv","cheese_serv",
        "daiswt_serv","wine_serv","beer_serv","liq_serv","tea_serv","coff_serv","h2o_serv","sugbev_serv",
        "lowbev_serv","citf_serv","fruj_serv","berr_serv","appear_serv","ban_serv","rais_serv",
        "peach_serv","othfr_serv","lfveg_serv","cruveg_serv","stachveg_serv","nightveg_serv",
        "legu_serv","othveg_serv","tompro_serv","apiveg_serv","fries_serv","wgrain_serv",
        "rgrain_serv","snack_serv","pnuts_serv","otnuts_serv","choco_serv","swt_serv",
        "dess_serv","cond_serv","saldre_serv","crsoup_serv","pizza_serv","cokh_serv","cokr_serv")
name_mapping <- setNames(fd,food_group)
temp_lvs$dietary_score <- 0
for(i in 1:nrow(coef_df)) {
  
  target_var <- coef_df$Variable[i]
  weight <- coef_df$Coefficient[i]
  
  base_name <- sub("_log_z", "", target_var)
  
  original_col_name <- names(name_mapping)[which(name_mapping == base_name)]
  
  if(length(original_col_name) > 0 && original_col_name %in% names(temp_lvs)) {
    
    raw_data <- temp_lvs[[original_col_name]]
    
    log_data <- log(raw_data + 1)
    
    z_data <- as.numeric(scale(log_data))
    
    term_value <- z_data * weight
    
    temp_lvs$dietary_score <- temp_lvs$dietary_score + term_value
    
    print(paste("已计算:", base_name, "-> 原始列:", original_col_name, "| 权重:", weight))
    
  } else {
    print(paste("⚠️ 警告: 未在 temp_lvs 或 mapping 中找到对应列:", base_name))
  }
}
summary(temp_lvs$dietary_score)
#######################Some covaritaes######################
nhs <- read_sas("nhs.sas7bdat")
table(temp_lvs$blddate,temp_lvs$cohort)
summary(nhs$interval)
nhs_2010<-subset(nhs,nhs$interval==13)
table(nhs_2010$id %in% temp_lvs$id)
table(temp_lvs$cohort) #281 in 333
rm(nhs)
gc()
#HPFS
hpfs<-read_sas("hpfs.sas7bdat")
table(temp_lvs$blddate,temp_lvs$cohort)
table(hpfs$rtmnyr12)
table(hpfs$interval)
hpfs_2010<-subset(hpfs,hpfs$interval==14)
table(hpfs_2010$id %in% temp_lvs$id)#451 in 453
table(temp_lvs$cohort)
rm(hpfs)
gc()
#NHSII
nhs2<-read_sas("nhs2.sas7bdat")
table(temp_lvs$blddate,temp_lvs$cohort)
table(nhs2$retmo11)
table(nhs2$interval)
nhs2_2010<-subset(nhs2,nhs2$interval==11)
table(nhs2_2010$id %in% temp_lvs$id)#359 in 410
table(temp_lvs$cohort)
rm(nhs2)
gc()
hpfs_2010$ses<-hpfs_2010$nSES
var<-c("highbp","hightc","ses","marry")
nhs2_2010<-nhs2_2010 %>% dplyr::select(id,all_of(var))
nhs_2010<-nhs_2010 %>% dplyr::select(id,all_of(var))
hpfs_2010<-hpfs_2010 %>% dplyr::select(id,all_of(var))

nhs2_2010$id<-nhs2_2010$id*1000
nhs2_2010$id<-as.numeric(nhs2_2010$id)
nhs_2010$id<-nhs_2010$id*100
nhs_2010$id<-as.numeric(nhs_2010$id)
hpfs_2010$id<-hpfs_2010$id*10000
hpfs_2010$id<-as.numeric(hpfs_2010$id)
lvs_ffq_food<- rbind(nhs_2010,nhs2_2010)# 1466
lvs_ffq_food<- rbind(lvs_ffq_food,hpfs_2010)# 1466
lvs_ffq_food$id<-as.numeric(lvs_ffq_food$id)
###################
class(temp_lvs$id)
class(lvs_ffq_food$id)
temp_lvs$id_new <- as.numeric(as.character(temp_lvs$id))
temp_lvs$id_new <- ifelse(temp_lvs$cohort == "nhs1", temp_lvs$id_new * 100,
                          ifelse(temp_lvs$cohort == "nhs2", temp_lvs$id_new * 1000,
                                 temp_lvs$id_new * 10000))


temp_lvs<-merge(temp_lvs,lvs_ffq_food,by.x="id_new",by.y = "id",all.x = T)
############Biomarker###########
lvs_bio<-read.csv("/udd/n2gji/micro/data/basic_bio.csv")
temp_lvs<-merge(temp_lvs,lvs_bio[,c("id", "hba1c","tg","tc","hdlc","crp","scd14","lbp" ,"non_hdlc","tg.hdlc","non_hdlcq","hdlcq","tgq","crpq","hba1cq","mscore_p")],by="id",all.x = T)
x_val <- temp_lvs[, c("dietary_score"), drop = FALSE]
class(x_val)

y_val <- temp_lvs[, c("hba1c", "tg", "tc", "hdlc", "crp", "tg.hdlc")]

cor_matrix <- cor(y_val, x_val, use = "pairwise.complete.obs")

p_matrix <- matrix(NA, nrow = ncol(y_val), ncol = ncol(x_val),
                   dimnames = list(colnames(y_val), colnames(x_val)))

for (row_name in rownames(p_matrix)) {
  for (col_name in colnames(p_matrix)) {
    complete_obs_count <- sum(complete.cases(y_val[[row_name]], x_val[[col_name]]))
    if (complete_obs_count > 2) {
      test_result <- cor.test(y_val[[row_name]], x_val[[col_name]], use = "pairwise.complete.obs")
      p_matrix[row_name, col_name] <- test_result$p.value
    } 
  }
}

plot_data <- data.frame(
  Variable = rownames(cor_matrix),
  Correlation = cor_matrix[, 1],
  P_Value = p_matrix[, 1]
)

plot_data <- plot_data %>%
  mutate(
    Significance = case_when(
      is.na(P_Value) ~ "",
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    Direction = ifelse(Correlation > 0, "Positive", "Negative"),
    
    Variable_Full = recode(
      Variable,
      hba1c   = "HbA1c",
      tg      = "Triglycerides",
      tc      = "Total Cholesterol",
      hdlc    = "HDL-C",
      crp     = "CRP",
      tg.hdlc = "TG/HDL-C Ratio",
      .default = Variable
    )
  )

print(plot_data)

p<-ggplot(plot_data, aes(x = reorder(Variable_Full, Correlation), y = Correlation, fill = Direction)) +
  
  geom_col(width = 0.7) + 
  
  geom_text(
    aes(label = Significance, y = Correlation + 0.02 * sign(Correlation)), 
    vjust = ifelse(plot_data$Correlation > 0, 0, 1),
    size = 5,
    color = "black"
  ) +
  
  coord_flip() + 
  
  scale_fill_manual(values = c("Positive" = "#d6604d", "Negative" = "#4393c3")) +
  
  labs(
    title = "'",
    subtitle = "* p<0.05, ** p<0.01, *** p<0.001",
    x = "Biomarker / Diet Score",
    y = "Pearson Correlation Coefficient"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  
  ylim(min(plot_data$Correlation) - 0.05, max(plot_data$Correlation) + 0.05)
ggsave("ddp_biomarker.png", plot = p, width = 10, height = 8, dpi = 300, bg = "white")
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
met_list = as.character(colnames(temp_lvs[, 3:297]))
print(met_list)
results_lvs=data.frame(met=met_list)
met_aligned=results_lvs %>% pull(met)
temp_lvs$score_s<-scale(temp_lvs$dietary_score)
temp_lvs$highbp<-ifelse(is.na(temp_lvs$highbp),0,temp_lvs$highbp)
temp_lvs$hightc<-ifelse(is.na(temp_lvs$hightc),0,temp_lvs$hightc)
temp_lvs$marry<-ifelse(is.na(temp_lvs$marry),1,temp_lvs$marry)
temp_lvs$ses<-ifelse(is.na(temp_lvs$ses),median(temp_lvs$ses,na.rm = T),temp_lvs$ses)
for (i in seq_along(met_aligned)) {
  
  formula_str <- paste0("`", met_aligned[i], "` ~ score_s + cohort + ageyr + factor(bmi_group) + factor(neversmoking) + factor(calor_g) + factor(act_g)+ highbp + hightc + marry + ses + factor(alco_g)")
  
  fit = lm(as.formula(formula_str), data = temp_lvs)
  summ = summary(fit)
  
  results_lvs[i, "beta_sc_adj"] = round(summ$coefficients[2, 1], 4)
  results_lvs[i, "se_sc_adj"]   = summ$coefficients[2, 2]
  results_lvs[i, "p_sc_adj"]    = summ$coefficients[2, 4]
  results_lvs[i, "n_sc_adj"]    = nrow(model.frame(fit))
}

head(results_lvs)
results_lvs<- results_lvs %>% mutate(
  sc_intake=case_when(beta_sc_adj>0 & p_sc_adj<0.05/295  ~ 'positive',
                      beta_sc_adj<0 & p_sc_adj<0.05/295  ~ 'negative',
                      T ~ 'insignificant')
)
load("/udd/n2gji/f_lvs_match.RData")
results_lvs<-merge(results_lvs, f_lvs_match, by = 'met', all = FALSE, sort = TRUE)
results_lvs$p_sdr<-p.adjust(results_lvs$p_sc_adj,method = "BH")
write.csv(results_lvs, "/udd/n2gji/micro/dietpt/db2/output/score_met_lvs.csv")
results_vol <- results_lvs
alpha_adjusted <- 0.05 / 295

# # sc
# install.packages("ggrepel")
library(ggrepel)
results_vol <- transform(results_vol, color = ifelse(sc_intake == "positive", "#EB4232",
                                                     ifelse(sc_intake == "negative", "#1597A5", "#d8d8d8")),
                         significance = ifelse(sc_intake == "positive" | sc_intake == "negative", "*", ""))

significant_metabolites <- subset(results_vol, p_sc_adj < alpha_adjusted )

w <- ggplot(results_vol, 
            aes(x = beta_sc_adj, y = -log10(p_sc_adj), 
                color = color, label = significance)) +
  geom_point() +
  geom_hline(yintercept = -log10(alpha_adjusted), linetype = "dashed", color = "black") + 
  geom_vline(xintercept = c(-0.01, 0.01), linetype = "dashed", color = "black") + 
  theme_minimal() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold")) +
  labs(
    x = "Adjusted-beta",
    y = "-log10(p-value)",
    color = "Metabolites") +
  scale_color_manual(values = c("#EB4232" = "#EB4232", 
                                "#1597A5" = "#1597A5", 
                                "#d8d8d8" = "#d8d8d8"),
                     labels = c("#EB4232"="Positive", 
                                "#1597A5"="Negative", 
                                "#d8d8d8"="Insignificant"),
                     breaks = c("#EB4232", "#1597A5", "#d8d8d8"))

Up <- filter(results_vol, sc_intake == 'positive') %>% distinct(met, .keep_all = T) %>% top_n(10, -log10(p_sc_adj))

Down <- filter(results_vol, sc_intake == 'negative') %>% distinct(met, .keep_all = T) %>% top_n(10, -log10(p_sc_adj))

nudge_x_up = 0.01 - Up$beta_sc_adj 
nudge_x_down = -0.01 - Down$beta_sc_adj

p3 <- w + 
  geom_point(data = Up,aes(x = beta_sc_adj, y = -log10(p_sc_adj)),
             color = '#1597A5', size = 7.5, alpha = 0.2) +
  geom_text_repel(data = Up,aes(x = beta_sc_adj, y = -log10(p_sc_adj), label = metabolite_name),
                  seed = 2024,color = 'black',show.legend = FALSE, 
                  min.segment.length = 0,
                  segment.linetype = 3, 
                  nudge_x = nudge_x_up, 
                  direction = "y", 
                  hjust = 0, 
                  force = 3,
                  force_pull = 2,
                  size = 4,
                  max.overlaps = Inf,
                  fontface = "bold")
p4 <- p3 + 
  geom_point(data = Down,aes(x = beta_sc_adj, y = -log10(p_sc_adj)),
             color = '#b4b4d8', size = 7.5, alpha = 0.2) +
  geom_text_repel(data = Down,aes(x = beta_sc_adj, y = -log10(p_sc_adj), label = metabolite_name),
                  seed = 2024,color = 'black',show.legend = FALSE, 
                  min.segment.length = 0,
                  segment.linetype = 3, 
                  nudge_x = nudge_x_down, 
                  direction = "y", 
                  hjust = 0.5, 
                  force = 3,
                  force_pull = 2,
                  size = 4,
                  max.overlaps = Inf,
                  fontface = "bold")
p4
ggsave(p4,file="met_volcano.pdf")
##########Plot#############
#Part 1: Metabolite and score
#Part 2: Individual food and metabolite
#Part 3: Metabolite and cardiometabolite biomarker
results_lvs_sub<-subset.data.frame(results_lvs,results_lvs$p_sc_adj<0.05/295)
names(results_lvs_sub)
mapping_df <- enframe(name_mapping, name = "Original_Name", value = "Current_Name")
coef_df$new_name<-coef_df$Variable
coef_df$new_name<-sub("_log_z", "", coef_df$new_name)
coef_df <- coef_df %>%
  left_join(mapping_df, by = c("new_name" = "Current_Name"))
met_list=as.character(results_lvs_sub$met)
food_list=as.character(coef_df$Original_Name)

total_rows <- length(food_list) * length(met_list)
results_food <- data.frame(
  Food = character(total_rows),
  Metabolite = character(total_rows),
  Beta = numeric(total_rows),
  SE = numeric(total_rows),
  P_value = numeric(total_rows),
  N_obs = numeric(total_rows),
  stringsAsFactors = FALSE
)
covariates <- "+ cohort + ageyr + factor(bmi_group) + factor(neversmoking) + factor(calor_g) + factor(act_g)+highbp+hightc+marry+ses+factor(alco_g)"

row_idx <- 1 

for (f_item in food_list) {
  
  temp_lvs$current_food_s <- as.numeric(scale(temp_lvs[[f_item]]))
  
  for (m_item in met_list) {
    
    formula_str <- paste0(m_item, " ~ current_food_s ", covariates)
    
    fit <- lm(as.formula(formula_str), data = temp_lvs)
    summ <- summary(fit)
    
    beta <- summ$coefficients[2, 1]
    se   <- summ$coefficients[2, 2]
    pval <- summ$coefficients[2, 4]
    n_obs <- nrow(model.frame(fit))
    
    results_food[row_idx, "Food"]       <- f_item
    results_food[row_idx, "Metabolite"] <- m_item
    results_food[row_idx, "Beta"]       <- round(beta, 4)
    results_food[row_idx, "SE"]         <- round(se, 4)
    results_food[row_idx, "P_value"]    <- pval
    results_food[row_idx, "N_obs"]      <- n_obs
    
    row_idx <- row_idx + 1
  }
}
head(results_food)
results_food<-merge(results_food, f_lvs_match, by.x = 'Metabolite',by.y="met", all = FALSE, sort = TRUE)
results_food$P_FDR <- p.adjust(results_food$P_value, method = "fdr")

met_list <- as.character(results_lvs_sub$met)
bio_list <- c("hba1c", "tg", "tc", "hdlc", "crp", "tg.hdlc")

total_rows <- length(bio_list) * length(met_list)
results_bio <- data.frame(
  Biomarker = character(total_rows),
  Metabolite = character(total_rows),
  Beta = numeric(total_rows),
  SE = numeric(total_rows),
  P_value = numeric(total_rows),
  N_obs = numeric(total_rows),
  Error_Msg = character(total_rows),
  Used_Covars = character(total_rows),
  stringsAsFactors = FALSE
)
potential_covars <- c("cohort", "ageyr","bmi_group", "neversmoking", "calor_g", "act_g","highbp","hightc","marry","ses","alco_g")
covar_formula_map <- c(
  "cohort" = "cohort", 
  "ageyr" = "ageyr", 
  "bmi_group" = "factor(bmi_group)",
  "neversmoking" = "factor(neversmoking)",
  "calor_g" = "factor(calor_g)",
  "act_g" = "factor(act_g)",
  "marry"="marry",
  "ses"="ses",
  "hightc"="hightc",
  "highbp"="highbp",
  "alco_g"="factor(alco_g)"
)
row_idx <- 1 

for (b_item in bio_list) {
  
  raw_data <- temp_lvs[[b_item]]
  
  valid_indices <- !is.na(raw_data)
  
  if (sum(valid_indices) < 30) { 
    print(paste("Skipping", b_item, "- N too small"))
    row_idx <- row_idx + length(met_list)
    next 
  }
  
  valid_covars_list <- c()
  
  for (cov in potential_covars) {
    cov_data <- temp_lvs[[cov]][valid_indices]
    
    n_unique <- length(unique(na.omit(cov_data)))
    
    if (n_unique >= 2) {
      valid_covars_list <- c(valid_covars_list, covar_formula_map[cov])
    } else {
    }
  }
  
  if (length(valid_covars_list) > 0) {
    final_covariates_str <- paste("+", paste(valid_covars_list, collapse = " + "))
  } else {
    final_covariates_str <- ""
  }
  
  temp_lvs$bio_list_s <- as.numeric(scale(raw_data))
  
  for (m_item in met_list) {
    
    formula_str <- paste0(m_item, " ~ bio_list_s ", final_covariates_str)
    
    fit <- tryCatch({
      lm(as.formula(formula_str), data = temp_lvs, na.action = na.omit)
    }, error = function(e) {
      results_bio[row_idx, "Error_Msg"] <<- as.character(e$message)
      return(NULL)
    })
    
    if (!is.null(fit)) {
      summ <- summary(fit)
      if (nrow(summ$coefficients) >= 2) {
        results_bio[row_idx, "Beta"]    <- round(summ$coefficients[2, 1], 4)
        results_bio[row_idx, "SE"]      <- round(summ$coefficients[2, 2], 4)
        results_bio[row_idx, "P_value"] <- summ$coefficients[2, 4]
        results_bio[row_idx, "N_obs"]   <- nrow(model.frame(fit))
        results_bio[row_idx, "Error_Msg"] <- "Success"
        results_bio[row_idx, "Used_Covars"] <- final_covariates_str
      } else {
        results_bio[row_idx, "Error_Msg"] <- "Coef NA"
      }
    }
    
    results_bio[row_idx, "Biomarker"]  <- b_item
    results_bio[row_idx, "Metabolite"] <- m_item
    
    row_idx <- row_idx + 1
  }
  print(paste("完成:", b_item, "| 剔除了常量协变量, 剩余协变量数:", length(valid_covars_list)))
}
results_bio$P_FDR <- p.adjust(results_bio$P_value, method = "fdr")
results_bio<-merge(results_bio, f_lvs_match, by.x = 'Metabolite',by.y="met", all = FALSE, sort = TRUE)
unique(results_food$Food)
food_group_labels <- c(
  "beer07_serv"    = "Beers",
  "liq07_serv"    = "Liquors",
  "sseaf07"    = "Shell seafood",
  "cokh07"    = "Home-made cookie",
  "cokr07"    = "Ready-made cookie",
  
  "wine07"     = "Wine",
  "crsoup07"   = "Cream Soups",
  "yogurt07"   = "Yogurt",
  "fries07"    = "French Fries",
  "liver07"    = "Liver",
  "sugbev07"   = "Sugar-Sweetened Beverages",
  "lowbev07"   = "Low-Calorie Beverages",
  "tea07"      = "Tea",
  "rais07"     = "Raisins and Prunes",
  "marg07"     = "Margarine",
  "othfr07"    = "Other Fruits",
  "coff07"     = "Coffee",
  "procm07"    = "Processed Red Meats",
  "fruj07"     = "Fruit Juice",
  "eggs07"     = "Eggs",
  "othveg07"   = "Other Vegetables",
  "appear07"   = "Apples and Pears",
  "dess07"     = "Desserts",
  "rmeat07"    = "Unprocessed Red Meat",
  "h2o07"      = "Water",
  "cruveg07"   = "Cruciferous Vegetables",
  "wgrain07"   = "Whole Grains",
  "saldre07"   = "Salad Dressings",
  "nightveg07" = "Nightshade and Cucurbitaceae Vegetables",
  "lfveg07"    = "Leafy Vegetables",
  "stachveg07" = "Starchy Vegetables",
  "rgrain07"   = "Refined Grains"
)

biomarker_labels <- c(
  "hba1c"="HbA1c",
  "tg.hdlc"="TG/HDL-C",
  "hdlc"="HDL-C",
  "tg"="Triglycerides",
  "tc"="Total cholesterol",
  "crp"="CRP"
)

meta_info <- results_lvs_sub %>%
  arrange(desc(beta_sc_adj)) %>%
  dplyr::select(metabolite_name, super_class_metabolon, beta_sc_adj) %>%
  distinct(metabolite_name, .keep_all = TRUE)
mat1_beta_raw <- results_lvs_sub %>%
  dplyr::select(metabolite_name, beta_sc_adj) %>%
  distinct(metabolite_name, .keep_all = TRUE) %>% 
  column_to_rownames("metabolite_name") %>%
  as.matrix()

met_order <- rownames(mat1_beta_raw)[order(mat1_beta_raw[, 1], decreasing = TRUE)]
mat1_beta <- mat1_beta_raw[met_order, , drop = FALSE]
colnames(mat1_beta) <- "DDP"

results_food_renamed <- results_food %>%
  mutate(Food = dplyr::recode(Food, !!!food_group_labels)) %>%
  filter(Food %in% food_group_labels) 
mat2_beta <- results_food_renamed %>%
  dplyr::select(metabolite_name, Food, Beta) %>%
  pivot_wider(names_from = Food, values_from = Beta) %>%
  column_to_rownames("metabolite_name") %>%
  as.matrix()

mat2_p <- results_food_renamed %>%
  dplyr::select(metabolite_name, Food, P_value) %>%
  pivot_wider(names_from = Food, values_from = P_value) %>%
  column_to_rownames("metabolite_name") %>%
  as.matrix()

mat2_beta <- mat2_beta[met_order, ]
mat2_p <- mat2_p[met_order, ]
food_order_original <- coef_df$Original_Name[order(coef_df$Coefficient, decreasing = TRUE)]

food_order_new <- food_group_labels[food_order_original]
food_order_new 
food_order_new <- food_order_new[!is.na(food_order_new)]
valid_food_cols <- intersect(food_order_new, colnames(mat2_beta))

mat2_beta <- mat2_beta[, valid_food_cols]
mat2_p <- mat2_p[, valid_food_cols]

results_bio_renamed <- results_bio 
if(length(biomarker_labels) > 0){
  results_bio_renamed <- results_bio_renamed %>%
    mutate(Biomarker = dplyr::recode(Biomarker, !!!biomarker_labels))
}

mat3_beta <- results_bio_renamed %>%
  dplyr::select(metabolite_name, Biomarker, Beta) %>%
  pivot_wider(names_from = Biomarker, values_from = Beta) %>%
  column_to_rownames("metabolite_name") %>%
  as.matrix()
mat3_beta <- mat3_beta[met_order, ]

mat3_p <- results_bio_renamed %>%
  dplyr::select(metabolite_name, Biomarker, P_value) %>%
  pivot_wider(names_from = Biomarker, values_from = P_value) %>%
  column_to_rownames("metabolite_name") %>%
  as.matrix()
mat3_p <- mat3_p[met_order, ]

max_1 <- max(abs(mat1_beta), na.rm = TRUE)
col_fun_1 = colorRamp2(c(-max_1, 0, max_1), c("navy", "white", "firebrick3"))

max_2 <- max(abs(mat2_beta), na.rm = TRUE)
col_fun_2 = colorRamp2(c(-max_2, 0, max_2), c("#1b7837", "white", "#762a83"))

max_3 <- max(abs(mat3_beta), na.rm = TRUE)
col_fun_3 = colorRamp2(c(-max_3, 0, max_3), c("#b2182b", "white", "#2166ac"))

meta_info_sorted <- results_lvs_sub %>%
  dplyr::select(metabolite_name, super_class_metabolon) %>%
  distinct(metabolite_name, .keep_all = TRUE) %>%
  column_to_rownames("metabolite_name")

sorted_classes <- meta_info_sorted[met_order, "super_class_metabolon"]
sorted_classes <- ifelse(is.na(sorted_classes), "NA", sorted_classes)
unique_classes <- sort(unique(sorted_classes))
n_colors <- length(unique_classes)
unique_classes

my_colors <- brewer.pal(n = n_colors, name = "Set1") 
names(my_colors) <- unique_classes


row_ha = rowAnnotation(
  Class = sorted_classes,
  col = list(Class = my_colors), 
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.5, "cm")
)

cell_fun_stars = function(p_mat) {
  function(j, i, x, y, width, height, fill) {
    p_val = p_mat[i, j]
    if (!is.na(p_val) & p_val < 0.001) {
      grid.text("***", x, y, gp = gpar(fontsize = 8))
    } else if (!is.na(p_val) & p_val < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = 8))
    } else if (!is.na(p_val) & p_val < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 8))
    }
  }
}

ht1 = Heatmap(mat1_beta, 
              name = "Beta (EDID)", 
              col = col_fun_1,
              cluster_rows = FALSE,
              cluster_columns = FALSE, 
              left_annotation = row_ha, 
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 10, fontface = "bold"),
              column_names_rot = 45,
              rect_gp = gpar(col = "white", lwd = 1), 
              width = unit(1.5, "cm"))
ht1

ht2 = Heatmap(mat2_beta, 
              name = "Beta (Food)", 
              col = col_fun_2,
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              show_row_names = FALSE,  
              column_names_rot = 45,
              column_names_gp = gpar(fontsize = 9),
              rect_gp = gpar(col = "white", lwd = 1),
              cell_fun = cell_fun_stars(mat2_p),
              width = unit(10, "cm"))
ht2

ht3 = Heatmap(mat3_beta, 
              name = "Beta (Bio)", 
              col = col_fun_3,
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              show_row_names = FALSE, 
              column_names_rot = 45,
              column_names_gp = gpar(fontsize = 9),
              rect_gp = gpar(col = "white", lwd = 1),
              cell_fun = cell_fun_stars(mat3_p),
              width = unit(6, "cm"))
ht3

ht_list = ht1 + ht2 + ht3
ht_list
draw(ht_list, 
     row_title = "Metabolites",
     merge_legend = FALSE)
pdf("LVS_met.pdf",width = 15,height = 8)
draw(ht_list, 
     row_title = "Metabolites",
     merge_legend = FALSE)
dev.off()
##########Add T2D forest plot##################
load("/udd/n2gji/micro/nut_final/t2d/data/nhs_db_cc.RData")
results<-NULL
run_clogit <- function(data, score_vars, covariates = NULL) {
  results <- data.frame()
  
  for (score_var in score_vars) {
    formula_str <- paste("db ~", score_var)
    if (!is.null(covariates)) {
      formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
    }
    formula_str <- paste(formula_str, "+ strata(matchid) ")
    formula <- as.formula(formula_str)
    
    model_summary <- clogit(formula, data = data) %>%
      broom::tidy() %>%
      filter(grepl(score_var, term)) %>%
      transmute(
        Variable = score_var,
        Level = term,
        OR = exp(estimate),
        Lower_95CI = exp(estimate - 1.96 * std.error),
        Upper_95CI = exp(estimate + 1.96 * std.error),
        P_Value = p.value
      )
    
    results <- bind_rows(results, model_summary)
  }
  
  return(results)
}
met_vars<-as.character(results_lvs_sub$met[results_lvs_sub$met%in% colnames(nhs_db_cohort)])
covariates <- c("pa_metbaseq", "aheibaseq", "smoke", "dbfh","factor(hbcbase)", "factor(htnbase)","factor(bmibase3cat)","factor(calorbaseq)","factor(alcoq)")
results <- run_clogit(nhs_db_cohort,met_vars, covariates)
results <- run_clogit(nhs_db_cohort, met_vars, covariates) %>%
  mutate(FDR_P = p.adjust(P_Value, method = "BH"))

id_name_map <- results_lvs_sub %>% 
  dplyr::select(met, metabolite_name) %>% 
  distinct()

results_t2d_mapped <- results %>%
  left_join(id_name_map, by = c("Variable" = "met"))

forest_data <- data.frame(metabolite_name = met_order) %>%
  left_join(results_t2d_mapped, by = "metabolite_name")
vals <- c(forest_data$Lower_95CI, forest_data$Upper_95CI)
vals <- vals[!is.na(vals) & is.finite(vals)]

if(length(vals) > 0) {
  x_min <- floor(min(vals) * 10) / 10
  x_max <- ceiling(max(vals) * 10) / 10
  if(x_min > 0.9) x_min <- 0.9
  if(x_max < 1.1) x_max <- 1.1
  x_limits <- c(x_min, x_max)
} else {
  x_limits <- c(0.5, 1.5) 
}

forest_ha <- rowAnnotation(
  "T2D OR (95% CI)" = anno_empty(border = FALSE, width = unit(3, "cm")),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
)

final_ht_list <- ht1 + ht2 + ht3 + forest_ha
final_ht_list
pdf("lvs_plot_final.pdf", width = 18, height = 15)
final_ht_list <- ht1 + ht2 + ht3 + forest_ha

draw(final_ht_list, 
     row_title = "Metabolites",
     merge_legend = FALSE)

decorate_annotation("T2D OR (95% CI)", {
  
  n_rows <- nrow(forest_data)
  
  pushViewport(viewport(xscale = x_limits, yscale = c(0.5, n_rows + 0.5)))
  
  grid.lines(x = unit(1, "native"), 
             y = unit(c(0, n_rows + 1), "native"), 
             gp = gpar(lty = 2, col = "grey70", lwd = 1))
  
  grid.xaxis(at = c(x_limits[1], 1, x_limits[2]), 
             gp = gpar(fontsize = 7))
  
  for (i in 1:n_rows) {
    row_dat <- forest_data[i, ]
    
    if (!is.na(row_dat$OR) && !is.na(row_dat$Lower_95CI) && !is.na(row_dat$FDR_P)) {
      
      pt_col <- "black"
      
      if (row_dat$FDR_P < 0.05) {
        if (row_dat$OR > 1) {
          pt_col <- "indianred"
        } else {
          pt_col <- "steelblue"
        }
      }
      
      y_pos <- n_rows - i + 1
      
      grid.segments(x0 = unit(row_dat$Lower_95CI, "native"), 
                    x1 = unit(row_dat$Upper_95CI, "native"), 
                    y0 = unit(y_pos, "native"), 
                    y1 = unit(y_pos, "native"), 
                    gp = gpar(col = pt_col, lwd = 1))
      
      tick_h <- 0.12 
      
      grid.segments(x0 = unit(row_dat$Lower_95CI, "native"), 
                    x1 = unit(row_dat$Lower_95CI, "native"),
                    y0 = unit(y_pos - tick_h, "native"), 
                    y1 = unit(y_pos + tick_h, "native"),
                    gp = gpar(col = pt_col, lwd = 1))
      
      grid.segments(x0 = unit(row_dat$Upper_95CI, "native"), 
                    x1 = unit(row_dat$Upper_95CI, "native"),
                    y0 = unit(y_pos - tick_h, "native"), 
                    y1 = unit(y_pos + tick_h, "native"),
                    gp = gpar(col = pt_col, lwd = 1))
      
      grid.points(x = unit(row_dat$OR, "native"), 
                  y = unit(y_pos, "native"), 
                  pch = 15,
                  size = unit(2, "mm"), 
                  gp = gpar(col = pt_col))
    }
  }
  
  popViewport()
})
dev.off()
###########Build metabolite score in LVS##############
set.seed(123)
dt_nut<-temp_lvs
#-------------------------------------------------------------------------------
# 					                         Training
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# NOTE 2: 10-fold cross validation to obtain the optimal lambda
#         Leave-one-out approach to calculate the signature
# Use cv.glmnet
#-------------------------------------------------------------------------------#
temp_lvs$score_s
met_list = as.character(colnames(temp_lvs[, 3:297]))
met_list<-met_list[met_list%in%colnames(nhs_db_cohort)]#228 in 295
stage1_mx1 <- dt_nut %>%  dplyr::select(all_of(met_list),score_s)
stage1a <- dt_nut %>%dplyr::select(id)
record = data.frame(accumulateid=NA,newid=NA,NoMetabs=NA,accumulate_cor=NA)

selected_metabolites_list = vector("list", nrow(stage1_mx1))
coefficients_data <- data.frame(matrix(NA, nrow = nrow(dt_nut), ncol = 229))
colnames(coefficients_data) <- c("intercept", colnames(stage1_mx1)[1:228])

set.seed(123)

for (i in 1:dim(stage1_mx1)[1]) {
  
  Training = cv.glmnet(as.matrix(stage1_mx1[-i,c(1:228)]), stage1_mx1[-i,"score_s"], nfolds=10, alpha=0.5, family="gaussian")
  lambda_min_10F = Training$lambda.min
  
  cmin = coef(Training, s=lambda_min_10F)
  
  coefficients_data[i, ] <- cmin[, 1]
  
  metabsmin = data.frame(cmin[which(cmin[,1]!=0),])
  names(metabsmin) = "Test.min"
  
  selected_metabolites = rownames(metabsmin)[which(metabsmin$Test.min != 0)]
  
  selected_metabolites_list[[i]] <- selected_metabolites
  
  if(i==1) {
    stage1_mx1[i,dim(stage1_mx1)[2]+1] = predict(Training, as.matrix(stage1_mx1[i,c(1:228)]), type="response", s=lambda_min_10F)
  } else if(i>1) {
    stage1_mx1[i,dim(stage1_mx1)[2]] = predict(Training, as.matrix(stage1_mx1[i,c(1:228)]), type="response", s=lambda_min_10F)
  }
  
  record[i,"accumulateid"] = i
  record[i,"newid"] = stage1a[i,"id"]
  record[i,"NoMetabs"] = dim(metabsmin)[1]-1
  
  if(i<3) {
    record[i,"accumulate_cor"] = NA
  } else if(i>=3) {
    tmp = stage1_mx1[1:i,]
    record[i,"accumulate_cor"] = cor(tmp[,"score_s"], tmp[,dim(tmp)[2]])
  }
}
coef_mean <- as.data.frame(colMeans(coefficients_data, na.rm = T))
coef_mean <- cbind(hmdb_id = rownames(coef_mean), coef_mean)
save(coef_mean,selected_metabolites_list,file="./met_score.RData")
all_selected = unlist(selected_metabolites_list)
metabolite_frequencies = table(all_selected)
sorted_metabolites = as.data.frame(sort(metabolite_frequencies, decreasing = TRUE))
colnames(sorted_metabolites) <- c("hmdb_id", "Freq") 
total_iters <- length(selected_metabolites_list) 
freq_threshold <- total_iters * 0.90 
metab_freq <- subset(sorted_metabolites, Freq >= freq_threshold & hmdb_id != "intercept" & hmdb_id != "(Intercept)")
coef_mean2 <- subset(coef_mean, hmdb_id != "intercept" & hmdb_id != "(Intercept)")
metab_sig <- inner_join(metab_freq, coef_mean2, by = "hmdb_id")
metab_sig$intercept <- coef_mean[coef_mean$hmdb_id %in% c("intercept", "(Intercept)"), 2]
metab_sig2 <- merge(metab_sig, f_lvs_match, by = "hmdb_id",by.y = "met")
# 接下去你原本的 write.csv 等后续代码...
write.csv(metab_sig2, file="./met_score")
metab_sig2<-read.csv("./met_score")
#-------------------------------------------------------------------------------
# 					                         Testing
#-------------------------------------------------------------------------------
# Perform linear regression to predict dietary intake (y) based on the score

# First 39 metabolites
stage1_metab <- stage1_mx1 %>% dplyr::select(metab_sig2$hmdb_id)
score1<- stage1_metab %>% as.matrix() %*% metab_sig2$colMeans.coefficients_data..na.rm...T.
intercept <- metab_sig2$intercept
score <- score1 + intercept
stage1_mx1$score <- score


print(length(score))
print(length(dt_nut))

dt_nut$score <- score
dt_nut$score_met_s<-scale(dt_nut$score)

cor.test(dt_nut$dietary_score, dt_nut$score_met_s, method="pearson", conf.level=0.05)#0.71
cor_test <- cor.test(dt_nut$dietary_score, dt_nut$score_met_s, method="spearman", conf.level=0.05)#0.59

r_value <- formatC(cor_test$estimate, digits = 2, format = "f")
p_value <- ifelse(cor_test$p.value < 0.001, 
                  "< 0.001", 
                  formatC(cor_test$p.value, digits = 3, format = "f"))
library(ggpubr)

p <- ggplot(dt_nut, aes(x = score_met_s, y = dietary_score)) +
  
  geom_point(
    shape = 21,
    size = 3,
    fill = "#4E79A7",
    color = "white",
    alpha = 0.8,
    stroke = 0.6
  ) +
  
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    color = "#E15759",
    fill = "#F28E2B",
    alpha = 0.2,
    linewidth = 1.2,
    se = TRUE
  ) +
  
  annotate(
    "text",
    x = min(dt_nut$score_met_s),
    y = max(dt_nut$dietary_score),
    label = paste0("r = ", r_value, "\np ", p_value),
    hjust = 0, vjust = 1,
    color = "#444444",
    size = 4.8,
    fontface = "bold"
  ) +
  
  labs(
    x = "Rank Score",
    y = "Dietary Score",
  ) +
  theme_pubr() +
  theme(
    plot.background = element_rect(fill = "#F5F7F9", color = NA),
    panel.background = element_rect(fill = "white", color = "#EAEAEA", linewidth = 1),
    plot.title = element_text(face = "bold", size = 16, color = "#2E4660"),
    plot.subtitle = element_text(color = "#555555", size = 11),
    plot.caption = element_text(color = "#888888"),
    axis.title = element_text(color = "#444444", face = "bold"),
    panel.grid = element_line(color = "#EEEEEE")
    
  )

print(p)
ggsave("scatterplot_metscore.png", plot = p, width = 10, height = 8, dpi = 300, bg = "white")
names(metab_sig2)

data_plot <- metab_sig2
names(data_plot)[names(data_plot) == "colMeans.coefficients_data..na.rm...T."] <- "Mean_Coefficient"
data_plot$super_class_metabolon<-ifelse(is.na(data_plot$super_class_metabolon),"NA",data_plot$super_class_metabolon)
table(data_plot$super_class_metabolon,useNA="always")

p<-ggplot(data_plot, aes(x = Mean_Coefficient, 
                         y = reorder(metabolite_name, Mean_Coefficient), 
                         fill = super_class_metabolon)) +
  geom_col(width = 0.7) +
  
  labs(x = "Mean Coefficients", 
       y = "Metabolite Name", 
       fill = "Class type") +
  
  theme_classic() +
  
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90", linetype = "dashed")
  ) +
  
  scale_fill_brewer(palette = "Set3") 
p
ggsave("coef_metscore.png", plot = p, width = 10, height = 8, dpi = 300, bg = "white")

names(dt_nut)
x_val <- dt_nut[, c("score_met_s"), drop = FALSE]
class(x_val)
dt_nut$dietary_score
y_val <- dt_nut[, c("hba1c","tg","tc","hdlc","crp", "tg.hdlc")]

cor_matrix <- cor(y_val, x_val, use = "pairwise.complete.obs")

p_matrix <- matrix(NA, nrow = ncol(y_val), ncol = ncol(x_val),
                   dimnames = list(colnames(y_val), colnames(x_val)))

for (row_name in rownames(p_matrix)) {
  for (col_name in colnames(p_matrix)) {
    complete_obs_count <- sum(complete.cases(y_val[[row_name]], x_val[[col_name]]))
    if (complete_obs_count > 2) {
      test_result <- cor.test(y_val[[row_name]], x_val[[col_name]], use = "pairwise.complete.obs")
      p_matrix[row_name, col_name] <- test_result$p.value
    } 
  }
}

plot_data <- data.frame(
  Variable = rownames(cor_matrix),
  Correlation = cor_matrix[, 1],
  P_Value = p_matrix[, 1]
)

plot_data <- plot_data %>%
  mutate(
    Significance = case_when(
      is.na(P_Value) ~ "",
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    Direction = ifelse(Correlation > 0, "Positive", "Negative"),
    
    Variable_Full = recode(
      Variable,
      hba1c   = "HbA1c",
      tg      = "Triglycerides",
      tc      = "Total Cholesterol",
      hdlc    = "HDL-C",
      crp     = "CRP",
      tg.hdlc = "TG/HDL-C Ratio",
      .default = Variable
    )
  )

print(plot_data)

p<-ggplot(plot_data, aes(x = reorder(Variable_Full, Correlation), y = Correlation, fill = Direction)) +
  
  geom_col(width = 0.7) + 
  
  geom_text(
    aes(label = Significance, y = Correlation + 0.02 * sign(Correlation)), 
    vjust = ifelse(plot_data$Correlation > 0, 0, 1),
    size = 5,
    color = "black"
  ) +
  
  coord_flip() + 
  
  scale_fill_manual(values = c("Positive" = "#d6604d", "Negative" = "#4393c3")) +
  
  labs(
    x = "Biomarker / Diet Score",
    y = "Pearson Correlation Coefficient"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  ylim(min(plot_data$Correlation) - 0.05, max(plot_data$Correlation) + 0.05)
p
ggsave("metscore_biomarker.png", plot = p, width = 10, height = 8, dpi = 300, bg = "white")