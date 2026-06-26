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
library(ggpubr)
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
library(Maaslin2)
var<-NULL
load("/udd/n2gji/micro/data/ToBeUsed.Taxon.RData")
load(file="/udd/n2gji/micro/nut_final/data/Final_Data_2nd.RData")
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
  
  # 1. 获取目标变量名 (例如: marg_serv_log_z) 和 权重
  target_var <- coef_df$Variable[i]
  weight <- coef_df$Coefficient[i]
  
  # 2. 去掉 "_log_z" 后缀，得到中间名 (例如: marg_serv)
  base_name <- sub("_log_z", "", target_var)
  
  # 3. 在 name_mapping 中查找对应的原始列名 (例如: marg07_serv)
  # name_mapping 的结构是: 原始名 = 中间名
  original_col_name <- names(name_mapping)[which(name_mapping == base_name)]
  
  # 4. 检查该列是否存在于数据库中
  if(length(original_col_name) > 0 && original_col_name %in% names(temp_lvs)) {
    
    # 提取原始数据
    raw_data <- temp_lvs[[original_col_name]]
    
    # --- 关键步骤：数据转换 ---
    # 因为系数对应的是 log_z，我们需要对原始数据做同样的转换
    # 1. 对数转换 (加一个小常数避免 log(0))
    log_data <- log(raw_data + 1)
    
    # 2. Z-score 标准化 (减均值除以标准差)
    # 注意：scale() 返回的是矩阵，需要转回 numeric
    z_data <- as.numeric(scale(log_data))
    
    # 5. 累加分数 (标准化后的值 * 权重)
    # 处理可能的 NA 值 (如果某行为 NA，则结果为 NA，或者你可以设为 0)
    term_value <- z_data * weight
    
    # 将计算结果加到总分中
    temp_lvs$dietary_score <- temp_lvs$dietary_score + term_value
    
    print(paste("已计算:", base_name, "-> 原始列:", original_col_name, "| 权重:", weight))
    
  } else {
    print(paste("⚠️ 警告: 未在 temp_lvs 或 mapping 中找到对应列:", base_name))
    # 注意：你的映射表中似乎缺少 "rcer_serv" (Refined Cereal)，这可能会触发警告
  }
}
summary(temp_lvs$dietary_score)
temp_lvs$score_s<-scale(temp_lvs$dietary_score)
temp_mlvs<-subset(temp_lvs,temp_lvs$cohort=="hpfs")
ToBeUsed.Taxon<-merge(ToBeUsed.Taxon,temp_mlvs[,c("id", "dietary_score","score_s",
                                                  "procm07","rmeat07","liver07","seaf07","sseaf07","poul07","eggs07","butter07",
                                                  "marg07","milk07","yogurt07","cheese07","daiswt07","wine07","beer07_serv","liq07_serv","tea07",
                                                  "coff07","h2o07","sugbev07","lowbev07","citf07","fruj07","berr07","appear07","ban07",
                                                  "rais07","peach07","othfr07","lfveg07","cruveg07","stachveg07",
                                                  "nightveg07","legu07","othveg07","tompro07","apiveg07",
                                                  "fries07","wgrain07","rgrain07","snack07","pnuts07","otnuts07",
                                                  "choco07","swt07","dess07","cond07","saldre07","crsoup07","pizza07","cokh07","cokr07")],by="id",all.x = T)
ToBeUsed.Taxon<-as.data.frame(ToBeUsed.Taxon)
Microbiome <- grep("^s_", colnames(ToBeUsed.Taxon), value = TRUE)
Phenotypes = setdiff(names(ToBeUsed.Taxon),Microbiome)
rownames(ToBeUsed.Taxon) <- ToBeUsed.Taxon$id
AdjVars = c("ageyr","totMETs_paq","calor_fo_dr_wtavg","smoke_bld","alco_bld",
            "probio_2m_fec","antibio_12m_fec","colsc_2m_fec","acid_2m_fec",
            "stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4","stooltype_fec.5","stooltype_fec.6")#Note: not adjusing for BMI
Maaslin2(input_data=ToBeUsed.Taxon[,Microbiome],
         input_metadata=ToBeUsed.Taxon[,Phenotypes],
         output=paste("/udd/n2gji/micro/dietpt/final/msl"), 
         min_abundance = 0.0001,
         min_prevalence = 0.1, 
         normalization = "NONE",  #Should set to NONE because species data were already normalized by TSS
         transform = "AST", 
         analysis_method = "LM",
         # random_effects = "SampleID", 
         #No random effect here because there was one-time measurement of metabolomic data and all 4 stool measurements were averaged
         fixed_effects = c("dietary_score",AdjVars), 
         correction = "BH",
         standardize = TRUE, 
         cores = 8,
         plot_heatmap=FALSE,
         plot_scatter=FALSE)
msl<-read_tsv("/udd/n2gji/micro/dietpt/final/msl/all_results.tsv")
msl<-subset(msl,msl$value=="dietary_score")
msl<-merge(msl,Annot_Taxon,by.x = "feature",by.y = "taxon_new",all.x = T )
msl_sub<-subset(msl,msl$qval<0.20)
write.csv(msl_sub,file = "/udd/n2gji/micro/dietpt/db2/output/msl/msl_sub.csv")
msl_sub<-read.csv("/udd/n2gji/micro/dietpt/db2/output/msl/msl_sub.csv")
########Rank score##########
sig_dta <- data.frame(
  feature = msl_sub$feature,
  coef = msl_sub$coef
)
rank_score <- function(dta, sig_dta, varname, feature_col = "feature") {
  # sig_dta must contain: <feature_col>, <varname> (e.g., "coef")
  if (!feature_col %in% names(sig_dta))
    stop("feature_col '", feature_col, "' not found in sig_dta. Have: ", paste(names(sig_dta), collapse = ", "))
  if (!varname %in% names(sig_dta))
    stop("varname '", varname, "' not found in sig_dta. Have: ", paste(names(sig_dta), collapse = ", "))
  
  # compact table of feature ids + weights
  sig_tbl <- data.frame(
    feature = sig_dta[[feature_col]],
    weight  = sig_dta[[varname]],
    stringsAsFactors = FALSE
  )
  sig_tbl <- sig_tbl[!is.na(sig_tbl$feature) & !is.na(sig_tbl$weight), , drop = FALSE]
  sig_tbl <- sig_tbl[!duplicated(sig_tbl$feature), , drop = FALSE]
  
  feats <- intersect(sig_tbl$feature, colnames(dta))
  if (length(feats) == 0L) {
    dta$rank.score <- 0
    dta$rank.score.binary <- 0
    dta$positive.rank <- 0
    dta$negative.rank <- 0
    dta$coef.weighted.score <- 0
    dta$coef.weighted.score.binary <- 0
    return(dta)
  }
  sig_tbl <- sig_tbl[match(feats, sig_tbl$feature), , drop = FALSE]
  
  # feature matrix
  X <- as.matrix(dta[, feats, drop = FALSE])
  storage.mode(X) <- "double"
  
  # per-feature ranks with zero-handling (zeros/NA -> 0; else rank - n_zeros)
  R <- apply(X, 2, function(col) {
    r  <- suppressWarnings(rank(col, ties.method = "average", na.last = "keep"))
    zc <- sum(col == 0, na.rm = TRUE)
    ifelse(is.na(col) | col == 0, 0, r - zc)
  })
  if (is.null(dim(R))) {        # keep matrix shape when only 1 feature
    dim(R) <- c(nrow(dta), 1L)
    colnames(R) <- feats
    rownames(R) <- rownames(dta)
  }
  
  w   <- sig_tbl$weight
  names(w) <- sig_tbl$feature
  w   <- w[colnames(R)]
  sgn <- sign(w)
  
  # scores
  rank.score              <- as.vector(R %*% sgn)   # unweighted signed rank sum
  coef.weighted.score     <- as.vector(R %*% w)     # coefficient-weighted rank sum
  contrib                 <- sweep(R, 2, sgn, `*`)
  positive.rank           <- if (any(sgn > 0)) rowSums(contrib[, sgn > 0, drop = FALSE]) else rep(0, nrow(dta))
  negative.rank           <- if (any(sgn < 0)) rowSums(contrib[, sgn < 0, drop = FALSE]) else rep(0, nrow(dta))
  rank.score.binary       <- as.integer(rank.score > stats::median(rank.score, na.rm = TRUE))
  coef.weighted.binary    <- as.integer(coef.weighted.score > stats::median(coef.weighted.score, na.rm = TRUE))
  
  # attach to dta
  dta$rank.score                 <- rank.score
  dta$rank.score.binary          <- rank.score.binary
  dta$positive.rank              <- positive.rank
  dta$negative.rank              <- negative.rank
  dta$coef.weighted.score        <- coef.weighted.score
  dta$coef.weighted.score.binary <- coef.weighted.binary
  dta
}
ToBeUsed.Taxon <- rank_score(ToBeUsed.Taxon, sig_dta, varname = "coef", feature_col = "feature")
cor.test(ToBeUsed.Taxon$rank.score,ToBeUsed.Taxon$dietary_score,method = "spearman")

cor_test <- cor.test(ToBeUsed.Taxon$rank.score, ToBeUsed.Taxon$dietary_score,method = "pearson")
r_value <- formatC(cor_test$estimate, digits = 2, format = "f")
p_value <- ifelse(cor_test$p.value < 0.001, 
                  "< 0.001", 
                  formatC(cor_test$p.value, digits = 3, format = "f"))
p <- ggplot(ToBeUsed.Taxon, aes(x = rank.score, y = dietary_score)) +
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
    x = min(ToBeUsed.Taxon$rank.score),
    y = max(ToBeUsed.Taxon$dietary_score),
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
ggsave("scatterplot_rankscore_notadjbmi.png", plot = p, width = 10, height = 8, dpi = 300, bg = "white")
msl_sub$feature
sgb<-read.csv("/udd/hpxil/3.source/SGB_annotation.csv")
msl_sub$new_name<-c("Negativibacillus sp000435195","GGB52130 SGB14966","UBA11774 sp003507655","Roseburia Inulinivorans",
                    "Butyrivibrio Crossotus",
                    "Ruminococcus Torques",
                    "Lawsonibacter Asaccharolyticus")
###############################Micro N############################
load("/udd/n2gji/micro/data/ReadIn_microN_taxon.RData")
dim( MicroN_Prescient_Microb_s_avg) # 897 2356
dim( MicroN_B2B_Microb_s_avg) # 844 2275
dim( MicroN_T2D_Microb_s_avg) # 233 1699
dim( MicroN_CVD_Microb_s_avg) # 152 1586
dim( MicroN_GDM_Microb_s_avg) # 610 2109
dim( MicroN_KidneyStone_Microb_s_avg) # 914 2346
dim( MicroN_Diverticulitis_lab20200_Microb_s_avg) # 236 1815
dim( MicroN_Diverticulitis_lab20246_20247_Microb_s_avg) # 343 1966
dim( MicroN_inc_polyp_batch1_Microb_s_avg) # 327 1916
dim( MicroN_inc_polyp_batch2_Microb_s_avg) # 506 1970
dim( MicroN_Cognitive_Microb_s_avg) # 1582 2542
get_id_and_source <- function(df, source_name) {
  df_out <- df
  if ("id" %in% colnames(df_out)) {
    
  } else {
    
    df_out <- df_out %>% rownames_to_column(var = 'id')
  }
  df_out %>% 
    mutate(source = source_name) %>% 
    dplyr::select(id, source) %>% 
    
    mutate(id = as.character(id))
}
MicroN_Pool_id <- bind_rows(
  get_id_and_source(MicroN_Prescient_Microb_s_avg,               'MicroN_Prescient'),
  get_id_and_source(MicroN_B2B_Microb_s_avg,                     'MicroN_B2B'),
  get_id_and_source(MicroN_T2D_Microb_s_avg,                     'MicroN_T2D'),
  get_id_and_source(MicroN_CVD_Microb_s_avg,                     'MicroN_CVD'),
  get_id_and_source(MicroN_KidneyStone_Microb_s_avg,             'MicroN_KidneyStone'),
  get_id_and_source(MicroN_GDM_Microb_s_avg,                     'MicroN_GDM'),
  get_id_and_source(MicroN_Diverticulitis_lab20200_Microb_s_avg, 'MicroN_Diverticulitis_lab20200'),
  get_id_and_source(MicroN_Diverticulitis_lab20246_20247_Microb_s_avg, 'MicroN_Diverticulitis_lab20246_20247'),
  get_id_and_source(MicroN_inc_polyp_batch1_Microb_s_avg,        'MicroN_inc_polyp_batch1'),
  get_id_and_source(MicroN_inc_polyp_batch2_Microb_s_avg,        'MicroN_inc_polyp_batch2'),
  get_id_and_source(MicroN_Cognitive_Microb_s_avg,               'MicroN_Cognitive')
)
dim(MicroN_Pool_id)
# 6644    2

MicroN_Dup_id <- MicroN_Pool_id %>% 
  dplyr::count(id,name='n') %>% 
  filter(!is.na(id),n>1) %>% 
  pull(id)
MicroN_Dup_id # n=714
MicroN_Dup_id_df <- MicroN_Pool_id %>% 
  filter(id%in%MicroN_Dup_id)
# n=1035

# ------------------- Rules for Duplicated Samples -----------------------------

MicroN_Pool_id_unique <- MicroN_Pool_id %>% 
  distinct(id) %>% 
  pull(id)

length(MicroN_Pool_id_unique)
# 6116 --> unique samples
# ------------------- Link with Fecal Sample Qx --------------------------------
micron_t2d <-read.csv("/udd/n2gji/data/nhs2_t2d.csv")
stool_cov<-read.csv("/proj/n2dats/n2dat02/subStudy/MicroN/stool_collection_qx/nhs2_micron.csv")
stool_cov<-stool_cov %>%
  dplyr::filter (stool_cov$id %in% MicroN_Pool_id$id)
micron_t2d<-merge(stool_cov,micron_t2d,by="id",all.x = T)

library(lubridate)
library(dplyr)

# Month map (locale-independent)
.mon_map <- setNames(1:12, tolower(month.abb))
.mon_map 
# .mon_map
# jan feb mar apr may jun jul aug sep oct nov dec 
# 1   2   3   4   5   6   7   8   9  10  11  12 
# Convert "mon yy" or a numeric string like "1434" to months since Jan 1900 (Jan 1900 = 0)
to_months1900 <- function(x, cutoff = 70, min_ok = 0L, max_ok = 3000L) {
  x_clean <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(x_clean))
  
  # Case 1: already numeric months since 1900 (e.g., "1434")
  is_num <- grepl("^[0-9]{3,4}$", x_clean)
  if (any(is_num)) {
    val <- suppressWarnings(as.integer(x_clean[is_num]))
    # sanity range (adjust max_ok if you need later years)
    val[is.na(val) | val < min_ok | val > max_ok] <- NA_integer_
    out[is_num] <- val
  }
  
  # Case 2: "mon yy" like "apr 20"
  is_mon_yy <- grepl("^[a-z]{3}\\s+\\d{2}$", x_clean)
  if (any(is_mon_yy)) {
    mon3 <- sub("^([a-z]{3}).*$", "\\1", x_clean[is_mon_yy])
    yy2  <- as.integer(sub("^[a-z]{3}\\s+(\\d{2})$", "\\1", x_clean[is_mon_yy]))
    
    mon_num <- unname(.mon_map[mon3])
    # Expand 2-digit year: 70–99 -> 1900s; 00–69 -> 2000s
    year <- ifelse(yy2 >= cutoff, 1900L + yy2, 2000L + yy2)
    
    months <- (year - 1900L) * 12L + (mon_num - 1L)
    # Keep only plausible range
    months[is.na(mon_num) | is.na(yy2) | months < min_ok | months > max_ok] <- NA_integer_
    out[is_mon_yy] <- months
  }
  
  out
}

# Apply to the derived MicroN data frame
micron_t2d <- micron_t2d %>%
  mutate(retmo19_clean = to_months1900(retmo19))

# Optional quick checks:
table(micron_t2d$retmo19, useNA = "ifany")
summary(micron_t2d$retmo19_clean)
micron_t2d$retmo19_clean<-ifelse(is.na(micron_t2d$retmo19_clean),1441,micron_t2d$retmo19_clean)##Using median to replace missing value
#apr 20 apr 21 aug 20 dec 19 feb 20 feb 21 jan 20 jan 21 jul 20 jun 20 mar 20 mar 21 may 20 may 21 nov 19 sep 20   <NA> 
# 12     76      3     13    238   1218     27   1562     38      2      7    398      4     27      8      3     12     73 
# summary(micron_t2d$retmo19_clean)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1438    1440    1441    1441    1441    1456 


priority_vars<- c("yr","dtdxdb","dtdxdb2","ninsu17","metf17","jard17","invk17", "sita17","ohypo17","insul17",
                  "insul19","metf19","jard19", "invk19","farx19","sita19","ohypo19","ninsu19")

summary(micron_t2d$yr)
MicroN_IncidentT2D_nonExc <- micron_t2d %>% 
  dplyr::filter(id%in%MicroN_Pool_id_unique) %>% 
  #filter(!chrons19==1) %>% # limited to no baseline disease 
  # treat sentinel as missing if applicable
  dplyr::mutate(dtdxdb2 = na_if(dtdxdb2, 9999),
                yr_num = suppressWarnings(as.numeric(yr))) %>%
  # priority score = number of non-missing fields among your variables
  dplyr::mutate(nonmiss = rowSums(across(all_of(priority_vars), ~ !is.na(.))),
                tie_confirmed = !is.na(dtdxdb2),   # prefer record with confirmed onset date
                tie_sr       = !is.na(dtdxdb)) %>% # then prefer record with SR onset date
  dplyr::arrange(id,
                 desc(nonmiss),     # 1) most complete record
                 desc(yr_num),      # 2) latest year (if yr is 1..5)
                 desc(tie_confirmed),
                 desc(tie_sr)) %>%
  distinct(id, .keep_all = TRUE) %>%
  dplyr::select(-nonmiss, -yr_num, -tie_confirmed, -tie_sr) %>% #3721
  dplyr::mutate(
    fecalmonths=case_when(
      yr==1 ~ (2018-1900)*12+mo,
      yr==2 ~ (2019-1900)*12+mo,
      yr==3 ~ (2020-1900)*12+mo,
      yr==4 ~ (2021-1900)*12+mo,
      yr==5 ~ (2022-1900)*12+mo,
      TRUE ~ NA_real_),
    
    fecalyears = case_when(
      yr == 1 ~ 2018,
      yr == 2 ~ 2019,
      yr == 3 ~ 2020,
      yr == 4 ~ 2021,
      yr == 5 ~ 2022,
      TRUE    ~ NA_real_ ),
    
    # turn 9999 into NA (cleaner than case_when)
    dbcase_sr = case_when(
      !is.na(dtdxdb) ~ 1,
      TRUE ~ 0
    ),
    dtdxdb2_clean = na_if(dtdxdb2, 9999),
    dbcase_confirmed = case_when(
      !is.na(dtdxdb2_clean) ~ 1,
      TRUE ~ 0
    )) %>%
  dplyr::mutate(
    metf17=case_when(
      metf17=='yes'~ 1,
      TRUE ~ 0
    ),
    insul17=case_when(
      insul17=='yes'~ 1,
      TRUE ~ 0
    ),
    jard17=case_when(
      jard17=='yes'~ 1,
      TRUE ~ 0
    ), 
    invk17=case_when(
      invk17=='yes'~ 1,
      TRUE ~ 0
    ),
    sita17=case_when(
      sita17=='yes'~ 1,
      TRUE ~ 0
    ),
    ohypo17=case_when(
      ohypo17=='yes'~ 1,
      TRUE ~ 0
    ),
    ninsu17=case_when(
      ninsu17=='yes'~ 1,
      TRUE ~ 0
    ),
    metf19=case_when(
      metf19=='yes'~ 1,
      TRUE ~ 0
    ),
    insul19=case_when(
      insul19=='yes'~ 1,
      TRUE ~ 0
    ),
    jard19=case_when(
      jard19=='yes'~ 1,
      TRUE ~ 0
    ),
    invk19=case_when(
      invk19=='yes'~ 1,
      TRUE ~ 0
    ),
    sita19=case_when(
      sita19=='yes'~ 1,
      TRUE ~ 0
    ),
    farx19=case_when(
      farx19=='yes'~ 1,
      TRUE ~ 0),
    ninsu19=case_when(
      ninsu19=='yes'~ 1,
      TRUE ~ 0),
    ohypo19=case_when(
      ohypo19=='yes'~ 1,
      TRUE ~ 0
    )) %>% 
  dplyr::mutate(
    # classify relative to diagnosis time
    fecalqx_T2D_BeforeAfter = case_when(
      dbcase_confirmed  == 1 & !is.na(dtdxdb2_clean) & !is.na(fecalmonths) & fecalmonths >= dtdxdb2_clean ~ "oldT2Dconfirmed",
      dbcase_confirmed  == 1 & !is.na(dtdxdb2_clean) & !is.na(fecalmonths) & fecalmonths <  dtdxdb2_clean ~ "newT2Dconfirmed",
      dbcase_confirmed  == 1 &  is.na(dtdxdb2_clean)                                                    ~ "case_confirmed_unknown_onset",
      dbcase_confirmed  == 0                                                                             ~ "non_cases",
      TRUE                                                                                    ~ NA_character_
    ),
    fecalqx_T2D_BeforeAfter_Yr = case_when(
      dbcase_confirmed  == 1 & !is.na(dtdxdb2_clean) & !is.na(fecalyears) & fecalyears >= (dtdxdb2_clean/12 + 1900)+1  ~ "oldT2Dconfirmed",
      dbcase_confirmed  == 1 & !is.na(dtdxdb2_clean) & !is.na(fecalyears) & fecalyears <  (dtdxdb2_clean/12 + 1900)+1  ~ "newT2Dconfirmed",
      dbcase_confirmed  == 1 &  is.na(dtdxdb2_clean)                                                    ~ "case_confirmed_unknown_onset",
      dbcase_confirmed  == 0                                                                             ~ "non_cases",
      TRUE                                                                                    ~ NA_character_
    ),
    fecalqx_srT2D_BeforeAfter = case_when(
      dbcase_sr == 1 & !is.na(dtdxdb) & !is.na(fecalmonths) & fecalmonths >= dtdxdb ~ "oldT2Dsr",
      dbcase_sr == 1 & !is.na(dtdxdb) & !is.na(fecalmonths) & fecalmonths <  dtdxdb  ~ "newT2Dsr",
      dbcase_sr == 1 &  is.na(dtdxdb)                                               ~ "case_sr_unknown_onset",
      dbcase_sr == 0                                                                ~ "non_cases",
      TRUE                                                                          ~ NA_character_
    ),
    fecalqx_srT2D_BeforeAfter_Yr = case_when(
      dbcase_sr == 1 & !is.na(dtdxdb) & !is.na(fecalyears) & fecalyears >= (dtdxdb/12 + 1900)+1   ~ "oldT2Dsr",
      dbcase_sr == 1 & !is.na(dtdxdb) & !is.na(fecalyears) & fecalyears <  (dtdxdb/12 + 1900)+1  ~ "newT2Dsr",
      dbcase_sr == 1 &  is.na(dtdxdb)                                               ~ "case_sr_unknown_onset",
      dbcase_sr == 0                                                                ~ "non_cases",
      TRUE                                                                          ~ NA_character_
    )
  ) %>% 
  dplyr::mutate(
    fecalqx_srT2D_BeforeAfter_Pool = case_when (
      fecalqx_srT2D_BeforeAfter == 'oldT2Dsr' | (metf17+insul17+jard17+invk17+sita17+ohypo17+ninsu17>0) ~ 'oldT2Dsr',
      fecalqx_srT2D_BeforeAfter == 'newT2Dsr' | (metf17+insul17+jard17+invk17+sita17+ohypo17+ninsu17==0)&(metf19+insul19+jard19+sita19+farx19+invk19+ohypo19+ninsu19>0) ~ 'newT2Dsr',
      fecalqx_srT2D_BeforeAfter == 'case_sr_unknown_onset' ~ 'case_sr_unknown_onset',
      fecalqx_srT2D_BeforeAfter == 'non_cases' ~ 'non_cases',
      TRUE  ~ NA_character_
    ),
    fecalqx_srT2D_BeforeAfter_Yr_Pool = case_when (
      fecalqx_srT2D_BeforeAfter_Yr == 'oldT2Dsr' | (metf17+insul17+jard17+invk17+sita17+ohypo17+ninsu17>0) ~ 'oldT2Dsr',
      fecalqx_srT2D_BeforeAfter_Yr == 'newT2Dsr' | (metf17+insul17+jard17+invk17+sita17+ohypo17+ninsu17==0)&(metf19+insul19+jard19+sita19+farx19+invk19+ohypo19+ninsu19>0) ~ 'newT2Dsr',
      fecalqx_srT2D_BeforeAfter_Yr == 'case_sr_unknown_onset' ~ 'case_sr_unknown_onset',
      fecalqx_srT2D_BeforeAfter_Yr == 'non_cases' ~ 'non_cases',
      TRUE  ~ NA_character_
    ),
    fecalqx_T2D_BeforeAfter_Pool = case_when (
      fecalqx_T2D_BeforeAfter == 'oldT2Dconfirmed' | (metf17+insul17+jard17+invk17+sita17+ohypo17+ninsu17>0) ~ 'oldT2Dconfirmed',
      fecalqx_T2D_BeforeAfter == 'newT2Dconfirmed' | (metf17+insul17+jard17+invk17+sita17+ohypo17+ninsu17==0)&(metf19+insul19+jard19+sita19+farx19+invk19+ohypo19+ninsu19>0) ~ 'newT2Dconfirmed',
      fecalqx_T2D_BeforeAfter == 'case_confirmed_unknown_onset' ~ 'case_confirmed_unknown_onset',
      fecalqx_T2D_BeforeAfter == 'non_cases' ~ 'non_cases',
      TRUE  ~ NA_character_
    ),
    fecalqx_T2D_BeforeAfter_Yr_Pool = case_when (
      fecalqx_T2D_BeforeAfter_Yr == 'oldT2Dconfirmed' | (metf17+insul17+jard17+invk17+sita17+ohypo17+ninsu17>0) ~ 'oldT2Dconfirmed',
      fecalqx_T2D_BeforeAfter_Yr == 'newT2Dconfirmed' | (metf17+insul17+jard17+invk17+sita17+ohypo17+ninsu17==0)&(metf19+insul19+jard19+sita19+farx19+invk19+ohypo19+ninsu19>0) ~ 'newT2Dconfirmed',
      fecalqx_T2D_BeforeAfter_Yr == 'case_confirmed_unknown_onset' ~ 'case_confirmed_unknown_onset',
      fecalqx_T2D_BeforeAfter_Yr == 'non_cases' ~ 'non_cases',
      TRUE  ~ NA_character_
    )
  ) %>% 
  #dplyr::rename(age=age19,
  #             calor=calor19n,
  #            alco=alco19n,
  #           ahei=nAHEI19a,
  #          act=act17m) %>% 
  #mutate(smk = case_when(
  # smk19==1|smk19==2 ~ "never/past",  # Never smoked or past
  #smk19==3 ~ 'current',
  #TRUE ~ NA_character_                                           # All other cases
  #)) %>% 
  dplyr::mutate(bristol_cats=case_when(
    hardst==1 ~ 0,
    lumpst==1 ~ 1,
    crckst==1 ~ 2,
    softst==1 ~ 3,
    blobst==1 ~ 4,
    mushst==1 ~ 5,
    wtryst==1 ~ 6,
    TRUE ~ NA_integer_
  )) %>% 
  dplyr::mutate(medss=case_when(
    col2mn==2 ~ 1, # 2=yes N=10, 1=no N=881, 3=pt, N=3
    TRUE ~ 0
  ))# %>% 
# mutate(wt19=case_when(
# !is.na(wt19) ~ wt19,
# is.na(wt19) ~ (bmi19*(height89*0.0254)^2)*2.20462,
#  TRUE ~ NA_real_
# )) 
dim(MicroN_IncidentT2D_nonExc)
# 6116  362
table(MicroN_IncidentT2D_nonExc$fecalqx_T2D_BeforeAfter)
table(MicroN_IncidentT2D_nonExc$fecalqx_T2D_BeforeAfter_Pool)
table(MicroN_IncidentT2D_nonExc$fecalqx_T2D_BeforeAfter_Yr)
table(MicroN_IncidentT2D_nonExc$fecalqx_T2D_BeforeAfter_Yr_Pool)
table(MicroN_IncidentT2D_nonExc$fecalqx_srT2D_BeforeAfter)#Need
table(MicroN_IncidentT2D_nonExc$fecalqx_srT2D_BeforeAfter_Yr)
table(MicroN_IncidentT2D_nonExc$fecalqx_srT2D_BeforeAfter_Pool)#Need
table(MicroN_IncidentT2D_nonExc$fecalqx_srT2D_BeforeAfter_Yr_Pool)
###############Starting orgnize nut related microbiome data####################
template_df <- data.frame(taxon_new = msl_sub$feature)
ordered_annots <- left_join(template_df, Annot_Taxon, by = "taxon_new")
target_names <- ordered_annots %>% pull(name)
target_names
#####Prescient#########
target_codes_df <- data.frame(name= target_names)
renamed_df_pre <- left_join(target_codes_df, prescient_annot_taxon,by="name")
renamed_species_pre <- renamed_df_pre %>%
  pull(taxon_new)
renamed_species_pre#30
#######B2B#########
target_codes_df <- data.frame(name= target_names)
renamed_df_b2b <- left_join(target_codes_df, B2B_annot_taxon,by="name")
renamed_species_b2b <- renamed_df_b2b %>%
  pull(taxon_new)#30
#######CVD#########
target_codes_df <- data.frame(name= target_names)
renamed_df_cvd <- left_join(target_codes_df, MicroN_CVD_annot_key, by = c("name" = "MicroN_CVD_name"))
renamed_species_cvd <- renamed_df_cvd %>%
  pull(MicroN_CVD_taxon_new)#30
renamed_species_cvd
#######T2D#########
target_codes_df <- data.frame(name= target_names)
renamed_df_t2d <- left_join(target_codes_df, MicroN_T2D_annot_key, by = c("name" = "MicroN_T2D_name"))
renamed_species_t2d <- renamed_df_t2d %>%
  pull(MicroN_T2D_taxon_new)#16 in 19
renamed_species_t2d# 26
#######Diverticulitis_lab20246_20247#########
target_codes_df <- data.frame(name= target_names)
renamed_df_diver_1 <- left_join(target_codes_df, Diverticulitis_lab20246_20247_annot_taxon ,by="name")
renamed_species_diver_1<- renamed_df_diver_1 %>%
  pull(taxon_new)#30
renamed_species_diver_1
#######Diverticulitis_lab20200#########
target_codes_df <- data.frame(name= target_names)
renamed_df_diver_2 <- left_join(target_codes_df,Diverticulitis_lab20200_annot_taxon ,by="name")
renamed_species_diver_2<- renamed_df_diver_2 %>%
  pull(taxon_new)#30
renamed_species_diver_2
#######Kidney Stone#########
target_codes_df <- data.frame(name= target_names)
renamed_df_KidneyStone <- left_join(target_codes_df, KidneyStone_annot_taxon,by="name")
renamed_species_KidneyStone<- renamed_df_KidneyStone %>%
  pull(taxon_new)#30
renamed_species_KidneyStone
#######GDM#########
target_codes_df <- data.frame(name= target_names)
renamed_df_GDM<- left_join(target_codes_df, GDM_annot_taxon,by="name")
renamed_species_GDM<- renamed_df_GDM %>%
  pull(taxon_new)#26
renamed_species_GDM
#######Poly1#########
target_codes_df <- data.frame(name= target_names)
renamed_df_poly1<- left_join(target_codes_df, MicroN_inc_polyp_batch1_annot_key, c("name" ="MicroN_inc_polyp_batch1_name"))
renamed_species_poly1<- renamed_df_poly1 %>%
  pull(MicroN_inc_polyp_batch1_taxon_new)#6 in 8
renamed_species_poly1#30
#######Poly2#########
target_codes_df <- data.frame(name= target_names)
MicroN_inc_polyp_batch2_annot_key$MicroN_inc_polyp_batch2_name
renamed_df_poly2<- left_join(target_codes_df,MicroN_inc_polyp_batch2_annot_key, c("name" ="MicroN_inc_polyp_batch2_name"))
renamed_species_poly2<- renamed_df_poly2%>%
  pull(MicroN_inc_polyp_batch2_taxon_new)#6 in 8
renamed_species_poly2#26
#######Cognitive#########
target_codes_df <- data.frame(name= target_names)
renamed_df_cognitive<- left_join(target_codes_df,MicroN_Cognitive_annot_key,c("name" ="MicroN_Cognitive_name"))
renamed_species_cognitive<- renamed_df_cognitive %>%
  pull(MicroN_Cognitive_taxon_new)
renamed_species_cognitive#29
########################rename####################
MicroN_Prescient_Microb_s_avg$id<-rownames(MicroN_Prescient_Microb_s_avg)
MicroN_B2B_Microb_s_avg$id<-rownames(MicroN_B2B_Microb_s_avg)
MicroN_T2D_Microb_s_avg$id<-rownames(MicroN_T2D_Microb_s_avg)
MicroN_CVD_Microb_s_avg$id<-rownames(MicroN_CVD_Microb_s_avg)
MicroN_Diverticulitis_lab20200_Microb_s_avg$id<-rownames(MicroN_Diverticulitis_lab20200_Microb_s_avg)
MicroN_KidneyStone_Microb_s_avg$id<-rownames(MicroN_KidneyStone_Microb_s_avg)
MicroN_GDM_Microb_s_avg$id<-rownames(MicroN_GDM_Microb_s_avg)
MicroN_Diverticulitis_lab20246_20247_Microb_s_avg$id<-rownames(MicroN_Diverticulitis_lab20246_20247_Microb_s_avg)
MicroN_inc_polyp_batch1_Microb_s_avg$id<-rownames(MicroN_inc_polyp_batch1_Microb_s_avg)
MicroN_inc_polyp_batch2_Microb_s_avg$id<-rownames(MicroN_inc_polyp_batch2_Microb_s_avg)
MicroN_Cognitive_Microb_s_avg$id<-rownames(MicroN_Cognitive_Microb_s_avg)
#######7 sub data###########
library(dplyr)

#===========================================================
# 1. 定义一个通用的处理函数
#===========================================================
process_microbiome_data <- function(raw_data, map_df, old_taxon_col, target_std_names) {
  
  # A. 确保id列存在
  if(!"id" %in% colnames(raw_data)) {
    raw_data$id <- rownames(raw_data)
  }
  
  # B. 建立映射向量 (从 旧代码 -> msl_sub$feature标准名)
  # 注意：这里基于前提 map_df 的行顺序是对应 target_std_names (msl_sub$feature) 的顺序
  # (因为你是用 ordered target_names join 得到的 map_df)
  
  # 获取当前数据集对应的旧代码列
  current_codes <- map_df[[old_taxon_col]]
  
  # 创建映射字典: names=旧代码, values=新代码(msl_sub$feature)
  # 这一步会自动把 NA (缺失的菌) 包含在内
  code_map <- setNames(target_std_names, current_codes)
  
  # C. 移除旧代码为 NA 的项 (即该数据集中不存在这个菌)
  valid_map <- code_map[!is.na(names(code_map))]
  
  # D. 从原始数据中提取存在的列
  # 只选择 id 和 那些在这个数据集中实际存在的菌代码
  cols_to_select <- c("id", names(valid_map))
  
  # 确保只选数据里真有的列 (防止map里有但数据里没有的异常情况)
  cols_to_select <- intersect(cols_to_select, colnames(raw_data))
  
  data_sub <- raw_data %>%
    dplyr::select(all_of(cols_to_select))
  
  # E. 重命名 (把旧代码 S_xxxx 变成 msl_sub$feature 里的名字)
  # rename_with 会根据 valid_map 自动查找并替换
  data_renamed <- data_sub %>%
    rename_with(~ valid_map[.], .cols = any_of(names(valid_map)))
  
  # F. 补齐缺失的列并设为 0
  # 找出 msl_sub$feature 中有哪些名字目前还没在 data_renamed 里
  missing_cols <- setdiff(target_std_names, colnames(data_renamed))
  
  if(length(missing_cols) > 0) {
    # 批量赋值为 0
    data_renamed[missing_cols] <- 0
  }
  
  # G. 最后整理顺序: id 在前，后面严格按照 msl_sub$feature 排序
  data_final <- data_renamed %>%
    dplyr::select(id, all_of(target_std_names))
  
  return(data_final)
}

#===========================================================
# 2. 批量应用函数到 7 个数据集
#===========================================================
# 1. Prescient (19/19)
prescient_taxon_total_sub <- process_microbiome_data(
  raw_data = MicroN_Prescient_Microb_s_avg,
  map_df = renamed_df_pre,
  old_taxon_col = "taxon_new",
  target_std_names = msl_sub$feature
)
prescient_taxon_total_sub
# 2. B2B (19/19)
B2B_taxon_total_sub <- process_microbiome_data(
  raw_data = MicroN_B2B_Microb_s_avg,
  map_df = renamed_df_b2b,
  old_taxon_col = "taxon_new",
  target_std_names = msl_sub$feature
)
# 3. CVD (19/19)
CVD_taxon_total_sub <- process_microbiome_data(
  raw_data = MicroN_CVD_Microb_s_avg,
  map_df = renamed_df_cvd,
  old_taxon_col = "MicroN_CVD_taxon_new",
  target_std_names = msl_sub$feature
)

# 4. T2D (16/19 - 有缺失值)
# 函数会自动检测 renamed_df_t2d 中的 NA，并在结果中将对应的 msl_sub$feature 列设为 0
T2D_taxon_total_sub <- process_microbiome_data(
  raw_data = MicroN_T2D_Microb_s_avg,
  map_df = renamed_df_t2d,
  old_taxon_col = "MicroN_T2D_taxon_new",
  target_std_names = msl_sub$feature
)

# 5. Diverticulitis (19/19)
Diverticulitis_lab20200_taxon_total_sub <- process_microbiome_data(
  raw_data = MicroN_Diverticulitis_lab20200_Microb_s_avg,
  map_df = renamed_df_diver_2,
  old_taxon_col = "taxon_new",
  target_std_names = msl_sub$feature
)
Diverticulitis_lab20200_taxon_total_sub
# 6. Kidney Stone (19/19)
KidneyStone_taxon_total_sub <- process_microbiome_data(
  raw_data = MicroN_KidneyStone_Microb_s_avg,
  map_df = renamed_df_KidneyStone,
  old_taxon_col = "taxon_new",
  target_std_names = msl_sub$feature
)
KidneyStone_taxon_total_sub
# 7. GDM (16/19 - 有缺失值)
GDM_taxon_total_sub <- process_microbiome_data(
  raw_data = MicroN_GDM_Microb_s_avg,
  map_df = renamed_df_GDM,
  old_taxon_col = "taxon_new",
  target_std_names = msl_sub$feature
)
GDM_taxon_total_sub
# 8. Diverticulitis_2
Diverticulitis_lab20246_20247_taxon_total_sub <- process_microbiome_data(
  raw_data = MicroN_Diverticulitis_lab20246_20247_Microb_s_avg,
  map_df = renamed_df_diver_1,
  old_taxon_col = "taxon_new",
  target_std_names = msl_sub$feature
)
Diverticulitis_lab20246_20247_taxon_total_sub

# 9. poly1
MicroN_inc_polyp_batch1_taxon_total_sub <- process_microbiome_data(
  raw_data = MicroN_inc_polyp_batch1_Microb_s_avg,
  map_df = renamed_df_poly1,
  old_taxon_col = "MicroN_inc_polyp_batch1_taxon_new",
  target_std_names = msl_sub$feature
)
MicroN_inc_polyp_batch1_taxon_total_sub
# 10. poly2
MicroN_inc_polyp_batch2_taxon_total_sub <- process_microbiome_data(
  raw_data = MicroN_inc_polyp_batch2_Microb_s_avg,
  map_df = renamed_df_poly2,
  old_taxon_col = "MicroN_inc_polyp_batch2_taxon_new",
  target_std_names = msl_sub$feature
)
MicroN_inc_polyp_batch2_taxon_total_sub
# 11. Cognitive function
MicroN_Cognitive_taxon_total_sub <- process_microbiome_data(
  raw_data =MicroN_Cognitive_Microb_s_avg,
  map_df = renamed_df_cognitive,
  old_taxon_col = "MicroN_Cognitive_taxon_new",
  target_std_names = msl_sub$feature
)
MicroN_Cognitive_taxon_total_sub

#===========================================================
# 3. 检查结果
#===========================================================

# 检查列名是否统一
print(names(T2D_taxon_total_sub))
print(names(prescient_taxon_total_sub))

# 检查 T2D 中自动补0的列是否生效 (比如检查几个应该为0的列)
# 假设 msl_sub$feature 里有但 T2D 原始数据里没有的，现在应该全是 0
# 我们可以查看 T2D 数据框的列数，应该等于 1(id) + 19(msl_sub$feature) = 20
print(dim(T2D_taxon_total_sub))
print(dim(GDM_taxon_total_sub))

#######1.
# 1. 创建一个包含所有数据框的命名列表
# 列表的名字 (e.g., "Prescient", "B2B") 将被用作 cohort 的值
df_list <- list(
  Prescient = prescient_taxon_total_sub,
  B2B = B2B_taxon_total_sub,
  CVD = CVD_taxon_total_sub,
  Diverticulitis_lab20200 = Diverticulitis_lab20200_taxon_total_sub,
  KidneyStone = KidneyStone_taxon_total_sub,
  T2D = T2D_taxon_total_sub,
  GDM = GDM_taxon_total_sub,
  poly1=MicroN_inc_polyp_batch1_taxon_total_sub,
  poly2=MicroN_inc_polyp_batch2_taxon_total_sub,
  Diverticulitis_lab20246_20247=Diverticulitis_lab20246_20247_taxon_total_sub,
  Cognitive=MicroN_Cognitive_taxon_total_sub
)

# 1. 准备阶段：将列表合并为长表，并暂存来源信息
combined_raw <- imap_dfr(df_list, ~{
  df <- .x
  # 强制将第一列重命名为 SampleID，确保合并时列名一致
  colnames(df)[1] <- "SampleID" 
  df$cohort_source <- .y # 暂存原始队列名
  return(df)
})
#######1.
# 1. 创建一个包含所有数据框的命名列表
# 列表的名字 (e.g., "Prescient", "B2B") 将被用作 cohort 的值
df_list <- list(
  Prescient = prescient_taxon_total_sub,
  B2B = B2B_taxon_total_sub,
  CVD = CVD_taxon_total_sub,
  Diverticulitis_lab20200 = Diverticulitis_lab20200_taxon_total_sub,
  KidneyStone = KidneyStone_taxon_total_sub,
  T2D = T2D_taxon_total_sub,
  GDM = GDM_taxon_total_sub,
  poly1=MicroN_inc_polyp_batch1_taxon_total_sub,
  poly2=MicroN_inc_polyp_batch2_taxon_total_sub,
  Diverticulitis_lab20246_20247=Diverticulitis_lab20246_20247_taxon_total_sub,
  Cognitive=MicroN_Cognitive_taxon_total_sub
)
# 1. 准备阶段：将列表合并为长表，并暂存来源信息
combined_raw <- imap_dfr(df_list, ~{
  df <- .x
  # 强制将第一列重命名为 SampleID，确保合并时列名一致
  colnames(df)[1] <- "SampleID" 
  df$cohort_source <- .y # 暂存原始队列名
  return(df)
})

# 2. 定义核心处理函数
process_global_mix <- function(df) {
  
  # --- 步骤 A: 留痕统计 ---
  total_rows <- nrow(df)
  unique_ids <- n_distinct(df$SampleID)
  
  cat("======= 数据合并与去重日志 =======\n")
  cat(sprintf("1. 原始总观测行数: %d\n", total_rows))
  cat(sprintf("2. 唯一 ID 数量:     %d\n", unique_ids))
  cat(sprintf("3. 涉及重复的行数:   %d\n", total_rows - unique_ids))
  
  # --- 步骤 B: 分组、平均与标记 Mix ---
  df_avg <- df %>%
    group_by(SampleID) %>%
    summarise(
      # 1. 数值列（菌群）：取平均
      across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      
      # 2. 来源列：判断逻辑
      # n_distinct(cohort_source) > 1 说明该ID出现在至少两个不同的队列中
      cohort = if(n_distinct(cohort_source) > 1) "Mix" else first(cohort_source)
    ) %>%
    ungroup()
  
  # --- 步骤 C: 打印 Mix 的统计信息 ---
  n_mix <- sum(df_avg$cohort == "Mix")
  cat(sprintf("4. 最终结果中标记为 'Mix' 的 ID 数量: %d\n", n_mix))
  cat("==================================\n")
  
  # --- 步骤 D: 数据转换 (Arcsin-Sqrt) ---
  # 警告：summarise 后，SampleID 是第1列，cohort 变成了最后一列。
  # 假设 SampleID 后紧接着就是那 8 个数值型的菌群列 (即 2:7)
  # 建议：运行 head(df_avg) 确认一下列的顺序是否符合预期
  
  df_avg[2:8] <- df_avg[2:8] / 100
  df_avg[2:8] <- asin(sqrt(df_avg[2:8]))
  
  #df_avg$MRS<- predict(final_model, newdata = df_avg[msl_sub$feature])
  
  return(df_avg)
}

# 3. 执行
final_combined_df_tidy <- process_global_mix(combined_raw)

#======= 数据合并与去重日志 =======
#  1. 原始总观测行数: 6644
#2. 唯一 ID 数量:     6116
#3. 涉及重复的行数:   528
#4. 最终结果中标记为 'Mix' 的 ID 数量: 507
#==================================
# 4. 检查结果
head(final_combined_df_tidy)
sig_dta <- data.frame(
  feature = msl_sub$feature,
  coef = msl_sub$coef
)
final_combined_df_tidy<- rank_score(final_combined_df_tidy, sig_dta, varname = "coef", feature_col = "feature")
mn_final1<-merge(MicroN_IncidentT2D_nonExc,final_combined_df_tidy,all.x = T,by.x="id",by.y = "SampleID")
mn_final1$MRS<-mn_final1$rank.score
summary(mn_final1$MRS)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.19895 -0.04768 -0.01239 -0.01511  0.02245  0.11585 

table(mn_final1$cohort)
stool_cov<-read.csv("/proj/n2dats/n2dat02/subStudy/MicroN/stool_collection_qx/nhs2_micron.042123.csv")
stool_cov2<-read.csv("/proj/n2dats/n2dat02/subStudy/MicroN/stool_collection_qx/nhs2_micron.csv")
stool_cov <- stool_cov %>%
  mutate(
    # 直接使用列名，无需 `$`
    stool_type = case_when(
      hardst == 1  ~ 1,
      lumpst == 1  ~ 2,
      crckst == 1  ~ 3,
      softst == 1  ~ 4,
      blobst == 1  ~ 5,
      mushst == 1  ~ 6,
      wtryst == 1  ~ 7,
      TRUE         ~ NA_real_ # 保持默认值为 NA
    )
  ) %>%
  mutate(
    # 根据 stool_type 创建新的分类变量 bs
    bs = case_when(
      stool_type %in% c(1, 2)  ~ 1,
      stool_type %in% c(3, 4)  ~ 2,
      stool_type %in% c(5, 6, 7) ~ 3,
      TRUE                     ~ 2 #众数
    )
  ) %>%
  mutate(
    # 创建哑变量，这里假设当 bs 为 NA 时，结果也为 NA，如果想为 0 请参考上面的建议
    bs_1 = ifelse(bs == 1, 1, 0),
    bs_2 = ifelse(bs == 2, 1, 0),
    bs_3 = ifelse(bs == 3, 1, 0)
  ) %>%
  mutate(
    # 将多个变量的转换合并在一个 mutate 调用中，更高效
    orab_c   = ifelse(is.na(orabnu) | orabnu != 1, 1, 0),
    injab    = ifelse(is.na(injabnu) | injabnu != 1, 1, 0),
    
    # 【已修正】添加了比较运算符 ==
    abuse    = ifelse(orab_c == 1 | injab == 1, 1, 0),
    
    preb_c   = ifelse(is.na(preb) | preb != 1, 0, 1), # 更简洁的写法，如果 preb=1 为 1, 其他(包括NA)为 0
    prob_c   = ifelse(is.na(prob) | prob != 1, 0, 1),
    col2mn_c = ifelse(is.na(col2mn) | col2mn == 1, 0, 1),
    
    # 转换年份
    yr2      = case_when(
      yr == 2 ~ 2019,
      yr == 3 ~ 2020,
      yr == 4 ~ 2021,
      yr == 5 ~ 2022,
      TRUE    ~ NA_real_
    )
  )
nhs2_cov<-read.csv("/udd/n2gji/micro/data/micron_cov.csv")
cov_final<-merge(stool_cov,nhs2_cov,by="id",all.x = T) #1741 169

cov_final <- cov_final%>%
  mutate(
    nuts = coalesce(nuts, 0), 
    nutsoft = if_else(nutsoft == 9, 1, nutsoft))%>%
  mutate(smkc2_new = if_else(smkc2 == 0, 0, 1),
         alcocat = coalesce(alcocat, 1))
source("/udd/hpyah/project/microb_metabo/MLVSfunctions.R")
cov_final$datem_qu = date2proc("yr2","mo","dd",cov_final,"months")#MLVSfunction
cov_final$age = (cov_final$datem_qu-cov_final$birthday)/12
cov_final$age = ifelse(is.na(cov_final$age),median(cov_final$age,na.rm = T),cov_final$age)
cov_final <- cov_final %>%
  mutate(smkc2_new = ifelse(smkc2 == 0, 0, 1)) %>%
  mutate(avg_calor= ifelse(is.na(avg_calor), median(avg_calor,na.rm=T), avg_calor)) %>%
  mutate(avg_ahei= ifelse(is.na(avg_ahei), median(avg_ahei,na.rm=T), avg_ahei)) %>%
  mutate(smkc2_new=ifelse(is.na(smkc2_new),0,smkc2_new))%>%
  mutate(bmi=case_when(
    is.na(bmi) ~ median(bmi,na.rm=TRUE), # impute 88 missing with median
    TRUE ~ bmi
  )) %>%
  mutate(bmiat=case_when(
    bmi<25 ~ 1,
    25<=bmi & bmi<30 ~ 2,
    bmi>=30 ~ 3))
mn_final<-merge(mn_final1,cov_final[,c("id","alco","alcocat","bmi","calor","avg_calor", "ahei","avg_ahei",
                                       "pa","avg_pa","smk","smkc","smkc2","smkc2_new","age","bmiat",
                                       "prob_c", "abuse", "col2mn_c",
                                       "bs_2", "bs_3","bs")],by="id",all.x = T)

mn_final<-mn_final%>% mutate(
  avg_calorq = ntile(avg_calor, 5),
  avg_paq    = ntile(avg_pa, 5),
  avg_aheiq=ntile(avg_ahei, 5)
)
############Define outcomes####################
table(MicroN_IncidentT2D_nonExc$fecalqx_srT2D_BeforeAfter)#Need
table(MicroN_IncidentT2D_nonExc$fecalqx_srT2D_BeforeAfter_Pool)#Need
mn_final <- mn_final %>%
  mutate(
    t2d_confirmed = case_when(
      fecalqx_T2D_BeforeAfter == "newT2Dconfirmed" ~ 1,
      TRUE ~ 0
    ),
    t2d_med_confirmed = case_when(
      fecalqx_T2D_BeforeAfter_Pool == "newT2Dconfirmed" ~ 1,
      TRUE ~ 0
    ),
    
    t2d_confirmed_p = case_when(
      fecalqx_T2D_BeforeAfter != "non_cases" ~ 1,
      TRUE ~ 0
    ),
    
    t2d_med_confirmed_p = case_when(
      fecalqx_T2D_BeforeAfter_Pool != "non_cases" ~ 1,
      TRUE ~ 0
    )
  )
table(mn_final$t2d_confirmed)#108
table(mn_final$t2d_med_confirmed)#157
table(mn_final$t2d_confirmed_p)#696
table(mn_final$t2d_med_confirmed_p)#884
mn_final<-subset(mn_final,!duplicated(mn_final$id))
mn_final_exdb<-mn_final %>%
  filter(mn_final$fecalqx_T2D_BeforeAfter != "oldT2Dconfirmed")
mn_final_exdbmed<-mn_final %>%
  filter(mn_final$fecalqx_T2D_BeforeAfter_Pool != "oldT2Dconfirmed")

#############Starting analysis######################
#1. Prevalent T2D#
mn_final_exdbmed <- mn_final_exdbmed %>%
  mutate(
    MRSt = factor(ntile(MRS, 3)))

mn_final_exdb <- mn_final_exdb %>%
  mutate(
    MRSt = factor(ntile(MRS, 3)))

mn_final <- mn_final %>%
  mutate(
    MRSt = factor(ntile(MRS, 3)))
#########Prevalence Db#################
mn_final$MRSs<-scale(mn_final$MRS)
mn_final[msl_sub$feature]<-scale(mn_final[msl_sub$feature])
mn_final <- mn_final %>%
  mutate(
    MRSq = factor(ntile(MRS, 4)))
table(mn_final$bmiat)
x_vars <- c("MRSs", "MRSq", msl_sub$feature)

y_vars <- c("t2d_confirmed_p", "t2d_med_confirmed_p")
mn_final$alcocat
# 定义协变量
cov1 <- c("age", "factor(cohort)")
cov2 <- c("age", "avg_pa", "avg_calor", "factor(alcocat)",
          "prob_c", "abuse", "col2mn_c", "bs_2", "bs_3", "factor(cohort)")
cov3 <- c("factor(bmiat)", "age", "avg_pa", "avg_calor", "factor(alcocat)", 
          "prob_c", "abuse", "col2mn_c", "bs_2", "bs_3", "factor(cohort)")

# 2. 模型定义
model_specs <- list(
  "Model_1_Basic" = cov1,
  "Model_2_Lifestyle" = cov2,
  "Model_3_Full" = cov3
)

final_results_table <- data.frame()

cat("Starting to run all models...\n")

for (y in y_vars) {
  for (x in x_vars) {
    for (model_name in names(model_specs)) {
      
      
      current_covs <- model_specs[[model_name]]
      
      
      formula_string <- paste(y, "~", x, "+", paste(current_covs, collapse = " + "))
      model_formula <- as.formula(formula_string)
      
      
      model <- glm(model_formula, data = mn_final, family = binomial())
      
      
      summary_coeffs <- summary(model)$coefficients
      matched_rows <- grep(paste0("^", x), rownames(summary_coeffs), value = TRUE)
      
      
      if (length(matched_rows) == 0 && x %in% rownames(summary_coeffs)) {
        matched_rows <- x
      }
      
      # 对找到的每一行（可能是1行，也可能是3行）进行提取
      for (term_name in matched_rows) {
        
        result_row <- summary_coeffs[term_name, ]
        
        current_result <- data.frame(
          model_name = model_name,
          outcome = y,
          predictor = x,          # 原始变量名 (例如 MRSq)
          term = term_name,       # 具体项名 (例如 MRSq2, MRSq3)
          estimate = result_row["Estimate"],
          std.error = result_row["Std. Error"],
          p.value = result_row["Pr(>|z|)"],
          nobs = nobs(model),
          OR = exp(result_row["Estimate"]), # 顺便算个OR
          CI_low = exp(result_row["Estimate"] - 1.96 * result_row["Std. Error"]),
          CI_high = exp(result_row["Estimate"] + 1.96 * result_row["Std. Error"])
        )
        
        final_results_table <- rbind(final_results_table, current_result)
      }
    }
  }
}

cat("All models finished!\n\n")
final_results_table <- final_results_table %>%
  select(outcome, predictor, term, model_name, OR, CI_low, CI_high, p.value, everything()) %>%
  arrange(outcome, predictor, model_name)
print(head(final_results_table))
final_results_table_p <- final_results_table
#########Cox regression##########
table(mn_final$dtdxdb2_clean)
mn_final_exdb$alcoft<-ifelse(is.na(mn_final_exdb$alcoft),median(mn_final_exdb$alcoft,na.rm = T),mn_final_exdb$alcoft)
mn_final$follow_data <- pmin(mn_final$dtdth, 1494, mn_final$dtdxdb2_clean, na.rm = TRUE)

mn_final_exdb$follow_data <- pmin(mn_final_exdb$dtdth, 1494, mn_final_exdb$dtdxdb2_clean, na.rm = TRUE)
mn_final_exdbmed$follow_data <- pmin(mn_final_exdbmed$dtdth, 1494, mn_final_exdbmed$dtdxdb2_clean, na.rm = TRUE)
summary(mn_final_exdb$fecalmonths)
mn_final_exdb$fecalmonths<-ifelse(is.na(mn_final_exdb$fecalmonths),median(mn_final_exdb$fecalmonths,na.rm = T),mn_final_exdb$fecalmonths)
mn_final_exdb$tdb2<-mn_final_exdb$follow_data-mn_final_exdb$fecalmonths
mn_final_exdb$MRSs<-scale(mn_final_exdb$MRS)
mn_final_exdb[msl_sub$feature]<-scale(mn_final_exdb[msl_sub$feature])
mn_final_exdb$avg_pa
mn_final_exdb <- mn_final_exdb %>%
  mutate(
    MRSq = factor(ntile(MRS, 4)))
x_vars <- c("MRSs", "MRSq", msl_sub$feature)
y_vars <- c("t2d_confirmed_p","t2d_med_confirmed_p")
summary(coxph(Surv(tdb2,t2d_confirmed) ~ MRSt +age+ factor(avg_paq)+ 
                factor(avg_calorq) + alcoft+ 
                prob_c + abuse + col2mn_c + bs_2 + bs_3 + strata(cohort)+factor(mn_final_exdb$bmiat), data = mn_final_exdb))



cov1 <- c("age", "strata(cohort)")
cov2 <- c("age", "factor(avg_paq)", "factor(avg_calorq)",  "factor(avg_aheiq)",
          "prob_c", "abuse", "col2mn_c", "bs_2", "bs_3", "strata(cohort)")
cov3 <- c("factor(bmiat)", "age",  "factor(avg_paq)", "factor(avg_calorq)",  "factor(avg_aheiq)",
          "prob_c", "abuse", "col2mn_c", "bs_2", "bs_3", "strata(cohort)")

model_specs <- list(
  "Model_1_Basic" = cov1,
  "Model_2_Lifestyle" = cov2,
  "Model_3_Full" = cov3
)
paste(cov3, collapse = " + ")
# ==========================================
# 2. 循环运行 Cox 模型
# ==========================================
results_list <- list()
summary(mn_final_exdb$tdb2)
cat("Starting Cox models...\n")
for (x in x_vars) {
  for (model_name in names(model_specs)) {
    
    covariates <- model_specs[[model_name]]
    
    formula_str <- paste("Surv(tdb2, t2d_confirmed) ~ `", x, "` + ", 
                         paste(covariates, collapse = " + "), sep = "")
    
    tryCatch({
      fit <- coxph(as.formula(formula_str), data = mn_final_exdb)
      fit_sum <- summary(fit)
      
      clean_x <- gsub("`", "", x)
      
      all_coef_names <- rownames(fit_sum$coefficients)
      
      
      target_rows <- grep(paste0("^", clean_x), all_coef_names, value = TRUE)
      
      if (length(target_rows) == 0) {
        target_rows <- grep(clean_x, all_coef_names, fixed = TRUE, value = TRUE)
      }
      
      for (term_name in target_rows) {
        
        coef_row <- fit_sum$coefficients[term_name, , drop = FALSE]
        conf_row <- fit_sum$conf.int[term_name, , drop = FALSE]
        
        res_row <- data.frame(
          Variable = x,          
          Term = term_name,      
          Model = model_name,
          HR = conf_row[1],        # exp(coef)
          CI_Lower = conf_row[3],  # lower .95
          CI_Upper = conf_row[4],  # upper .95
          P_Value = coef_row[5],   # Pr(>|z|)
          N = fit$n,
          Events = fit$nevent
        )
        list_key <- paste(term_name, model_name, sep = "_")
        results_list[[list_key]] <- res_row
      }
      
      # ------------------------------------------------------
      # 修改结束
      # ------------------------------------------------------
      
    }, error = function(e) {
      message(paste("Error processing:", x, "in", model_name, "-", e$message))
    })
  }
}

cat("Cox models finished!\n")

final_results_cox <- do.call(rbind, results_list)
rownames(final_results_cox) <- NULL # 重置行名


final_results_cox <- final_results_cox %>%
  arrange(Variable, Model, Term)


final_results_cox$HR <- round(final_results_cox$HR, 3)
final_results_cox$CI_Lower <- round(final_results_cox$CI_Lower, 3)
final_results_cox$CI_Upper <- round(final_results_cox$CI_Upper, 3)

final_results_cox$P_Value_Fmt <- ifelse(final_results_cox$P_Value < 0.001, "<0.001", 
                                        round(final_results_cox$P_Value, 3))


print(head(final_results_cox))
setwd("/udd/n2gji/micro/dietpt/final/data")
write.csv(final_results_cox, "cox_results_micro.csv", row.names = FALSE)
write.csv(final_results_table_p,file = "mn_prevalence.csv")
########################Heatmap#####################
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(RColorBrewer)

# ==============================================================================
# Part 1: Microbiome and score
# ==============================================================================
mapping_df <- enframe(name_mapping, name = "Original_Name", value = "Current_Name")

coef_df$new_name <- coef_df$Variable
coef_df$new_name <- sub("_log_z", "", coef_df$new_name)
coef_df <- coef_df %>%
  left_join(mapping_df, by = c("new_name" = "Current_Name"))

micro_list <- as.character(msl_sub$feature)
food_list <- as.character(coef_df$Original_Name)

total_rows_food <- length(food_list) * length(micro_list)
results_food <- data.frame(
  Food = character(total_rows_food),
  Metabolite = character(total_rows_food), # 实际上是微生物，为了命名一致保留
  Beta = numeric(total_rows_food),
  SE = numeric(total_rows_food),
  P_value = numeric(total_rows_food),
  N_obs = numeric(total_rows_food),
  stringsAsFactors = FALSE
)
ToBeUsed.Taxon <- as.data.frame(ToBeUsed.Taxon)
covariates <- "+ ageyr + totMETs_paq + calor_fo_dr_wtavg+alco_g + smoke_bld + bmi_bld + probio_2m_fec + antibio_12m_fec + colsc_2m_fec + acid_2m_fec + stooltype_fec.1 + stooltype_fec.2 + stooltype_fec.3 + stooltype_fec.4 + stooltype_fec.5 + stooltype_fec.6+alco_bld"

row_idx <- 1
for (f_item in food_list) {
  ToBeUsed.Taxon$current_food_s <- as.numeric(scale(ToBeUsed.Taxon[[f_item]]))
  for (m_item in micro_list) {
    ToBeUsed.Taxon$current_micro_s <- as.numeric(scale(suppressWarnings(asin(sqrt(ToBeUsed.Taxon[[m_item]])))))
    formula_str <- paste0("current_micro_s ~ current_food_s ", covariates)
    fit <- tryCatch({ lm(as.formula(formula_str), data = ToBeUsed.Taxon) }, error = function(e) NULL)
    
    if(!is.null(fit)) {
      summ <- summary(fit)
      results_food[row_idx, "Food"]       <- f_item
      results_food[row_idx, "Metabolite"] <- m_item
      results_food[row_idx, "Beta"]       <- round(summ$coefficients[2, 1], 4)
      results_food[row_idx, "SE"]         <- round(summ$coefficients[2, 2], 4)
      results_food[row_idx, "P_value"]    <- summ$coefficients[2, 4]
      results_food[row_idx, "N_obs"]      <- nrow(model.frame(fit))
      row_idx <- row_idx + 1
    }
  }
}


results_food <- results_food[results_food$Food != "", ]
results_food$P_FDR <- p.adjust(results_food$P_value, method = "fdr")

results_food <- merge(results_food, msl_sub[,c("feature","new_name")], by.x = 'Metabolite', by.y="feature", all = FALSE)
results_food
# ==============================================================================
# Part 2: Microbiome and biomarkers
# ==============================================================================
bio_list <- c("hba1c", "tg", "tc", "hdlc", "crp", "tg.hdlc")
covariates_list <- c("ageyr", "totMETs_paq", "calor_fo_dr_wtavg", "smoke_bld","bmi_bld","alc_bld",
                     "probio_2m_fec", "antibio_12m_fec", "colsc_2m_fec", "acid_2m_fec",
                     "stooltype_fec.1", "stooltype_fec.2", "stooltype_fec.3",
                     "stooltype_fec.4", "stooltype_fec.5", "stooltype_fec.6")

total_rows_bio <- length(bio_list) * length(micro_list)
results_bio <- data.frame(
  Biomarker = character(total_rows_bio),
  Metabolite = character(total_rows_bio),
  Beta = numeric(total_rows_bio),
  SE = numeric(total_rows_bio),
  P_value = numeric(total_rows_bio),
  N_obs = numeric(total_rows_bio),
  stringsAsFactors = FALSE
)

row_idx <- 1
for (b_item in bio_list) {
  if (!b_item %in% colnames(ToBeUsed.Taxon)) next 
  
  raw_bio_data <- ToBeUsed.Taxon[[b_item]]
  valid_indices <- !is.na(raw_bio_data)
  if (sum(valid_indices) < 30) next 
  
  valid_covars_list <- c()
  for (cov in covariates_list) {
    if (cov %in% colnames(ToBeUsed.Taxon) && length(unique(na.omit(ToBeUsed.Taxon[[cov]][valid_indices]))) >= 2) {
      valid_covars_list <- c(valid_covars_list, cov)
    }
  }
  final_covariates_str <- ifelse(length(valid_covars_list) > 0, paste0(" + ", paste(valid_covars_list, collapse = " + ")), "")
  
  ToBeUsed.Taxon$current_bio_s <- as.numeric(scale(raw_bio_data))
  
  for (m_item in micro_list) {
    if (!m_item %in% colnames(ToBeUsed.Taxon)) next
    ToBeUsed.Taxon$current_micro_s <- as.numeric(scale(suppressWarnings(asin(sqrt(ToBeUsed.Taxon[[m_item]])))))
    formula_str <- paste0("current_micro_s ~ current_bio_s", final_covariates_str)
    
    fit <- tryCatch({ lm(as.formula(formula_str), data = ToBeUsed.Taxon, na.action = na.omit) }, error = function(e) NULL)
    
    if (!is.null(fit) && "current_bio_s" %in% rownames(summary(fit)$coefficients)) {
      coefs <- summary(fit)$coefficients
      results_bio[row_idx, "Beta"]    <- round(coefs["current_bio_s", "Estimate"], 4)
      results_bio[row_idx, "SE"]      <- round(coefs["current_bio_s", "Std. Error"], 4)
      results_bio[row_idx, "P_value"] <- coefs["current_bio_s", "Pr(>|t|)"]
      results_bio[row_idx, "N_obs"]   <- nrow(model.frame(fit))
      results_bio[row_idx, "Biomarker"]  <- b_item
      results_bio[row_idx, "Metabolite"] <- m_item
      row_idx <- row_idx + 1
    }
  }
}

# 修正点：必须先去除预留的空行，再算FDR，否则FDR计算会被大量的0或NA带偏
results_bio <- results_bio[results_bio$Biomarker != "", ]
results_bio$P_FDR <- p.adjust(results_bio$P_value, method = "fdr")
results_bio <- merge(results_bio, msl_sub[,c("feature","new_name")], by.x = 'Metabolite', by.y="feature", all = FALSE)


# ==============================================================================
# Part 3 & 4: Data Preparation for Heatmap and Forest Plot
# ==============================================================================
# 重命名表
food_group_labels <- c("beer07_serv"="Beers", "cokh07"="Home-made cookies", "cokr07"="Ready-made cookies","sseaf07"="Shell seafood", "wine07"="Wine", "crsoup07"="Cream Soups", "yogurt07"="Yogurt", "fries07"="French Fries", "liver07"="Liver", "sugbev07"="Sugar-Sweetened Beverages", "lowbev07"="Low-Calorie Beverages", "tea07"="Tea", "rais07"="Raisins and Prunes", "marg07"="Margarine", "othfr07"="Other Fruits", "coff07"="Coffee", "procm07"="Processed Meats", "fruj07"="Fruit Juice", "eggs07"="Eggs", "othveg07"="Other Vegetables", "appear07"="Apples and Pears", "dess07"="Desserts", "rmeat07"="Red Meat", "h2o07"="Water", "cruveg07"="Cruciferous Vegetables", "wgrain07"="Whole Grains", "saldre07"="Salad Dressings", "nightveg07"="Nightshade and Cucurbitaceae Vegetables", "lfveg07"="Leafy Vegetables", "stachveg07"="Starchy Vegetables", "rgrain07"="Refined Grains")

biomarker_labels <- c("hba1c"="HbA1c", "tg.hdlc"="TG/HDL-C", "hdlc"="HDL-C", "tg"="Triglycerides", "tc"="Total cholesterol", "crp"="CRP")

# --- 1. 确定全局微生物顺序 (基于 MRS Score Beta) ---
mat1_beta_raw <- msl_sub %>%
  dplyr::select(new_name, coef) %>%
  distinct(new_name, .keep_all = TRUE) %>% 
  column_to_rownames("new_name") %>% as.matrix()
mat1_beta_raw
met_order <- rownames(mat1_beta_raw)[order(mat1_beta_raw[, 1], decreasing = TRUE)]

mat1_beta <- mat1_beta_raw[met_order, , drop = FALSE]
colnames(mat1_beta) <- "MRS"

# --- 2. 准备 Food 矩阵 ---
results_food_renamed <- results_food %>%
  mutate(Food = dplyr::recode(Food, !!!food_group_labels)) %>%
  filter(Food %in% food_group_labels) 

mat2_beta <- results_food_renamed %>% dplyr::select(new_name, Food, Beta) %>% pivot_wider(names_from = Food, values_from = Beta) %>% column_to_rownames("new_name") %>% as.matrix()
mat2_p    <- results_food_renamed %>% dplyr::select(new_name, Food, P_value) %>% pivot_wider(names_from = Food, values_from = P_value) %>% column_to_rownames("new_name") %>% as.matrix()

valid_food_cols <- intersect(food_group_labels[coef_df$Original_Name[order(coef_df$Coefficient, decreasing = TRUE)]], colnames(mat2_beta))
mat2_beta <- mat2_beta[met_order, valid_food_cols, drop=FALSE]
mat2_p    <- mat2_p[met_order, valid_food_cols, drop=FALSE]

# --- 3. 准备 Biomarker 矩阵 ---
results_bio_renamed <- results_bio %>% mutate(Biomarker = dplyr::recode(Biomarker, !!!biomarker_labels))
mat3_beta <- results_bio_renamed %>% dplyr::select(new_name, Biomarker, Beta) %>% pivot_wider(names_from = Biomarker, values_from = Beta) %>% column_to_rownames("new_name") %>% as.matrix()
mat3_p    <- results_bio_renamed %>% dplyr::select(new_name, Biomarker, P_value) %>% pivot_wider(names_from = Biomarker, values_from = P_value) %>% column_to_rownames("new_name") %>% as.matrix()

mat3_beta <- mat3_beta[met_order, , drop=FALSE]
mat3_p    <- mat3_p[met_order, , drop=FALSE]


mn_forest <- final_results_table_p %>%
  dplyr::select(outcome, predictor, OR, CI_low, CI_high, p.value, model_name) %>%
  filter(outcome == "t2d_confirmed_p" & model_name == "Model_2_Lifestyle" & !predictor %in% c("MRSs", "MRSq")) %>%
  left_join(msl_sub[,c("feature","new_name")], by = c("predictor" = "feature")) %>%
  distinct(new_name, .keep_all = TRUE)

forest_df <- data.frame(new_name = met_order) %>%
  left_join(mn_forest, by = "new_name")
rownames(forest_df) <- forest_df$new_name

mn_cox <- final_results_cox %>%
  filter(Model == "Model_2_Lifestyle" & !Variable %in% c("MRSs", "MRSq")) %>% 
  left_join(msl_sub[, c("feature", "new_name")], by = c("Variable" = "feature")) %>%
  distinct(new_name, .keep_all = TRUE)


cox_df <- data.frame(new_name = met_order) %>%
  left_join(mn_cox, by = "new_name")
cox_df <- cox_df %>%
  mutate(
    
    is_invalid = is.infinite(CI_Upper) | HR == 0 | is.na(HR),
    HR_plot = ifelse(is_invalid, NA, HR),
    CI_L_plot = ifelse(is_invalid, NA, CI_Lower),
    CI_U_plot = ifelse(is_invalid, NA, CI_Upper)
  )

rownames(cox_df) <- cox_df$new_name


hr_min <- min(cox_df$CI_L_plot, na.rm = TRUE) * 0.9
hr_max <- max(cox_df$CI_U_plot, na.rm = TRUE) * 1.1



map_to_x_hr <- function(val) { (val - hr_min) / (hr_max - hr_min) }


axis_anno_hr <- columnAnnotation(axis_hr = anno_empty(border = FALSE, height = unit(1, "cm")))


# ==============================================================================
# Part 5:  (Heatmaps + Forest Plot)
# ==============================================================================


cell_fun_stars = function(p_mat) {
  function(j, i, x, y, width, height, fill) {
    p_val = p_mat[i, j]
    if (!is.na(p_val) & p_val < 0.001) grid.text("***", x, y, gp = gpar(fontsize = 8))
    else if (!is.na(p_val) & p_val < 0.01) grid.text("**", x, y, gp = gpar(fontsize = 8))
    else if (!is.na(p_val) & p_val < 0.05) grid.text("*", x, y, gp = gpar(fontsize = 8))
  }
}


col_fun_1 = colorRamp2(c(-max(abs(mat1_beta), na.rm=T), 0, max(abs(mat1_beta), na.rm=T)), c("navy", "white", "firebrick3"))
col_fun_2 = colorRamp2(c(-max(abs(mat2_beta), na.rm=T), 0, max(abs(mat2_beta), na.rm=T)), c("#1b7837", "white", "#762a83"))
col_fun_3 = colorRamp2(c(-max(abs(mat3_beta), na.rm=T), 0, max(abs(mat3_beta), na.rm=T)), c("#b2182b", "white", "#2166ac"))

ht1 = Heatmap(mat1_beta, name = "Coef (MRS)", col = col_fun_1, cluster_rows = F, cluster_columns = F, 
              row_names_side = "left", row_names_gp = gpar(fontsize = 7,fontface = "italic"),
              column_names_gp = gpar(fontsize = 9, fontface = "bold"), column_names_rot = 45, width = unit(1.5, "cm"))
ht1
ht2 = Heatmap(mat2_beta, name = "Beta (Food)", col = col_fun_2, cluster_rows = F, cluster_columns = F,
              show_row_names = F, column_names_rot = 45, column_names_gp = gpar(fontsize = 9),
              cell_fun = cell_fun_stars(mat2_p), width = unit(10, "cm"))

ht3 = Heatmap(mat3_beta, name = "Beta (Bio)", col = col_fun_3, cluster_rows = F, cluster_columns = F,
              show_row_names = F, column_names_rot = 45, column_names_gp = gpar(fontsize = 9),
              cell_fun = cell_fun_stars(mat3_p), width = unit(5, "cm"))


or_min <- min(forest_df$CI_low, na.rm = TRUE) * 0.9
or_max <- max(forest_df$CI_high, na.rm = TRUE) * 1.1


map_to_x <- function(val) { (val - or_min) / (or_max - or_min) }


axis_anno <- columnAnnotation(axis = anno_empty(border = FALSE, height = unit(1, "cm")))

ht4 = Heatmap(
  matrix(forest_df$OR, ncol = 1), # 虚拟矩阵结构
  name = "T2D Risk",
  width = unit(4, "cm"),
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 9, fontface="bold"), column_names_rot = 0,
  rect_gp = gpar(type = "none"), # 隐藏背景块，变成纯白
  show_heatmap_legend = FALSE,
  bottom_annotation = axis_anno, # 占位用于画坐标轴
  cell_fun = function(j, i, x, y, width, height, fill) {
    or_val <- forest_df$OR[i]
    ci_l   <- forest_df$CI_low[i]
    ci_u   <- forest_df$CI_high[i]
    pval   <- forest_df$p.value[i]
    
    if (!is.na(or_val)) {
      
      pt_col <- ifelse(pval < 0.05, ifelse(or_val > 1, "firebrick3", "navy"), "grey50")
      
      
      left <- as.numeric(x) - as.numeric(width)/2
      w <- as.numeric(width)
      x_or <- left + w * map_to_x(or_val)
      x_l  <- left + w * map_to_x(ci_l)
      x_u  <- left + w * map_to_x(ci_u)
      x_1  <- left + w * map_to_x(1) # OR=1 的参考线
      
      
      grid.lines(c(x_1, x_1), c(as.numeric(y)-as.numeric(height)/2, as.numeric(y)+as.numeric(height)/2), 
                 gp = gpar(col = "grey70", lty = 2))
      
      
      grid.lines(c(x_l, x_u), c(y, y), gp = gpar(col = pt_col, lwd = 1.5))
      
      
      grid.points(x_or, y, pch = 16, size = unit(2, "mm"), gp = gpar(col = pt_col))
    }
  }
)
head(cox_df)
ht5 = Heatmap(
  matrix(cox_df$HR_plot, ncol = 1), 
  name = "T2D HR (Cox)",       
  width = unit(4, "cm"),
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 9, fontface="bold"), column_names_rot = 0,
  rect_gp = gpar(type = "none"), 
  show_heatmap_legend = FALSE,
  bottom_annotation = axis_anno_hr, 
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    hr_val <- cox_df$HR_plot[i]
    ci_l   <- cox_df$CI_L_plot[i]
    ci_u   <- cox_df$CI_U_plot[i]
    pval   <- cox_df$P_Value[i]
    
    
    if (!is.na(hr_val)) {
      
      ci_u_plot <- min(ci_u, hr_max) 
      ci_l_plot <- max(ci_l, hr_min)
      
      pt_col <- ifelse(!is.na(pval) & pval < 0.05, 
                       ifelse(hr_val > 1, "firebrick3", "navy"), "grey50")
      
      left <- as.numeric(x) - as.numeric(width)/2
      w <- as.numeric(width)
      x_hr <- left + w * map_to_x_hr(hr_val)
      x_l  <- left + w * map_to_x_hr(ci_l_plot)
      x_u  <- left + w * map_to_x_hr(ci_u_plot)
      x_1  <- left + w * map_to_x_hr(1) 
      
      
      grid.lines(c(x_1, x_1), c(as.numeric(y)-as.numeric(height)/2, as.numeric(y)+as.numeric(height)/2), 
                 gp = gpar(col = "grey70", lty = 2))
      
      
      grid.lines(c(x_l, x_u), c(y, y), gp = gpar(col = pt_col, lwd = 1.5))
      
      grid.points(x_hr, y, pch = 18, size = unit(2.5, "mm"), gp = gpar(col = pt_col))
    }
  }
)

ht_list = ht1 + ht2 + ht3 + ht4 + ht5
ht_list

setwd("/udd/n2gji/micro/dietpt/final/data")
pdf("MLVS_microbiome_forest_and_cox.pdf", width = 18, height = 5)
draw(ht_list, row_title = "Microbiome", merge_legend = FALSE)


decorate_annotation("axis", {
  pushViewport(viewport(xscale = c(or_min, or_max)))
  grid.xaxis(gp = gpar(fontsize = 8))
  grid.lines(x = unit(1, "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 2, col = "grey50")) 
  popViewport()
})

decorate_annotation("axis_hr", {
  pushViewport(viewport(xscale = c(hr_min, hr_max)))
  grid.xaxis(gp = gpar(fontsize = 8))
  grid.lines(x = unit(1, "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 2, col = "grey50")) 
  popViewport()
})
dev.off()
