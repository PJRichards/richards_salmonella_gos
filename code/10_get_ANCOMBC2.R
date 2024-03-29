################################################################################
# 
# get ANCOM-BC2 analysis
# (adapted from Huang Lin 
# https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
# Thank you
#
################################################################################

set.seed(1216)

# packages
library("ANCOMBC")
library("readr")
library("tidyr")
library("dplyr")
library("stringr")
library("TreeSummarizedExperiment")
library("tibble")
library("ggplot2")
library("cowplot")
library("purrr")

# set variables
min.prevalence = 0.1    # prevalence filtering
min.depth = 50          # community depth filtering
tax_level = "species"

# paths
shared.path <- "results/mothur/sal.trim.contigs.sort.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.pick.shared"
tax.path <- "results/mothur/sal.trim.contigs.sort.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.cons.taxonomy"

# output path
outpath <- "results/ANCOMBC_"

### Functions ############################

fig_fix <- function(df, diff_var, lfc_var) {
  df %>%
    filter({{diff_var}} == 1) %>% 
    arrange(desc({{lfc_var}})) %>%
    mutate(direct = ifelse({{lfc_var}} > 0, "Positive LFC","Negative LFC"),
           taxon = factor(taxon, levels = .$taxon),
           direct = factor(direct, levels = c("Positive LFC", "Negative LFC"))
    )
}

fig_plot <- function(df, lfc_var, se_var, title_var) {
  df %>%
    ggplot(aes(x = taxon, y = {{lfc_var}}, fill = direct)) + 
    geom_bar(stat = "identity", width = 0.7, color = "black", 
             position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = {{lfc_var}} - {{se_var}}, 
                      ymax = {{lfc_var}} + {{se_var}}), 
                  width = 0.2, position = position_dodge(0.05), 
                  color = "black") + 
    labs(x = NULL, y = "Log fold change", title = title_var) + 
    scale_fill_discrete(name = NULL) +
    scale_color_discrete(name = NULL) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 80, hjust = 1, size = 6))
}

sens_fix <- function(sens_df, df, diff_var, newname, oldname) {
  sens_df %>%
    dplyr::rename({{newname}} := {{oldname}}) %>% 
    select(taxon, {{newname}}) %>%
    left_join(df, by = "taxon") %>% 
    mutate({{diff_var}} := recode({{diff_var}} * 1, 
                                  `1` = "Significant",
                                  `0` = "Nonsignificant"))
}

sens_plot <- function(df, sens_var, diff_var) {
  df %>%
    ggplot(aes(x = taxon, y = {{sens_var}}, color = {{diff_var}})) +
    geom_point() +
    scale_color_brewer(palette = "Dark2", name = NULL) +
    labs(x = NULL, y = "Sensitivity Score") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2))
}

log_corr.fix <- function(df) {
  as_tibble(x = df, rownames="taxon") %>% 
    pivot_longer(cols = -taxon, 
                 names_to = "community", 
                 values_to = "log(bc-abundance)") %>% 
    mutate(pen = paste(str_sub(community,3,10),str_sub(community,-1L,-1L)),
           age = str_sub(community,11,12 ),
           cohort = str_sub(community,3,9))
}

# read in data
shared.raw <- read_tsv(file = shared.path, col_names = TRUE)
tax.raw <- read_tsv(file = tax.path)

# drop unnecessary 'label' and 'numOtus' columns
shared <- shared.raw %>% 
                      select(-c(label, numOtus)) 

# drop columns with min reads cutoff
abund.shared <- shared %>% 
                        select(where(~ any(.x >= min.depth)))

# drop OTUs present in 10% communities or below
prev.shared <- abund.shared %>% 
                          select(which(colSums(. != 0)>=nrow(.)*min.prevalence))


# change orientation for analysis
prev.shared.long  <- prev.shared %>% 
                            pivot_longer(cols = -Group,names_to = "OTU") %>% 
                            pivot_wider(names_from = Group, values_from = value)

# get taxonomy
tax <- prev.shared.long %>%
  select(OTU) |> 
  left_join(tax.raw |> select(-Size)) %>% 
  separate(col = Taxonomy, sep = ";", remove = TRUE, 
           into = c("kingdom", "phylum", "class", "order",
                    "family", "genus")) |> 
  mutate(species = paste0(str_extract(string = genus, 
                                      pattern = "^[^\\(]+"), " (", OTU, ")")) |> 
  column_to_rownames(var = "OTU")

# get metadata
design <- prev.shared %>% 
  separate(col = Group, sep = "_", remove = FALSE,
           into = c("trial", "diet", "challenge", "age")) %>% 
  mutate(age = as.integer(age),
         diet = as.factor(diet),
         challenge = if_else(challenge == "sal", "SE",
                             if_else(challenge == "unc", "mock", "fail")),
         challenge = as.factor(challenge)) %>% 
  column_to_rownames(var = "Group") %>% 
  select(c(diet, challenge, age))

# get matrix
matrix <- prev.shared.long %>% 
                            column_to_rownames(var = "OTU") %>% 
                            as.matrix()

# Continuous covariates “age”
# Categorical covariates: “diet”, “challenge”

# The group variable of interest: “challenge”
# Two groups: “SE”, “Mock”
# The reference group: “SE”

# construct of TreeSummarizedExperiment
tse <- TreeSummarizedExperiment(assays=list(counts=matrix),
                                rowData=tax,
                                colData=design
                                )

# Run ancombc2 function on 'late' trial
output <- ancombc2(data = tse, 
                   assay_name = "counts", tax_level = tax_level,
                   fix_formula = "age + diet + challenge", 
                   rand_formula = NULL,
                   p_adj_method = "holm", 
                   pseudo = 0, 
                   pseudo_sens = TRUE,
                   prv_cut = 0, 
                   lib_cut = 0, 
                   s0_perc = 0.05,
                   group = "challenge", 
                   struc_zero = TRUE, 
                   neg_lb = TRUE,
                   alpha = 0.05, 
                   n_cl = 4, 
                   verbose = TRUE,
                   global = FALSE, 
                   pairwise = TRUE, 
                   dunnet = TRUE, 
                   trend = FALSE,
                   iter_control = list(tol = 1e-2, max_iter = 20, 
                                       verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl(),
                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100)
                    )

# Structural zeros
tab_zero <- output$zero_ind

# Sensitivity scores
#tab_sens <- output$pseudo_sens_tab

# ANCOM-BC2 primary analysis
res_prim <- output$res

### Results for age ############################################################

df_age <- res_prim %>% select(taxon, ends_with("age")) 


df_fig_age <- fig_fix(df_age, diff_age, lfc_age)

fig_age.p <- fig_plot(df_fig_age, lfc_age, se_age, 
                      "Log fold changes as one unit increase of age")

# output
write_tsv(x = df_age, file = paste0(outpath,"age.tsv"))

pdf(file = paste0(outpath,"age_plot.pdf"), paper = "a4", height = 9, width = 8)
fig_age.p
dev.off()

### Results for diet ######################################################

# Original results
df_diet <- res_prim %>% select(taxon, contains("diet"))

df_fig_diet <- fig_fix(df_diet, diff_dietgos, lfc_dietgos)

fig_diet.p <- fig_plot(df_fig_diet, lfc_dietgos, se_dietgos, 
                           "Log fold changes associated with jGOS diet")

# outputs
write_tsv(x = df_diet, file = paste0(outpath,"diet.tsv"))

pdf(file = paste0(outpath,"diet_plot.pdf"), paper = "a4", height = 9, width = 8)
fig_diet.p 
dev.off()


### Results for challenge ######################################################

# Original results
df_challenge <- res_prim %>% select(taxon, contains("challenge")) 


df_fig_challenge <- fig_fix(df_challenge, diff_challengeSE, lfc_challengeSE)


fig_challenge.p <- fig_plot(df_fig_challenge, lfc_challengeSE, se_challengeSE, 
                                 "Log fold changes associated with SE challenge")

# outputs
write_tsv(x = df_challenge, file = paste0(outpath,"challenge.tsv"))

pdf(file = paste0(outpath,"challenge_plot.pdf"), paper = "a4", height = 9, width = 8)
fig_challenge.p
dev.off()


### Bias-corrected abundances##################################################
samp_frac <- output$samp_frac

# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 

# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn <- log(output$feature_table + 1)

# Adjust the log observed abundances
log_corr_abn <- t(t(log_obs_abn) - samp_frac)

log_corr_abn_df <- log_corr.fix(log_corr_abn)

# export
write_tsv(x = log_corr_abn_df, file =  paste0(outpath, "bc_abundance.txt"))

writeLines(capture.output(sessionInfo()), paste0(outpath, "sessioninfo.txt"))
