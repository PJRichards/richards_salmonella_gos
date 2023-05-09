################################################################################
#
# Figure S4.  Correlation plot showing significant relationships between 
# proportional OTU abundances and absolute concentration (ppm) of the 
# SCFA metabolome
# 
# script adapted from (with warm thanks):
# - https://dominicroye.github.io/en/2019/tidy-correlation-tests-in-r/
# - https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
# - https://stackoverflow.com/questions/24173194/remove-parentheses-and-text-within-from-strings-in-r/24173271 
# - https://stackoverflow.com/users/3063910/rich-scriven
#
# stats analysis with guidance from:
# - #https://stackoverflow.com/questions/59954619/executing-a-statistical-test-across-multiple-subsets-using-purrr-map
#
################################################################################

set.seed(1216)

# packages
library("readxl")
library("broom")
library("readr")
library("dplyr")
library("tidyr")
library("stringr")
library("purrr")
library("forcats")
library("ggplot2")


# inputs
shared_file <- "results/mothur/sal.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.pick.shared"
tax_file <- "results/mothur/sal.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.cons.taxonomy"
scfa_file <- "resources/zootechnical/late_scfa_concentrations.xls"

min.prevalence = 0.1
min.depth = 50  

# SCFA mol wt
propionate_MolWt <- 74.08
acetic_MolWt <- 60.05
isobutyric_MolWt <- 88.11
butyrate_MolWt <- 178.23
methylbutyric_MolWt <- 102.13
isovaleric_MolWt <- 102.13
valerate_MolWt <- 102.13
hexanoic_MolWt <- 116.16
lactic_MolWt <- 90.08
  
# output
figpath = "submission/figS4_OTU_SCFA_correspondance.pdf"

# formatting
taxonomy_list <- c("kingdom","phylum","class","order","family","genus")

scfa_names <- c(`Acetic acid` = "[Acetic Acid]", 
                Propionate = "[Propionic Acid]", 
                `Isobutyric acid` = "[isoButyric Acid]", 
                Butyrate = "[Butyric Acid]",
                `2-Methyl butyric acid` = "[2-Methyl butyric Acid]",  
                `Isovaleric acid` = "[isoValeric Acid]",  
                Valerate = "[Valeric Acid]",
                `Hexanoic acid` = "[Hexanoic Acid]", 
                `Lactic acid` = "[Lactic Acid]"
                 )         

days <- c("22","24","28","35")
access_colours = c("#f1eef6", "#df65b0", "#d7b5d8", "#ce1256")

# read in scfa data
scfa_xls <- read_excel(scfa_file, 
                       sheet = "late_scfa_-_results", 
                       range = "A2:K114")

# read in shared OTUs
OTU_raw <- read_tsv(shared_file)

# read in taxonomy data from mothur
tax_raw <- read_tsv(tax_file) 

# fix SCFA formatting and names
# covert isobutyric acid from ppb to ppm
scfa_raw <- scfa_xls %>% 
  rename(all_of(scfa_names)) %>% 
  slice(2:n()) %>% 
  select(-`Standard curve`) %>% 
  mutate(trial = str_sub(`Sample ID`,1,3),
         exposure = case_when(trial == "t42" ~ "l", 
                              trial == "t40" ~ "e"),
         cohort = case_when(str_sub(`Sample ID`, 5,5) == "1" ~ "ctl x Mock",
                            str_sub(`Sample ID`, 5,5) == "2" ~ "jGOS x Mock",
                            str_sub(`Sample ID`, 5,5) == "3" ~ "ctl x SE",
                            str_sub(`Sample ID`, 5,5) == "4" ~ "jGOS x SE"),
         treatment = case_when(str_sub(`Sample ID`, 5,5) == "1" ~ "ctl_unc",
                               str_sub(`Sample ID`, 5,5) == "2" ~ "gos_unc",
                               str_sub(`Sample ID`, 5,5) == "3" ~ "ctl_sal",
                               str_sub(`Sample ID`, 5,5) == "4" ~ "gos_sal"),
         age = case_when(trial == "t40" & str_sub(`Sample ID`, 7,7) == "0" ~ "8",
                         trial == "t40" & str_sub(`Sample ID`, 7,7) == "1" ~ "15",
                         trial == "t40" & str_sub(`Sample ID`, 7,7) == "2" ~ "22",
                         trial == "t40" & str_sub(`Sample ID`, 7,7) == "3" ~ "28",
                         trial == "t40" & str_sub(`Sample ID`, 7,7) == "4" ~ "35",
                         trial == "t42" ~ str_sub(`Sample ID`, 7,8)),
         replicate = str_sub(`Sample ID`, -1L,-1L),
         Group = paste0(exposure,"_",treatment,"_",age,"_",replicate)
  ) %>% 
  mutate_at(vars(all_of(c("Acetic acid", "Propionate", "Isobutyric acid", 
                          "Butyrate", "2-Methyl butyric acid",
                          "Isovaleric acid", "Valerate", "Hexanoic acid", 
                          "Lactic acid"))), as.numeric) %>% 
  mutate("Isobutyric acid" = `Isobutyric acid`/1000)

# clean data
# convert to mM
scfa_wide <- scfa_raw %>% 
                filter(str_detect(replicate,"\\d")) %>% 
                drop_na() %>% 
                mutate(`Acetic acid` = ((`Acetic acid`/1000)/acetic_MolWt)*1000,
                       Propionate = ((Propionate/1000)/propionate_MolWt)*1000,
                      `Isobutyric acid` = ((`Isobutyric acid`/1000)/isobutyric_MolWt)*1000,
                      `Butyrate` = ((`Butyrate`/1000)/butyrate_MolWt)*1000,
                      `2-Methyl butyric acid` = ((`2-Methyl butyric acid`/1000)/methylbutyric_MolWt)*1000,
                      `Isovaleric acid` = ((`Isovaleric acid`/1000)/isovaleric_MolWt)*1000,
                      `Valerate` = ((`Valerate`/1000)/valerate_MolWt)*1000,
                      `Hexanoic acid` = ((`Hexanoic acid`/1000)/hexanoic_MolWt)*1000,
                      `Lactic acid` = ((`Lactic acid`/1000)/lactic_MolWt)*1000,
                      age = fct_relevel(age,days))

# code check
# how many communities?
scfa_wide %>% 
          group_by(cohort, age) %>% 
          summarise(n = n()) %>% 
          pivot_wider(values_from = n,
                      names_from = cohort)  

# keep SCFA of interest
scfa_wide_key <- scfa_wide %>% 
                            select(cohort, age, replicate, Group, 
                                   Propionate, Butyrate, Valerate)


# format taxonomy data
tax <- tax_raw %>% 
                separate(Taxonomy, taxonomy_list, ";") %>% 
                mutate(id = paste0(str_replace(genus, "\\(.*\\)", ""), " (", OTU, ")")) %>%  
                select(-Size)

# read in OTU data
shared_raw <- read_tsv(file = shared_file, col_names = TRUE)

# drop superfluous columns
shared <- shared_raw %>% 
                       select(-c(label, numOtus))


# drop columns with min reads below cutoff
abund_shared <- shared %>% 
                    select(where(~ any(.x >= min.depth)))


# drop OTUs present in <10% communities 
prev_shared <- abund_shared %>% 
                            select(which(colSums(. != 0)>=nrow(.)*min.prevalence))

# merge taxonomy
# calculate sequencing depth
prev_shared_taxon <- prev_shared %>% 
                          pivot_longer(cols = -Group,
                                       names_to = "OTU",
                                       values_to = "reads") %>% 
                          pivot_wider(names_from = Group,
                                      values_from = reads) %>% 
                          left_join(tax) %>% 
                          bind_rows(summarise(.,
                                          across(where(is.numeric), sum),
                                          across(where(is.character), ~"Total")))

# make data 'tidy'
prev_shared_taxon_long <- prev_shared_taxon %>% 
                                          pivot_longer(cols = where(is.double),
                                                       names_to = "Group",
                                                       values_to = "reads") 

# how many communities?
prev_shared_taxon_long %>% 
                mutate(variable = str_sub(Group, 3, 12)) %>% 
                group_by(variable) %>% 
                summarise(n = n_distinct(Group)) 

# get relative abundance
OTU <- prev_shared_taxon_long %>% 
                      select(-c(OTU, kingdom, phylum, class, order, family, genus)) %>% 
                      pivot_wider(names_from = id,
                                  values_from = reads) %>% 
                      rowwise() %>% 
                      mutate(across(where(is.numeric),  ~100*(. / Total))) %>% 
                      select(-Total)

# code check
# does RA add up to 100%?
OTU %>% 
        ungroup() %>% 
        mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>% 
        pull(Total)
        
# bind analytic data together (microbiota & metabolome)
OTU_scfa <- left_join(scfa_wide_key, OTU, by = "Group") %>% 
                pivot_longer(cols = c(Propionate, Butyrate, Valerate), 
                             names_to = "SCFA", 
                             values_to = "ppm") %>% 
                pivot_longer(-c(cohort, age, replicate, Group, SCFA, ppm), 
                             names_to = "OTU", 
                             values_to = "RA") 

OTU_scfa_nest <- OTU_scfa %>% 
                            group_nest(age, SCFA, OTU)

# calculate correlations
OTU_corr <- OTU_scfa_nest %>% 
                      mutate(test = map(data, ~ cor.test(.x$RA, .x$ppm, 
                                                         method = "spearman")), 
                            tidied = map(test, tidy)      ) %>% 
                      unnest(tidied, .drop = TRUE) %>% 
                      mutate(estimate_size = abs(estimate)   ) %>% 
                      group_by(age,SCFA) %>% 
                      mutate(q = p.adjust(`p.value`, method='BH'),
                             estimate_fill = if_else(q <= 0.05, estimate, 0),
                             order = str_extract(OTU, "\\d+") %>% 
                                              str_remove("^0+") %>% 
                                                as.numeric(),
                             age = paste0(age," (",as.numeric(as.character(age))-20,")")
                             )

# plot correlation that meet 0.05 probability threshold
p <- OTU_corr %>% filter(q <= 0.05) %>% 
              ggplot(aes(x = SCFA, y = fct_reorder(OTU, order, .desc = TRUE))) + 
                    geom_point(aes(size=estimate_size, colour=estimate, 
                                   fill=estimate_fill), shape = 21) +
                    facet_grid(cols = vars(age), rows = NULL, drop = FALSE)

# format plot
p.fmt <- p +
            theme_bw() +
            theme(plot.margin=unit(c(5.5, 5.5, 5.5, 25), "points"),
                  panel.grid.major.x = element_blank(),
                  axis.text.x = element_text(colour = "black", size = 8, 
                                             angle = 45, hjust = 1), 
                  axis.text.y = element_text(colour = "black", size = 8), 
                  legend.text = element_text(size = 10, colour ="black")
                  ) +  
            scale_fill_gradient2(trans = 'reverse') +
            scale_colour_gradient2(trans = 'reverse') +
            geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5)) +
            guides(fill = guide_colorbar(reverse = TRUE), 
                   colour = "none",
                   size="none") +
            labs(fill= "Correlation\ncoefficient", y = "OTU", x = "SCFA") 


# draw figure                      
pdf(file = figpath, width=7.5, height=9, paper = "a4")
p.fmt 
dev.off()
