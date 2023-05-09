################################################################################
#
# Figure 4. Relationship of microbiota to SCFA profile. 
# 
# script adapted from (with warm thanks):
# - https://dominicroye.github.io/en/2019/tidy-correlation-tests-in-r/
# - https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
# - https://stackoverflow.com/users/3063910/rich-scriven
#
# stats analysis with guidance from:
# - https://stackoverflow.com/questions/59954619/executing-a-statistical-test-across-multiple-subsets-using-purrr-map
#
################################################################################

set.seed(1216)

# packages
library("readxl")
library("broom")
library("cowplot")
library("readr")
library("dplyr")
library("tidyr")
library("stringr")
library("purrr")
library("ggplot2")
library("cowplot")
library("forcats")

# inputs
abund_file <- "results/ANCOMBC_bc_abundance.txt"
scfa_file <- "resources/zootechnical/late_scfa_concentrations.xls"

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
figpath = "submission/fig4_OTU_SCFA_correspondance.pdf"

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

# read in bias corrected taxon abundances
abund_raw <- read_tsv(abund_file)

# fix SCFA formatting and names
# covert isobutyric acid from ppb to ppm
scfa_raw <- scfa_xls %>% 
  rename(all_of(scfa_names)) %>% 
  slice(2:n()) %>% 
  select(-`Standard curve`) %>% 
  mutate(cohort = case_when(str_sub(`Sample ID`, 5,5) == "1" ~ "ctl x Mock",
                            str_sub(`Sample ID`, 5,5) == "2" ~ "jGOS x Mock",
                            str_sub(`Sample ID`, 5,5) == "3" ~ "ctl x SE",
                            str_sub(`Sample ID`, 5,5) == "4" ~ "jGOS x SE"),
         treatment = case_when(str_sub(`Sample ID`, 5,5) == "1" ~ "ctl_unc",
                               str_sub(`Sample ID`, 5,5) == "2" ~ "gos_unc",
                               str_sub(`Sample ID`, 5,5) == "3" ~ "ctl_sal",
                               str_sub(`Sample ID`, 5,5) == "4" ~ "gos_sal"),
         age = str_sub(`Sample ID`, 7,8),
         replicate = str_sub(`Sample ID`, -1L,-1L),
         Group = paste0(treatment,"_",age,"_",replicate)
  ) %>% 
  mutate_at(vars(all_of(c("Acetic acid", "Propionate", "Isobutyric acid", 
                          "Butyrate", "2-Methyl butyric acid",
                          "Isovaleric acid", "Valerate", "Hexanoic acid", 
                          "Lactic acid"))), as.numeric) %>% 
  mutate("Isobutyric acid" = `Isobutyric acid`/1000)

# clean data
# convert to mM
# fix levels
scfa_wide <- scfa_raw %>% 
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
                      age = fct_relevel(age, days))

# code check
# how many communities?
print("how many communities?")
scfa_wide %>% 
          group_by(cohort, age) %>% 
          summarise(n = n()) %>% 
          pivot_wider(values_from = n,
                      names_from = cohort)


######### conc ################################################################

# make formatting 'tidy'
scfa_long <- scfa_wide %>% 
                          select(-c(`Sample ID`,replicate,treatment)) %>% 
                          pivot_longer(-c(cohort, age, Group), 
                                       names_to = "SCFA", 
                                       values_to = "Concentration") 

# perform kruskal test for significance
kw <- scfa_long %>% 
                  select(-Group) %>% 
                  group_by(age, SCFA) %>% 
                  nest(data = c(cohort, Concentration)) %>% 
                  mutate(kruskal_raw = map(data, ~ kruskal.test(.x$Concentration, .x$cohort)),
                         kruskal = map(kruskal_raw, broom::tidy)) %>% 
                  unnest(kruskal)

# export kw stats
#kw %>% 
#      select(-c(data, kruskal_raw)) %>% 
#      write_tsv(., file = "results/kw.tsv")

pw <- kw %>% 
  filter(`p.value` < 0.05) %>% 
  select(c(age, SCFA, data)) %>%  
  mutate(wilcox_raw = map(data, ~ pairwise.wilcox.test(.x$Concentration, g = .x$cohort, 
                                                       p.adjust.method = "BH")),
         wilcox = map(wilcox_raw, tidy)) %>% 
  unnest(cols = wilcox) %>% 
  filter(`p.value` < 0.05) %>% 
  mutate(sig = case_when(`p.value` < 0.001 ~ "***",
                         `p.value` < 0.01 ~ "**",
                         `p.value` < 0.05 ~ "*"))

# export kw stats
#pw %>% 
#  select(-c(data, wilcox_raw)) %>% 
#  write_tsv(., file = "results/pw.tsv")


# get vector of SCFA with significant differences following late exposure
sig_pw <- pw %>% pull(SCFA) %>% unique()

# add placeholder/blank data to generate headspace for bars 
# indicating significance
placeholder <- tibble(SCFA = "Valerate",
                Concentration = 0.6,
                age = "22",
                cohort = "ctl x SE")

# plot SCFA concentrations
scfa_pw.p <- scfa_long %>% 
                      filter(SCFA %in% sig_pw) %>% 
                      ggplot(aes(x = age, y = Concentration, fill = cohort)) +
                        geom_bar(position="dodge", stat="summary", fun="median",
                                 colour = "black", width = 0.7)+
                        stat_summary(fun.data=median_hilow, geom = "errorbar", 
                                     position = position_dodge(width = 0.7),
                                     width = 0.2) +
                        facet_wrap(~SCFA, scales = "free_y") +
                        geom_blank(data=placeholder) + 
                        scale_y_continuous(name="Concentration (mM)",
                                           expand=expansion(c(0,0.15) )) +
                        xlab("Age (d)") +
                        scale_fill_manual(values = access_colours) +
                        theme_bw() +
                        theme(plot.margin=unit(c(5.5, 5.5, 5.5, 25), "points"),
                              panel.grid = element_blank(),
                              axis.text = element_text(colour = "black")
                              ) 


# pull out p values
# butyrate
But_GOSxSE_ctlxMock_28 <- pw %>% filter(SCFA=="Butyrate", group1=="jGOS x SE", group2=="ctl x Mock", age=="28") %>% pull(sig)

# propionate
Prop_GOSxSE_ctlxSE_28 <- pw %>% filter(SCFA=="Propionate", group1=="jGOS x SE", group2=="ctl x SE", age=="28") %>% pull(sig)
Prop_GOSxSE_ctlxMock_28 <- pw %>% filter(SCFA=="Propionate", group1=="jGOS x SE", group2=="ctl x Mock", age=="28") %>% pull(sig)
Prop_GOSxSE_GOSxMock_28 <- pw %>% filter(SCFA=="Propionate", group1=="jGOS x SE", group2=="jGOS x Mock", age=="28") %>% pull(sig)
Prop_ctlxSE_ctlxMock_28 <- pw %>% filter(SCFA=="Propionate", group1=="ctl x SE", group2=="ctl x Mock", age=="28") %>% pull(sig)
Prop_GOSxMock_ctlxSE_28 <- pw %>% filter(SCFA=="Propionate", group1=="jGOS x Mock", group2=="ctl x SE", age=="28") %>% pull(sig)

Prop_GOSxSE_ctlxMock_35 <- pw %>% filter(SCFA=="Propionate", group1=="jGOS x SE", group2=="ctl x Mock", age=="35") %>% pull(sig)
Prop_GOSxSE_GOSxMock_35 <- pw %>% filter(SCFA=="Propionate", group1=="jGOS x SE", group2=="jGOS x Mock", age=="35") %>% pull(sig)
Prop_ctlxSE_ctlxMock_35 <- pw %>% filter(SCFA=="Propionate", group1=="ctl x SE", group2=="ctl x Mock", age=="35") %>% pull(sig)
Prop_GOSxMock_ctlxSE_35 <- pw %>% filter(SCFA=="Propionate", group1=="jGOS x Mock", group2=="ctl x SE", age=="35") %>% pull(sig)

# valerate
Val_GOSxSE_ctlxMock_28 <- pw %>% filter(SCFA=="Valerate", group1=="jGOS x SE", group2=="ctl x Mock", age=="28") %>% pull(sig)
Val_GOSxSE_ctlxSE_28 <- pw %>% filter(SCFA=="Valerate", group1=="jGOS x SE", group2=="ctl x SE", age=="28") %>% pull(sig)
Val_GOSxSE_GOSxMock_28 <- pw %>% filter(SCFA=="Valerate", group1=="jGOS x SE", group2=="jGOS x Mock", age=="28") %>% pull(sig)

Val_ctlxSE_ctlxMock_35 <- pw %>% filter(SCFA=="Valerate", group1=="ctl x SE", group2=="ctl x Mock", age=="35") %>% pull(sig)
Val_GOSxMock_ctlxSE_35 <- pw %>% filter(SCFA=="Valerate", group1=="jGOS x Mock", group2=="ctl x SE", age=="35") %>% pull(sig)
Val_GOSxSE_ctlxMock_35 <- pw %>% filter(SCFA=="Valerate", group1=="jGOS x SE", group2=="ctl x Mock", age=="35") %>% pull(sig)
Val_GOSxSE_GOSxMock_35 <- pw %>% filter(SCFA=="Valerate", group1=="jGOS x SE", group2=="jGOS x Mock", age=="35") %>% pull(sig)


# add p values to plot
scfa_pw_sig.p <- 
  ggdraw(scfa_pw.p) +
  # butyrate
      draw_line(x=c(0.221,0.248), y=0.85, color="black", size=0.5) +  
        draw_label(But_GOSxSE_ctlxMock_28, x=0.236, y=0.855, size = 15) +
  # propionate
      draw_line(x=c(0.4625,0.4895), y=0.85, color="black", size=0.5) +  
        draw_label(Prop_GOSxSE_ctlxMock_28, x=0.476, y=0.855, size = 15) +
      draw_line(x=c(0.472,0.4895), y=0.81, color="black", size=0.5) +  
        draw_label(Prop_GOSxSE_ctlxSE_28, x=0.48075, y=0.815, size = 15) +
      draw_line(x=c(0.48,0.4895), y=0.77, color="black", size=0.5) +  
        draw_label(Prop_GOSxSE_GOSxMock_28, x=0.48475, y=0.775, size = 15) +
      draw_line(x=c(0.472,0.48), y=0.73, color="black", size=0.5) +  
        draw_label(Prop_GOSxMock_ctlxSE_28, x=0.476, y=0.735, size = 15) +
      draw_line(x=c(0.4625,0.472), y=0.69, color="black", size=0.5) +  
        draw_label(Prop_ctlxSE_ctlxMock_28, x=0.46725, y=0.695, size = 15) +
  
      draw_line(x=c(0.5135,0.54), y=0.85, color="black", size=0.5) +  
        draw_label(Prop_GOSxSE_ctlxMock_35, x=0.52675, y=0.855, size = 15) +
      draw_line(x=c(0.531,0.54), y=0.77, color="black", size=0.5) +  
        draw_label(Prop_GOSxSE_GOSxMock_35, x=0.5355, y=0.775, size = 15) +
      draw_line(x=c(0.5225,0.531), y=0.73, color="black", size=0.5) +  
        draw_label(Prop_GOSxMock_ctlxSE_35, x=0.52675, y=0.735, size = 15) +
      draw_line(x=c(0.5135,0.5225), y=0.69, color="black", size=0.5) +  
        draw_label(Prop_ctlxSE_ctlxMock_35, x=0.518, y=0.695, size = 15) +

  # valerate
     draw_line(x=c(0.718,0.745), y=0.85, color="black", size=0.5) +  
       draw_label(Val_GOSxSE_ctlxMock_28, x=0.7315, y=0.855, size = 15) +
    draw_line(x=c(0.727,0.745), y=0.81, color="black", size=0.5) +  
      draw_label(Val_GOSxSE_ctlxSE_28, x=0.736, y=0.815, size = 15) +
    draw_line(x=c(0.736,0.745), y=0.77, color="black", size=0.5) +  
      draw_label(Val_GOSxSE_GOSxMock_28, x=0.7405, y=0.775, size = 15) +

      draw_line(x=c(0.769,0.796), y=0.85, color="black", size=0.5) +  
        draw_label(Val_GOSxSE_ctlxMock_35, x=0.7825, y=0.855, size = 15) +
      draw_line(x=c(0.787,0.796), y=0.77, color="black", size=0.5) +  
        draw_label(Val_GOSxSE_GOSxMock_35, x=0.7915, y=0.775, size = 15) +
      draw_line(x=c(0.778,0.787), y=0.73, color="black", size=0.5) +  
        draw_label(Val_GOSxMock_ctlxSE_35, x=0.7825, y=0.735, size = 15) +
        draw_line(x=c(0.769,0.778), y=0.69, color="black", size=0.5) +  
      draw_label(Val_ctlxSE_ctlxMock_35, x=0.7735, y=0.695, size = 15) 
        
 
### microbiota x SCFA correspondance ###########################################


# check how many communities?
# a full house :)
print("how many communities?")
abund_raw %>% 
          group_by(cohort, age) %>% 
          summarise(n = n_distinct(community)) %>% 
          pivot_wider(values_from = n,
                      names_from = cohort
                      )

# format reads x taxa
abund <- abund_raw %>% 
                    pivot_wider(names_from = taxon,
                                values_from = `log(bc-abundance)`) %>% 
                    mutate(age = paste0(age," (",age-20,")"),
                           community = str_sub(community, 3,-1L) )


# bind analytic data together (microbiota & metabolome)
abund_scfa <- left_join(
                      scfa_wide %>% select(cohort, replicate, Group, all_of(sig_pw)),
                      abund %>% select(-c(pen, cohort)),
                      by = c("Group"="community")) %>% 
            pivot_longer(cols = c(Propionate, Butyrate, Valerate), 
                         names_to = "SCFA", 
                         values_to = "ppm") %>% 
            pivot_longer(-c(cohort, age, replicate, Group, SCFA, ppm), 
                         names_to = "taxon", 
                         values_to = "bc") 

abund_scfa_nest <- abund_scfa %>% 
                              group_nest(age, SCFA, taxon)


# define correlations 
abund_corr <- abund_scfa_nest %>% 
                  mutate(test = map(data, ~ cor.test(.x$bc, .x$ppm, 
                                                     method = "spearman")), 
                         tidied = map(test, tidy)      ) %>% 
                  unnest(tidied) %>% 
                  mutate(estimate_size = abs(estimate)   ) %>% 
                  group_by(age,SCFA) %>% 
                  mutate(q = p.adjust(`p.value`, method='BH'),
                         estimate_fill = if_else(q <= 0.05, estimate, 0),
                         order = str_extract(taxon, "\\d+") %>% 
                                    str_remove("^0+") %>% 
                                    as.numeric())

# get vector of taxa that show sig correlations
# at any time
keep <- abund_corr %>% 
                      filter(q <= 0.05) %>% 
                      pull(taxon) %>% 
                      unique()

# plot sig correlated taxa only
abund_corr.p <- abund_corr %>% filter(taxon %in% keep) %>% 
                    ggplot(aes(x = SCFA, y =  fct_reorder(taxon, order, .desc = TRUE))) + 
                        geom_point(aes(size=estimate_size, colour=estimate, fill=estimate_fill),
                        shape = 21) + 
                        theme_bw() +
                        theme(plot.margin=unit(c(5.5, 10, 5.5, 25), "points"),
                              panel.grid.major.x = element_blank(),
                              axis.text.x = element_text(colour = "black", size = 8, 
                                             angle = 45, hjust = 1), 
                              axis.text.y = element_text(colour = "black", size = 8), 
                              legend.text = element_text(size = 10, colour ="black")
                              ) +  
                        scale_fill_gradient2(trans = 'reverse') +
                        scale_colour_gradient2(trans = 'reverse') +
                        geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5)) +
                        facet_grid(cols = vars(age), rows = NULL, drop = FALSE) +
                        guides(fill = guide_colorbar(reverse = TRUE), 
                               colour = "none",
                               size="none") +
                        labs(fill= "Correlation\ncoefficient", y = "taxon", x = "SCFA") 


# merge plots
all.p <- plot_grid(scfa_pw_sig.p, abund_corr.p, 
                   nrow = 2, labels = c('A','B'), 
                   rel_heights = c(0.4,0.6))


# draw figure                      
pdf(file = figpath, width=7.5, height=8, paper = "a4")
all.p
dev.off()
