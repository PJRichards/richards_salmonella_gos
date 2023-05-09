#########################################################################
#
# Figure S1. Sequencing effort for all communities.
#
#########################################################################

# packages
library("readr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("scales")
library("stringr")
library("RColorBrewer")

set.seed(1216)

# input
raref.path <- "results/mothur/sal.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.groups.rarefaction"

group_meta_id <- c("Trial", "Diet", "Challenge", "Age", "Replicate")

# output
fig_path <- "submission/figS1_rarefaction.pdf"

# read in data
raref_raw <- read_tsv(file=raref.path)


raref <- raref_raw %>% 
                  select(numsampled, starts_with("0"))


# format names
names(raref)[2:length(raref)] <- substring(names(raref), 6)[2:length(raref)]


# tidy tibble
# exclude NTC
raref_tidy <- raref %>% 
                      pivot_longer(names_to = "community", 
                                   values_to = "OTUs", 
                                   -numsampled) %>%
                      drop_na(OTUs) |> 
                      filter(str_detect(community,
                                        pattern = "^l_[a-z]{3}_[a-z]{3}_\\d{2}_\\d"))

# extract metadata
raref_meta <- raref_tidy %>% 
                    separate(community, group_meta_id, "_", remove = FALSE) %>% 
                    mutate(Cohort = case_when(paste(Diet, Challenge) == "ctl unc" ~ "ctl x Mock",
                                              paste(Diet, Challenge) == "gos unc" ~ "jGOS x Mock",
                                              paste(Diet, Challenge) == "ctl sal" ~ "ctl x SE",
                                              paste(Diet, Challenge) == "gos sal" ~ "jGOS x SE"),
                           Age = paste0(Age, " (", as.numeric(Age)-20,")")
                          )
  

# plot data
p <- raref_meta %>% 
                  filter(Trial == "l") %>% 
                  ggplot() +
                        geom_line(aes(x=numsampled / 1000, y = OTUs, 
                                      group = community, colour = Replicate)) +
                        theme_bw() + 
                        scale_color_brewer(palette = "Dark2") +
                        theme(plot.margin = unit(c(0.5, 0, 0, 1.75), "cm"),
                              axis.text = element_text(size=8)
                              ) +
                        labs(y = "Number of OTUS", x = "Number of reads sampled (1000s)") +
                        facet_grid(vars(Age), vars(Cohort)) +
                        ylim(0, 300)





# print plot
pdf(fig_path, paper = "a4r", height = 6, width = 9)
p 
dev.off()




