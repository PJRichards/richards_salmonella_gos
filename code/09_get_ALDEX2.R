#################################################################
#
# Discriminative OTUs determined using ALDEx2
#
#################################################################

# packages  
library("readr")
library("stringr")
library("dplyr")
library("tidyr")
library("tibble")
library("purrr")
library("ALDEx2")

set.seed(1216)

min.prevalence = 0.1
min.depth = 50  

# paths
shared.path <- "results/mothur/sal.trim.contigs.sort.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.pick.shared"
tax.path <- "results/mothur/sal.trim.contigs.sort.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.cons.taxonomy"

# output
tsvpath <- "results/aldex2_glmtest.tsv"
figpath <- "results/aldex2_plot.pdf"
sessioninfopath <- "results/aldex2_sessioninfo.txt"


# read in data
shared.raw <- read_tsv(file = shared.path, col_names = TRUE)
tax.raw <- read_tsv(file = tax.path)


shared <- shared.raw %>% 
                      dplyr::select(-c(label, numOtus)) 
              

# drop columns with min reads cutoff
abund.shared <- shared %>% 
                    dplyr::select(where(~ any(.x >= min.depth)))


# drop OTUs present in 10% communities or below
prev.shared <- abund.shared %>% 
                      dplyr::select(which(colSums(. != 0)>=nrow(.)*min.prevalence))
          

# change orientation for analysis
prev.shared.long  <- prev.shared %>% 
                                   pivot_longer(cols = -Group,names_to = "OTU") %>% 
                                   pivot_wider(names_from = Group, values_from = value) 

# get taxonomy
tax <- prev.shared.long %>%
                            dplyr::select(OTU) %>% 
                            left_join(tax.raw %>% dplyr::select(-Size)) %>% 
                            separate(col = Taxonomy, sep = ";", remove = TRUE, 
                                     into = c("kingdom", "phylum", "class", "order",
                                              "family", "genus")) %>% 
                            mutate(species = paste0(str_extract(string = genus, 
                                   pattern = "^[^\\(]+"), " (", OTU, ")")) %>% 
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
                             dplyr::select(c(diet, challenge, age))
                        


# get matrix
matrix <- prev.shared.long %>% 
                              column_to_rownames(var = "OTU") %>% 
                              as.matrix()


## late DAA
mm <- model.matrix(~ challenge + diet + age, design)

x <- aldex.clr(matrix, 
               mm, 
               mc.samples=128, denom="all", verbose=T)

glm.test <- aldex.glm(x, mm)

glm.effect <- aldex.glm.effect(x)


# export outputs
pdf(file = figpath, paper = "a4r")
ALDEx2::aldex.plot(glm.effect[["challengeSE"]], test="effect", cutoff=2)
sig <- glm.test[,20]<0.05
points(glm.effect[["challengeSE"]]$diff.win[sig],
       glm.effect[["challengeSE"]]$diff.btw[sig], col="blue")
dev.off()

write_tsv(x = as_tibble(glm.test, rownames = "OTU"), file = tsvpath)
writeLines(capture.output(sessionInfo()), sessioninfopath)


