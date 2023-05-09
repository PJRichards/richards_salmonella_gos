########################################################################
#
# Figure S2. Effect of *S.* Enteritidis challenge on expression 
# of genes connected to innate immune response.
#
########################################################################

set.seed(1216)

# packages
library("readxl")
library("dplyr") 
library("tidyr") 
library("ggplot2")
library("cowplot")

# input
# and paths
text_size <- 7

GE_path <- "resources/zootechnical/SE_challenge_host_gene_expression.xlsx"
GE_sheet <- "T42 Cecal tonsil"

# output
figpath <- "submission/figS2_hostinnate.pdf"

# read in data
raw_GE <- read_excel(path = GE_path, sheet = GE_sheet)

# set column name for metadata
group_meta_id <- c("Trial", "Feed", "Challenge", "Age", "Replicate")

# Genes of interest
otherimm <- c("IL4", "STAT6", "STAT5B", "STAT3", "IL10", "NOS2", "STAT4", 
              "GATA3", "IL15", "TGIF1")

# format data
# make 'tidy'
GE <- raw_GE %>% 
              rename(group = '...1',
                     'IFN-gamma' = IFNG,
                     'IL-1*Beta' = IL1B,
                     'IL-4' = IL4,
                     'IL-10' = IL10,
                     'IL-15' = IL15) %>%
              pivot_longer(names_to = "gene", 
                           values_to = "log10(2-^dCt)",
                           -group) %>% 
              separate(group, group_meta_id, "_", remove = FALSE) %>%
              select(-Trial)

# fix 'feed' and 'challenge' names
GE <- GE %>%  
        mutate(Feed = if_else(Feed == "ctl", 
                        true = "Ctl", 
                        false = if_else(Feed == "gos", 
                                        true = "jGOS", 
                                        false = "fail")),
              Challenge = if_else(Challenge == "sal", 
                                  true = "SE", 
                                  false = if_else(Challenge == "unc", 
                                                  true = "Mock", 
                                                  false = "fail"))
                                  )

# invert log10 expression values
# drop NA values
alDE <-  GE %>% 
              mutate(cohort = paste0(Feed,"x",Challenge),
                     `2-^dCt` = 10^`log10(2-^dCt)`) %>% 
              drop_na(`2-^dCt`)

# Observations were grouped as lists according to both feed x challenge group 
# and age and the mean 2^-$\Delta$Ct^ levels for expression observed in 
# mock-challenged birds was determined for groups sustained on either control or
# jGOS-supplemented feed
alDE_wide <- alDE %>% 
                    pivot_wider(names_from = cohort,
                                values_from = `2-^dCt`,
                                -c(group,Replicate,Feed,Challenge,`log10(2-^dCt)`),
                                values_fn = list) %>% 
                    rowwise() %>% 
                    mutate(m.CtlxMock = mean(unlist(CtlxMock)),
                           m.jGOSxMock = mean(unlist(jGOSxMock)))

# calculate fold-change relative to mean, 
# The fold-change value was then calculated by dividing the 2^-$\Delta$Ct^ value 
# observed in each bird by the mean 2^-$\Delta$Ct for mock-challenged birds on the 
# same diet. 
FC <- alDE_wide %>% 
                mutate(FC.CtlxSE = list(CtlxSE/m.CtlxMock),
                       FC.CtlxMock = list(CtlxMock/m.CtlxMock),
                       FC.jGOSxSE = list(jGOSxSE/m.jGOSxMock),
                       FC.jGOSxMock = list(jGOSxMock/m.jGOSxMock)
                       )  %>% 
                select(-c(CtlxMock, jGOSxMock))  

FC.n <- FC %>% 
              mutate(FC.CtlxSE.n = length(unlist(FC.CtlxSE)),
                     FC.jGOSxSE.n = length(unlist(FC.jGOSxSE))) %>% 
              group_by(gene) %>% 
              count(FC.CtlxSE.n,FC.jGOSxSE.n) 

keep <- FC.n %>% filter(n == 4) %>% pull(gene) %>% unique()

# keep only genes with n = 7  
# get stats
FC.stats <- 
    FC %>% 
        filter(gene %in% keep) %>% 
        select(c(Age, gene, FC.CtlxSE, FC.jGOSxSE)) %>% 
        group_by(Age, gene) %>% 
        mutate(p_value = t.test(unlist(FC.CtlxSE), unlist(FC.jGOSxSE))$p.value,
               t_value = t.test(unlist(FC.CtlxSE), unlist(FC.jGOSxSE))$statistic) %>%                      
        group_by(Age) %>% 
        mutate(p_adjust = p.adjust(p_value, method = "BH")) %>% 
        rowwise() %>% 
        mutate(FC.CtlxSE.n = length(unlist(FC.CtlxSE)),
               FC.jGOSxSE.n = length(unlist(FC.jGOSxSE)),
               sig = if_else(p_adjust < 0.001, "***", 
                          if_else(p_adjust < 0.01, "**",
                               if_else(p_adjust < 0.05, "*", "ns"))))

# make data 'long' for plotting
FC.long <- FC %>%  
                select(-c(CtlxSE, jGOSxSE, m.CtlxMock, m.jGOSxMock)) %>% 
                pivot_longer(names_to = "Feed",
                             values_to = "expression",
                             cols = c("FC.CtlxSE", "FC.CtlxMock", 
                                      "FC.jGOSxSE", "FC.jGOSxMock")) %>% 
                unnest(cols = c(expression)) %>% 
                filter(Feed != "FC.CtlxMock", Feed != "FC.jGOSxMock")
  

# plot gene expression data for cecal
# annotate with significant differences
GE_SE.p <- FC.long %>% 
                    filter(gene %in% otherimm) %>% 
                    ggplot(aes(x = Age, y = expression, fill = Feed)) +
                      geom_boxplot(colour = "black", outlier.shape = NA, lwd = 0.3) +
                      geom_point(position=position_jitterdodge(jitter.width = 0.15),
                                 size =0.4) +
                      geom_hline(yintercept=1, size=0.3) +
                      theme_bw() + 
                      theme(legend.position = "bottom") +
                      facet_wrap(~gene, scales = "free_y", ncol = 3, 
                                 labeller = label_parsed) +
                      scale_y_continuous(name="Fold-change relative to\nmock-challenged birds", 
                                         expand=expansion(c(0,0.3) ) ) +
                      scale_shape_manual(values=c(19, 17))+
                      scale_fill_manual(values = c("#ffffff", "#f5e4ae")) +
                      theme(axis.text = element_text(colour = "black", size = text_size),
                            legend.position = "bottom",
                            panel.grid.minor = element_blank(), 
                            panel.grid.major = element_blank()
                            )


GE_SE.annot.p <- ggdraw(GE_SE.p) + 
  draw_line(x = c(0.742, 0.767), y=0.92, color="black", size=0.4) +
    draw_label(FC.stats %>% filter(Age==22, gene=="STAT3") %>% pull(p_adjust) %>% round(.,3), 
               size=7, x=0.7545, y=0.93) +
  draw_line(x = c(0.123, 0.148), y = 0.34, color = "black", size = 0.4) +
    draw_label(FC.stats %>% filter(Age==22, gene=="TGIF1") %>% pull(p_adjust) %>% round(.,3), 
              size=7, x = 0.1355, y = 0.35)


# export plots to pdf
pdf(file = figpath, paper = "a4", width=6.5, height = 8)
GE_SE.annot.p
dev.off()
