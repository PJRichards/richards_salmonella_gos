###############################################################################
#
# Figure 3. A distinct microbiota distinguishes *S.* Enteritidis-challenged 
# chickens raised on a GOS-supplemented diet
#
###############################################################################


# packages
library("readr")
library("cowplot")
library("FSA")
library("readxl")
library("dplyr")
library("tidyr")
library("ggplot2")
library("stringr")
library("purrr")
library("ggrepel")
library("forcats")

set.seed(1216)

# input
# and paths
text_size <- 6
header_size <- 8

GE_path <- "resources/zootechnical/SE_challenge_host_gene_expression.xlsx"
GE_sheet <- "T42 Cecal tonsil"

pcoa_path <- "results/mothur/sal.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.pick.braycurtis.0.03.lt.ave.pcoa.axes"
pcoa_loadings_path <- "results/mothur/sal.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.pick.braycurtis.0.03.lt.ave.pcoa.loadings"
alpha_path <- "results/mothur/sal.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.pick.groups.ave-std.summary"

aldex2_path <- "results/aldex2_glmtest.tsv"
ancombc_diet_path <- "results/ANCOMBC_diet.tsv"
ancombc_chal_path <- "results/ANCOMBC_challenge.tsv"

# output
figpath <- "submission/fig3_hostmicrobe_test.pdf"

# formatting
group_meta_id <- c("Trial", "Feed", "Challenge", "Age", "Replicate")
cohort_levels <- c("ctl unc","ctl sal","gos unc", "gos sal")
sal_imm <- c("IFN-gamma","'IL-17A'","IL-1*Beta","'IL-22'","'IL-8L1'","'IL-8L2'")


############ Bray Curtis PCOA ################

# read in data
pcoa_loadings <- read_tsv(pcoa_loadings_path)
pcoa_raw <- read_tsv(pcoa_path)

# format data
pcoa <- pcoa_raw %>% 
                separate(group, group_meta_id, "_") %>% 
                mutate(Feed = if_else(Feed == "gos", "jGOS", Feed),
                       Challenge = if_else(Challenge == "sal", "SE", 
                                     if_else(Challenge == "unc", "Mock", "fail")),
                       Cohort = paste(Feed,"x", Challenge),
                       Age = paste0(Age," (",as.numeric(Age)-20,")"))
                       
# get n
print("get n")
pcoa %>% group_by(Trial,Challenge,Age,Feed) %>% count()

# plot figure
pcoa.p <- 
  pcoa %>% 
  ggpubr::ggscatter(x = "axis1", y = "axis2",
                    color = "Cohort", shape = "Cohort", 
                    palette = c("#000000","#995c5c","#000000","#995c5c"),
                    ellipse = TRUE, ellipse.type = "confidence",
                    mean.point = FALSE, star.plot = FALSE,
                    ggtheme = theme_bw()) + 
  scale_shape_manual(values = c(16,16,17,17)) +
  facet_grid(cols = vars(Age)) +
  theme(aspect.ratio = 1,
        axis.text.y = element_text(colour = "black", size = text_size),
        axis.text.x = element_text(colour = "black", size = text_size),
        axis.title = element_text(colour = "black", size = header_size),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text.x = element_text(size = header_size),
        legend.title = element_text(size=header_size),
        legend.text = element_text(size=header_size)
  ) +
  xlab(paste0("Axis 1 (",pcoa_loadings %>% 
                filter(axis == 1) %>% 
                pull(loading) %>% 
                round(., 1),"%)")) +
  ylab(paste0("Axis 2 (",pcoa_loadings %>% 
                filter(axis == 2) %>% 
                pull(loading) %>% 
                round(., 1),"%)"))


############ Alpha Diversity ################

# read in data
alpha_raw <- read_tsv(alpha_path)

# format data
# order factors
alpha <- alpha_raw %>% 
                      filter(method == "ave") %>% 
                      separate(group, group_meta_id, "_") %>% 
                      mutate(Cohort = paste(Feed, Challenge),
                             Cohort = fct_relevel(Cohort, cohort_levels),
                             Diet = if_else(Feed == "gos", "GOS", Feed),
                             Challenge = if_else(Challenge == "sal", "SE", 
                                            if_else(Challenge == "unc", "Mock", 
                                                    "fail")),
                             feed = if_else(Diet == "ctl", "ctl",
                                            if_else(Diet == "GOS" & Age < 28,"+GOS","ctl")),
                             Age = paste0(Age," (",as.numeric(Age)-20,")"))


# are data normally distributed?
print("# are data normally distributed?")
alpha %>% 
          select(Feed, Challenge, Age, invsimpson) %>% 
          group_by(Feed, Challenge, Age) %>% 
          summarise(shapiro = shapiro.test(invsimpson)$p.value)


# get stats
dunn <- alpha %>% 
                group_nest(Age) %>% 
                mutate(dunntest = map(data, ~.x %>% dunnTest(invsimpson ~ Cohort, 
                                                             data=., method="bh")))


# plot alpha
invsimp.p <- alpha %>% 
                      ggplot(aes(x = Cohort, y = invsimpson, fill = feed)) +
                        geom_boxplot(outlier.shape = NA,lwd=0.3) +
                        geom_point(aes(shape = Challenge), size = 1,
                                   position = position_jitterdodge(0.2)) +
                        scale_y_continuous(limits = c(0, 35), 
                                           breaks = seq(0, 35, by = 5), 
                                           expand = c(0,0)) +
                        scale_x_discrete(labels = c("Mock","SE","Mock","SE")) +
                        facet_grid(cols = vars(Age)) +
                        theme_bw() +
                        theme(strip.text.x = element_text(size = header_size),
                              plot.margin=margin(t=5.5,r=15, b=5.5,l=20),
                              axis.text.y = element_text(colour = "black", 
                                                         size = text_size),
                              axis.text.x = element_text(colour = "black", 
                                                         size = text_size),
                              axis.title = element_text(colour = "black", 
                                                        size = header_size),
                              panel.grid.minor = element_blank(), 
                              panel.grid.major = element_blank(),
                              legend.position = "None")+
                        scale_fill_manual(values = c("#f5e4ae","#ffffff")) +
                        xlab("") + ylab("Inverse Simpson index") 


# print plot
invsimp_annotate.p <- ggdraw(invsimp.p) +
  draw_line(x = c(0.112, 0.165), y = 0.12, color = "black", size = 0.35) +
    draw_label("ctl", x = 0.1385, y = 0.08, size = 7) +
  draw_line(x = c(0.214, 0.267), y = 0.12, color = "black", size = 0.35) +
    draw_label("jGOS", x = 0.2405, y = 0.08, size = 7) +
  
  draw_line(x = c(0.337, 0.39), y = 0.12, color = "black", size = 0.35) +
    draw_label("ctl", x = 0.3635, y = 0.08, size = 7) +
  draw_line(x = c(0.439, 0.492), y = 0.12, color = "black", size = 0.35) +
    draw_label("jGOS", x = 0.4655, y = 0.08, size = 7) +
  
  draw_line(x = c(0.562, 0.615), y = 0.12, color = "black", size = 0.35) +
    draw_label("ctl", x = 0.5885, y = 0.08, size = 7) +
  draw_line(x = c(0.664, 0.717), y = 0.12, color = "black", size = 0.35) +
  draw_label("jGOS", x = 0.6905, y = 0.08, size = 7) +
  
  draw_line(x = c(0.787, 0.84), y = 0.12, color = "black", size = 0.35) +
    draw_label("ctl", x = 0.8135, y = 0.08, size = 7) +
  draw_line(x = c(0.889, 0.942), y = 0.12, color = "black", size = 0.35) +
    draw_label("jGOS", x = 0.9155, y = 0.08, size = 7) +
  
  # ctl sal - gos sal
  draw_line(x = c(0.84, 0.942), y = 0.745, color = "black", size = 0.35) +
    draw_label(round(dunn[[4,3]][[1]]$res[2,4],3), x = 0.891, y = 0.78, size = text_size) +
  # gos sal - gos unc
  draw_line(x = c(0.889, 0.942), y = 0.595, color = "black", size = 0.35) +
  draw_label(round(dunn[[4,3]][[1]]$res[6,4],3), x = 0.9155, y = 0.63, size = text_size) +
  # ctl unc - gos sal
  draw_line(x = c(0.787, 0.942), y = 0.67, color = "black", size = 0.35) +
  draw_label(round(dunn[[4,3]][[1]]$res[3,4],3), x = 0.8645, y = 0.705, size = text_size)


############ GE ################################################

# read in data
raw_GE <- read_excel(path = GE_path, sheet = GE_sheet)

# format data
# make 'tidy'
GE <- raw_GE %>% 
              rename(group = '...1', 
                    'IFN-gamma' = IFNG,
                    'IL-1*Beta' = IL1B,
                    "'IL-17A'" = IL17A,
                    "'IL-22'" = IL22,
                    "'IL-8L1'" = IL8L1,
                    "'IL-8L2'" = IL8L2) %>%
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
# and age 
# mean 2^-$\Delta$Ct^ levels for expression observed in mock-challenged birds 
# was determined for groups sustained on either control or jGOS-supplemented feed
alDE_wide <- alDE %>% 
                  pivot_wider(names_from = cohort,
                              values_from = `2-^dCt`,
                              -c(group,Replicate,Feed,Challenge,`log10(2-^dCt)`),
                              values_fn = list) %>% 
                  rowwise() %>% 
                  mutate(m.CtlxMock = mean(unlist(CtlxMock)),
                         m.jGOSxMock = mean(unlist(jGOSxMock)))


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

# keep only genes with n = 7 for all groups
keep <- FC.n %>% filter(n == 4) %>% pull(gene) %>% unique()

FC.trim <- FC %>% filter(gene %in% keep)

# get n 
print("# what is n for jGOS x SE?")
lengths(FC.trim$jGOSxSE) %>% unique()

print("# what is n for Control x SE?")
lengths(FC.trim$CtlxSE) %>% unique()

# get stats
FC.stats <- 
  FC.trim %>% 
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
              filter(Feed != "FC.CtlxMock", Feed != "FC.jGOSxMock") %>% 
              mutate(Feed = str_extract(Feed, "(?<=\\.).[a-tA-t]*"))


# plot gene expression data for cecal tissue
# annotate with significant differences
GE_SE.p <- 
    FC.long %>% 
        filter(gene %in% sal_imm) %>% 
        ggplot(aes(x = Age, y = expression, fill = Feed)) +
          geom_boxplot(colour = "black", outlier.shape = NA, lwd = 0.3) +
          geom_point(position=position_jitterdodge(jitter.width = 0.15),
                     size =0.4) +
          geom_hline(yintercept=1, linewidth=0.3) +
          theme_bw() + 
          facet_wrap(~gene, scales = "free_y", ncol = 3, 
                     labeller = label_parsed) +
          scale_y_continuous(name="Fold-change relative to\nmock-challenged birds", 
                             expand=expansion(c(0,0.3) ) ) +
          scale_shape_manual(values=c(19, 17))+
          scale_fill_manual(values = c("white", "#f5e4ae")) +
          theme(axis.text = element_text(colour = "black", size = text_size), 
                axis.title.y = element_text(colour = "black", size = header_size),
                strip.text.x = element_text(size = header_size),
                legend.position = "right",
                legend.title = element_text(size=header_size),
                legend.text = element_text(size=header_size),
                panel.grid.minor = element_blank(), 
                panel.grid.major = element_blank()
               )


GE_SE.annot.p <- ggdraw(GE_SE.p) + 
                  draw_line(x = c(0.098, 0.12), y = 0.375, 
                            color = "black", size = 0.35) +
                  draw_label(FC.stats %>% 
                                    filter(Age == 22, gene == "'IL-8L2'") %>% 
                                    pull(p_adjust) %>% 
                                    round(.,3), 
                             x = 0.109, y = 0.398, size = 6, color = "black") 


############ DAA ################################################

# read in data
aldex2_raw <- read_tsv(aldex2_path)
ancombc_diet_raw <- read_tsv(ancombc_diet_path)
ancombc_chal_raw <- read_tsv(ancombc_chal_path)


# format data
ancombc_diet_df <- ancombc_diet_raw %>% 
                        mutate("log10(ancombc.diet.padj)" = log10(q_dietgos))

ancombc_chal_df <- 
      ancombc_chal_raw %>% 
                mutate("log10(ancombc.challenge.padj)" = log10(q_challengeSE))

ancombc_df <- left_join(x = ancombc_diet_df, y = ancombc_chal_df) %>% 
                        mutate(OTU = str_extract(taxon, pattern = "Otu\\d{4}"))


aldex2_df <- 
  aldex2_raw %>% 
      mutate("log10(ALDEx2.challenge.padj)" = log10(`challengeSE:pval.holm`),
             "log10(ALDEx2.diet.padj)" = log10(`dietgos:pval.holm`))

df <- left_join(x = ancombc_df,  y = aldex2_df) %>% 
        rowwise() %>% 
        mutate(id = taxon,
               taxon = str_extract(id, pattern = "(?<=\\().+?(?=\\))") %>% str_remove(., "00"),
                chal.max = max(`log10(ancombc.challenge.padj)`,`log10(ALDEx2.challenge.padj)`),
               diet.max =  max(`log10(ancombc.diet.padj)`,`log10(ALDEx2.diet.padj)`),
               chal_fill = case_when(
                                    `challengeSE:pval.holm` >  0.05 | q_challengeSE >  0.05  ~ "p.adj > 0.05",
                                    `challengeSE:pval.holm` < 0.5 && lfc_challengeSE < 0 ~ "Negative LFC",
                                    `challengeSE:pval.holm` < 0.5 && lfc_challengeSE > 0 ~ "Positive LFC"),
                diet_fill = case_when(
                                    `dietgos:pval.holm` >  0.05 | q_dietgos > 0.05 ~ "p.adj > 0.05",
                                    `dietgos:pval.holm` < 0.5 && lfc_dietgos < 0 ~ "Negative LFC",
                                    `dietgos:pval.holm` < 0.5 && lfc_dietgos > 0 ~ "Positive LFC")
               )

# fix levels
df$chal_fill <- factor(df$chal_fill, levels = c("Positive LFC","Negative LFC","p.adj > 0.05"))
df$diet_fill <- factor(df$diet_fill, levels = c("Positive LFC","Negative LFC","p.adj > 0.05"))


# plot results                 
challenge_raw.p <- ggplot(df, aes(x=`log10(ALDEx2.challenge.padj)`, 
                                  y=`log10(ancombc.challenge.padj)`)) + 
  geom_point(aes(color = chal_fill)) + 
  scale_colour_manual(values = c("#D41159","#1A85FF","#808080")) +
  scale_x_reverse() + 
  scale_y_reverse() +
  geom_hline(yintercept=log10(0.05),linetype="dashed",linewidth=0.3) +
  geom_vline(xintercept=log10(0.05),linetype="dashed",linewidth=0.3) +
  theme_bw() +
  theme(aspect.ratio=1,
        plot.margin=margin(t=5.5,r=0, b=5.5,l=11),
        axis.text = element_text(colour = "black", size = text_size),
        axis.title = element_text(colour = "black", size = header_size)) +
  geom_text_repel(aes(label=ifelse(`log10(ancombc.challenge.padj)`< -20,
                             as.character(`taxon`),'')), size=2) +
  labs(x ="log10(ALDEx2.glm FDR)", 
       y = "log10(ANCOM-BC FDR)",
       color = "SE challenge associated FC")


# SE
challenge.p <- df |> filter(OTU != "Otu0010") |> 
  ggplot(aes(x=`log10(ALDEx2.challenge.padj)`, 
             y=`log10(ancombc.challenge.padj)`)) + 
  geom_point(aes(color = chal_fill)) + 
  scale_colour_manual(values = c("#D41159","#1A85FF","#808080"))  +
  scale_x_reverse(limits = c(0, -11.5), expand = c(0,0)) + 
  scale_y_reverse(limits = c(0, -11.5), expand = c(0,0)) +
  geom_hline(yintercept=log10(0.05),linetype="dashed", linewidth=0.3) +
  geom_vline(xintercept=log10(0.05),linetype="dashed", linewidth=0.3) +
  theme_bw() +
  theme(aspect.ratio=1,
        plot.margin=margin(t=5.5,r=0, b=5.5,l=11),
        axis.text = element_text(colour = "black", size = text_size),
        axis.title = element_text(colour = "black", size = header_size)) +
  geom_text_repel(aes(label=ifelse(chal.max < log10(0.05), 
                             as.character(`taxon`),'')), size=2) +
  labs(x ="log10(ALDEx2.glm FDR)", 
       y = "log10(ANCOM-BC FDR)",
       color = "SE challenge associated FC")


# diet
gos.p <- ggplot(df, aes(x=`log10(ALDEx2.diet.padj)`, 
                        y=`log10(ancombc.diet.padj)`)) + 
  geom_point(aes(color = diet_fill)) + 
  scale_colour_manual(values = c("#D41159","#1A85FF", "#808080")) +
  scale_x_reverse(limits = c(0, -8.5), expand = c(0,0)) + 
  scale_y_reverse(limits = c(0, -8.5), expand = c(0,0)) +
  geom_hline(yintercept=log10(0.05),linetype="dashed",linewidth=0.3) +
  geom_vline(xintercept=log10(0.05),linetype="dashed",linewidth=0.3 ) +
  theme_bw() +
  theme(aspect.ratio=1,
        plot.margin=margin(t=5.5,r=0, b=5.5,l=11),
        axis.text = element_text(colour = "black", size = text_size),
        axis.title = element_text(colour = "black", size = header_size)) +
  geom_text_repel(aes(label=ifelse(diet.max < log10(0.05), 
                             as.character(`taxon`),'')),
            size=2)+
  labs(x ="log10(ALDEx2.glm FDR)", 
       y = "log10(ANCOM-BC FDR)",
       color = "jGOS associated FC")


# generate legend
fix_legend <- gos.p + 
                  theme(legend.position="right",
                        legend.title=element_blank(),
                        legend.text = element_text(size=text_size)) +
                  guides(color = guide_legend(override.aes = list(size=4) ))

legend_scat <- get_legend(fix_legend)


DAA.p <- plot_grid(challenge_raw.p + theme(legend.position = "none", 
                                           axis.title = element_text(size = header_size)), 
                   challenge.p + theme(legend.position = "none", 
                                       axis.title = element_text(size =header_size)), 
                   gos.p + theme(legend.position = "none", 
                                 axis.title = element_text(size =header_size)),
                   legend_scat,
                   nrow = 1, labels = c('C','D','E'), 
                   rel_widths = c(0.29,0.29,0.29,0.13))


############ plot ################################################
p <- plot_grid(invsimp_annotate.p,
               pcoa.p,
               DAA.p,  
               GE_SE.annot.p, 
               ncol = 1, 
               labels = c("A","B", "","F"),
               rel_heights = c(0.8, 0.85,0.85,1.4))


# export plots to pdf
pdf(file = figpath, paper = "a4", height = 8)
p
dev.off()


