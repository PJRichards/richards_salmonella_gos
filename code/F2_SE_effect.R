##############################################################################
#
# Figure 2. GOS-supplemented feed hastens the clearance of *S.* Enteritidis 
# from sub-clinically colonised broiler chickens.
#
##############################################################################

set.seed(1216)

# packages
library("readr")
library("dplyr")
library("cowplot")
library("ggplot2")
library("purrr")
library("tidyr")
library("broom")

# input
SE_counts_path <- "resources/zootechnical/late_SE_counts.csv"
weight_raw_path <- "resources/zootechnical/late_BW_cull.tsv"

# output
figpath <- "submission/fig2_SE_effect.pdf"

############ Salmonella counts #################################################

# read in data
SE_counts_raw <- read_csv(SE_counts_path) 

# count n
print("count n")
SE_counts_raw %>% group_by(Age,Challenge, Diet) %>% count()

# fix feed nomenclature
SE_counts <- SE_counts_raw %>% 
                  mutate(feed = if_else(Diet == "ctl", "ctl",
                                  if_else(Diet == "GOS" & Age < 28,"+GOS","ctl")))

# make placeholder df to harmonize/clarify plot spacing
placeholder_for_count <- tibble(
                                Diet =  c("GOS","ctl","GOS","ctl",
                                          "GOS","ctl","GOS","ctl"), 
                                Age = rep(c(22, 24, 28, 35), each = 2),
                                Challenge = "Mock"
                                ) %>% 
                              mutate(Cohort = paste(Diet, sep = "_", Challenge)) 

# add placeholder count data to add 'Mock' data to facet
# add dpi figure to "age"
SE_counts_placeholder <- bind_rows(SE_counts, placeholder_for_count) %>% 
                              mutate(Age = paste0(Age," (",Age-20,")"))

# format data for stats
# calculate significance
SE_counts_stats <- SE_counts %>% 
                            select(Age, Diet, `Salmonella (log10 CFU)`) %>% 
                            group_nest(Age) %>% 
                            mutate(results = map(data, ~tidy(
                                      wilcox.test(`Salmonella (log10 CFU)` ~ Diet, 
                                                  #paired = FALSE, 
                                                  exact = FALSE,
                                                  data = .)))) %>% 
                            unnest(cols = results)

# plot late counts
counts.p <- 
  SE_counts_placeholder %>% 
      ggplot(aes(x = Cohort, y = `Salmonella (log10 CFU)`, fill = feed)) + 
      geom_boxplot(outlier.shape = NA) + 
      facet_grid(cols = vars(Age)) +
      geom_point(position = position_jitterdodge(0.2)) +
      geom_hline(linetype="dashed", yintercept=c(1, 1.7)) +
      scale_fill_manual(values = c("#f5e4ae","#ffffff")) +
      scale_y_continuous(limits = c(1, 5.5), breaks = seq(1, 5, by = 1)  ) +
      scale_x_discrete(name = "", labels = c("Mock","SE","Mock","SE")) +
      theme_bw() +
      theme(axis.text = element_text(colour = "black"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()
            ) +
      ylab(expression(italic("Salmonella")~textstyle("CFU/g (log"[10]*")"))) +
      facet_grid(cols = vars(Age), scales = "free_x") 


############ body weight #######################################################

# read in mass data
weight_raw <- read_tsv(weight_raw_path) 

# mass n
print("mass n")
weight_raw %>% group_by(Diet, Challenge, Age) %>% summarise(n = n())

# read in male performance objectives 2014
target <- tibble(trial = "target",
                 Diet = "target", 
                 Age = c("22 (2)", "24 (4)", "28 (8)", "35 (15)"), 
                 status = "target", 
                 mean_mass = c(1040, 1209, 1576, 2283),
                 mass_sd = NA, n = NA, cohort = "target")

# add dpi
# fix feed names
weight <- weight_raw  %>% 
            mutate(cohort = paste0(Diet, sep="_", Challenge),
                   Age = paste0(Age," (",Age-20,")"),
                   feed = if_else(Diet == "control", "ctl",
                               if_else(Diet == "GOS" & Age < 28,"+GOS","ctl")))

# get weight stats              
# perform kruskal test for significance
weight_kw <- weight %>% 
              select(-c(Group,Diet,Challenge, feed, bird_ID)) %>% 
              group_by(Age) %>% 
              nest(data = c(LBW, cohort)) %>% 
              mutate(kruskal_raw = map(data, ~ kruskal.test(.x$LBW, .x$cohort)),
                     kruskal = map(kruskal_raw, broom::tidy)) %>% 
              unnest(kruskal)

weight_pw <- 
    weight_kw %>% 
            filter(`p.value` < 0.05) %>% 
            select(-`p.value`) %>%  
            mutate(wilcox_raw = map(data, ~                                         
                                    pairwise.wilcox.test(.x$LBW, g = .x$cohort, 
                                                         p.adjust.method = "BH")),
                   wilcox = map(wilcox_raw, tidy)) %>% 
            unnest(cols = wilcox)

# draw mass plot                
weight.p <-
  weight %>% 
        ggplot(aes(x = cohort, y = LBW, fill = feed)) + 
          geom_boxplot(outlier.shape = NA) +
          geom_point(aes(shape = Challenge), position = position_jitterdodge(0.2)) +
          facet_grid(cols = vars(Age)) +
          scale_fill_manual(values = c("#f5e4ae","#ffffff")) +
          scale_y_continuous(name = "Total Mass (g)", limits = c(800,3200), 
                             breaks = seq(1000, 3000, by = 500), expand = c(0,0)) +
          scale_x_discrete(name = "", labels = c("Mock","SE","Mock","SE")) +
          theme_bw() +
          theme(axis.text = element_text(colour = "black"),
                panel.grid.minor = element_blank(), 
                panel.grid.major = element_blank()) + 
          geom_hline(data = target %>% filter(Age >= 22),  
                     aes(yintercept = mean_mass), linetype = "dashed")
  


        
# Print unannotated figure 
zootech.p <- plot_grid(counts.p + theme(legend.position="none"), 
                       weight.p + theme(legend.position="none"), 
                       ncol = 1, axis = "lr", align = "v", 
                       labels = c("A", "B", "C"))

# grab key stats
# pen mass
weight_SE_GOSvMock_ctl_24 <- weight_pw %>% 
                                filter(Age == "24 (4)", 
                                       group1 == "GOS_SE",
                                       group2 == "control_mock") %>% 
                                pull(`p.value`)

# grab key stats
# counts
SE_counts_stats_28 <- SE_counts_stats %>% filter(Age == 28) %>% pull(p.value)
SE_counts_stats_35 <- SE_counts_stats %>% filter(Age == 35) %>% pull(p.value)

# annotate plot with significance
zootech_annotate.p <- ggdraw(zootech.p) +
  # panel B groups
  draw_line(x = c(0.116, 0.17), y = 0.035, size = 0.43) +
    draw_label("ctl", x=0.143, y= 0.019, size=10) +
  draw_line(x = c(0.22, 0.273), y = 0.035, size = 0.43) +
    draw_label("jGOS", x=0.2465, y=0.019, size=10) +
  
  draw_line(x = c(0.346, 0.398), y = 0.035, size = 0.43) +
    draw_label("ctl", x=0.372, y=0.019, size=10) +
  draw_line(x = c(0.448, 0.501), y = 0.035, size = 0.43) +
    draw_label("jGOS", x=0.4745, y=0.019, size=10) +
  
  draw_line(x = c(0.573, 0.626), y = 0.035, size = 0.43) +
    draw_label("ctl", x=0.5995, y=0.019, size=10) +
  draw_line(x = c(0.677, 0.73), y = 0.035, size = 0.43) +
    draw_label("jGOS", x=0.7035, y=0.019, size=10) +
  
  draw_line(x = c(0.802, 0.855), y = 0.035, size = 0.43) +
    draw_label("ctl", x=0.8285, y=0.019, size=10) +
  draw_line(x = c(0.905, 0.958), y = 0.035, size = 0.43) +
    draw_label("jGOS", x=0.9315, y=0.019, size=10) +
  
  # panel A groups
  draw_line(x = c(0.116, 0.17), y = 0.535, size = 0.43) +
    draw_label("ctl", x=0.143, y=0.519, size=10) +
  draw_line(x = c(0.22, 0.273), y = 0.535, size = 0.43) +
    draw_label("jGOS", x=0.2465, y=0.519, size=10) +
  
  draw_line(x = c(0.346, 0.398), y = 0.535, size = 0.43) +
    draw_label("ctl", x=0.372, y=0.519, size=10) +
  draw_line(x = c(0.448, 0.501), y = 0.535, size = 0.43) +
    draw_label("jGOS", x=0.4757, y=0.519, size=10) +
  
  draw_line(x = c(0.573, 0.626), y = 0.535, size = 0.43) +
    draw_label("ctl", x=0.5995, y=0.519, size=10) +
  draw_line(x = c(0.677, 0.73), y = 0.535, size = 0.43) +
    draw_label("jGOS", x=0.7035, y=0.519, size=10) +
  
  draw_line(x = c(0.802, 0.855), y = 0.535, size = 0.43) +
    draw_label("ctl", x=0.8285, y=0.519, size=10) +
  draw_line(x = c(0.905, 0.958), y = 0.535, size = 0.43) +
    draw_label("jGOS", x=0.9315, y=0.519, size=10) +
  
  # panel B stats
  draw_line(x = c(0.346, 0.501), y = 0.35, size = 0.43) +
    draw_label(label = round(weight_SE_GOSvMock_ctl_24,3), x=0.4235, y=0.365, size=9) +

  # panel A stats
  draw_line(x = c(0.626, 0.73), y = 0.9, size = 0.43) +
    draw_label(label = round(SE_counts_stats_28,3), x=0.678, y=0.915, size=9) +
  
  draw_line(x = c(0.855, 0.958), y = 0.75, size = 0.43) +
    draw_label(label = round(SE_counts_stats_35,3), x=0.9065, y=0.765, size=9)
  
  
# print figure
pdf(file = figpath, paper = "a4", height = 6)
zootech_annotate.p 
dev.off()

