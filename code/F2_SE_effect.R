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
library("tidyr")
library("cowplot")
library("ggplot2")
library("purrr")
library("broom")

# input
SE_counts_path <- "resources/zootechnical/late_SE_counts.csv"
vilus_raw_path <- "resources/zootechnical/late_vilus_length.csv"
mass_raw_path <- "resources/zootechnical/late_BW_cull.tsv"

# output
figpath <- "submission/fig2_SE_effect.pdf"


############ Salmonella counts ###############################################

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
placeholder <- tibble(
                      Diet =  c("GOS","ctl","GOS","ctl","GOS","ctl","GOS","ctl"), 
                      Age = rep(c(22, 24, 28, 35), each = 2),
                      Challenge = "Mock"
                      ) %>% 
                mutate(Cohort = paste(Diet, sep = "_", Challenge)) 

# add placeholder count data to add 'Mock' data to facet
# add dpi figure to "age"
SE_counts_placeholder <- bind_rows(SE_counts, placeholder) %>% 
                              mutate(Age = paste0(Age," (",Age-20,")"))

# format data for stats
# calculate significance
SE_counts_stats <- SE_counts %>% 
                            select(Age, Diet, `Salmonella (log10 CFU)`) %>% 
                            group_nest(Age) %>% 
                            mutate(results = map(data, ~tidy(
                                      wilcox.test(`Salmonella (log10 CFU)` ~ Diet, 
                                                  paired = FALSE, 
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
      ylab(expression(italic("Salmonella")~textstyle("CFU (log"[10]*")"))) +
      facet_grid(cols = vars(Age), scales = "free_x") 


############ Growth performance ###############################################

# read in mass data
mass_raw <- read_tsv(mass_raw_path) 

# mass n
print("mass n")
mass_raw %>% group_by(Diet, Challenge, Age) %>% summarise(n = n())

# read in aviagen ross 308 male performance objectives 2014
target <- tibble(trial = "target",
                 Diet = "target", 
                 Age = c("22 (2)", "24 (4)", "28 (8)", "35 (15)"), 
                 status = "target", 
                 mean_mass = c(1040, 1209, 1576, 2283),
                 mass_sd = NA, n = NA, cohort = "target")

# add dpi
# fix feed names
mass <- mass_raw  %>% 
  mutate(cohort = paste0(Diet, sep="_", Challenge),
         Age = paste0(Age," (",Age-20,")"),
         feed = if_else(Diet == "control", "ctl",
                        if_else(Diet == "GOS" & Age < 28,"+GOS","ctl")))

# get mass stats              
mass_TukeyHSD <- mass %>% 
                    group_by(Age) %>% 
                    nest() %>%  
                    mutate(Tukey = map(data, 
                                       ~ TukeyHSD(aov(LBW ~ Challenge * Diet, 
                                                      data=.x)) %>% 
                                                      tidy())) %>% 
                    select(-data) %>% 
                    unnest(cols = Tukey) %>% 
                    filter(term == "Challenge:Diet")

# draw mass plot                
mass.p <-
  mass %>% 
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
  

############ Villus height #####################################################

# read in histo data
villus_len_raw <- read_csv(vilus_raw_path) 

# format histo data
villus_len <- villus_len_raw %>% 
                    group_by(pen, Diet, Challenge, Age, Cohort) %>%
                    summarise("VL" = mean(VL),
                              n = n()) %>% 
                    ungroup()
                   
# count villus n   
print("villus n")
villus_len %>% 
            group_by(Age, Diet, Challenge) %>% 
            summarise(n = n()) 

# add dpi
# fix feed names
villus_len_fmt <- villus_len %>% 
                    mutate(Age = paste0(Age," (",Age-20,")"),
                            feed = if_else(Diet == "ctl", "ctl",
                              if_else(Diet == "GOS" & Age < 28,"+GOS","ctl")))
                        
# get mass stats              
villus_TukeyHSD <- 
          villus_len_fmt %>% 
              group_by(Age) %>% 
              nest() %>%  
              mutate(Tukey = map(data, ~ TukeyHSD(aov(VL ~ Challenge * Diet, data=.x)) %>% 
                       tidy())) %>% 
              select(-data) %>% 
              unnest(cols = Tukey) %>% 
              filter(term == "Challenge:Diet")

# plot vilus height
villus_len.p <- 
  villus_len_fmt  %>% 
          ggplot(aes(x = Cohort, y = VL, fill = feed)) + 
            geom_boxplot(outlier.shape = NA) +
            geom_point(aes(shape = Challenge), position = position_jitterdodge(0.2)) +
            facet_grid(cols = vars(Age)) +
            scale_fill_manual(values = c("#f5e4ae","#ffffff")) +
            scale_x_discrete(name = "", labels = c("Mock","SE","Mock","SE")) +
            theme_bw() +
            theme(axis.text = element_text(colour = "black"),
                  panel.grid.minor = element_blank(), 
                  panel.grid.major = element_blank()) +
          ylab(label = "Villus height (\u00b5m)")
        
# Print unannotated figure 
zootech.p <- plot_grid(counts.p + theme(legend.position="none"), 
                       mass.p + theme(legend.position="none"), 
                       villus_len.p + theme(legend.position="none"),
                       ncol = 1, axis = "lr", align = "v", labels = c("A", "B", "C"))

# grab key stats
# pen mass
mass_SE_GOSvMock_ctl_22 <- mass_TukeyHSD %>% 
                                filter(Age == "22 (2)", 
                                       contrast == "SE:GOS-mock:control") %>% 
                                pull(`adj.p.value`)

mass_SE_GOSvMock_ctl_24 <- mass_TukeyHSD %>% 
                                filter(Age == "24 (4)", 
                                       contrast == "SE:GOS-mock:control") %>% 
                                pull(`adj.p.value`)

# grab key stats
# counts
SE_counts_stats_28 <- SE_counts_stats %>% filter(Age == 28) %>% pull(p.value)
SE_counts_stats_35 <- SE_counts_stats %>% filter(Age == 35) %>% pull(p.value)

# annotate plot with significance
zootech_annotate.p <- ggdraw(zootech.p) +
  # panel C groups
  draw_line(x = c(0.117, 0.17), y = 0.035, size = 0.43) +
    draw_label("ctl", x=0.1435, y= 0.019, size=10) +
  draw_line(x = c(0.218, 0.271), y = 0.035, size = 0.43) +
    draw_label("jGOS", x=0.2445, y=0.019, size=10) +
  
  draw_line(x = c(0.346, 0.399), y = 0.035, size = 0.43) +
    draw_label("ctl", x=0.3725, y=0.019, size=10) +
  draw_line(x = c(0.448, 0.501), y = 0.035, size = 0.43) +
    draw_label("jGOS", x=0.4757, y=0.019, size=10) +
  
  draw_line(x = c(0.57, 0.628), y = 0.035, size = 0.43) +
    draw_label("ctl", x=0.5965, y=0.019, size=10) +
  draw_line(x = c(0.676, 0.729), y = 0.035, size = 0.43) +
    draw_label("jGOS", x=0.7025, y=0.019, size=10) +
  
  draw_line(x = c(0.804, 0.857), y = 0.035, size = 0.43) +
    draw_label("ctl", x=0.8305, y=0.019, size=10) +
  draw_line(x = c(0.905, 0.958), y = 0.035, size = 0.43) +
    draw_label("jGOS", x=0.9315, y=0.019, size=10) +
  
  # panel B groups
  draw_line(x = c(0.117, 0.17), y = 0.363, size = 0.43) +
    draw_label("ctl", x=0.1435, y=0.347, size=10) +
  draw_line(x = c(0.218, 0.271), y = 0.363, size = 0.43) +
    draw_label("jGOS", x=0.2445, y=0.347, size=10) +
  
  draw_line(x = c(0.346, 0.399), y = 0.363, size = 0.43) +
    draw_label("ctl", x=0.3725, y=0.347, size=10) +
  draw_line(x = c(0.448, 0.501), y = 0.363, size = 0.43) +
    draw_label("jGOS", x=0.4757, y=0.347, size=10) +
  
  draw_line(x = c(0.57, 0.628), y = 0.363, size = 0.43) +
    draw_label("ctl", x=0.5965, y=0.347, size=10) +
  draw_line(x = c(0.676, 0.729), y = 0.363, size = 0.43) +
    draw_label("jGOS", x=0.7025, y=0.347, size=10) +
  
  draw_line(x = c(0.804, 0.857), y = 0.363, size = 0.43) +
    draw_label("ctl", x=0.8305, y=0.347, size=10) +
  draw_line(x = c(0.905, 0.958), y = 0.363, size = 0.43) +
    draw_label("jGOS", x=0.9315, y=0.347, size=10) +
  
  # panel A groups
  draw_line(x = c(0.117, 0.17), y = 0.7, size = 0.43) +
    draw_label("ctl", x=0.1435, y=0.684, size=10) +
  draw_line(x = c(0.218, 0.271), y = 0.7, size = 0.43) +
    draw_label("jGOS", x=0.2445, y=0.684, size=10) +
  
  draw_line(x = c(0.346, 0.399), y = 0.7, size = 0.43) +
    draw_label("ctl", x=0.3725, y=0.684, size=10) +
  draw_line(x = c(0.448, 0.501), y = 0.7, size = 0.43) +
    draw_label("jGOS", x=0.4757, y=0.684, size=10) +
  
  draw_line(x = c(0.57, 0.628), y = 0.7, size = 0.43) +
    draw_label("ctl", x=0.5965, y=0.684, size=10) +
  draw_line(x = c(0.676, 0.729), y = 0.7, size = 0.43) +
    draw_label("jGOS", x=0.7025, y=0.684, size=10) +
  
  draw_line(x = c(0.804, 0.857), y = 0.7, size = 0.43) +
    draw_label("ctl", x=0.8305, y=0.684, size=10) +
  draw_line(x = c(0.905, 0.958), y = 0.7, size = 0.43) +
    draw_label("jGOS", x=0.9315, y=0.684, size=10) +

  # panel B stats
  draw_line(x = c(0.117, 0.271), y = 0.58, size = 0.43) +
    draw_label(label = round(mass_SE_GOSvMock_ctl_22,3), x=0.194, y=0.595, size=9) +
  draw_line(x = c(0.346, 0.501), y = 0.58, size = 0.43) +
    draw_label(label = round(mass_SE_GOSvMock_ctl_24,3), x=0.4235, y=0.595, size=9) +

  # panel A stats
  draw_line(x = c(0.628, 0.729), y = 0.925, size = 0.43) +
    draw_label(label = round(SE_counts_stats_28,3), x=0.6785, y=0.94, size=9) +
  
  draw_line(x = c(0.857, 0.958), y = 0.8545, size = 0.43) +
    draw_label(label = round(SE_counts_stats_35,3), x=0.9075, y=0.8695, size=9)
  
  
# print figure
pdf(file = figpath, paper = "a4")
zootech_annotate.p 
dev.off()

