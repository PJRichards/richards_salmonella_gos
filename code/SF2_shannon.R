#######################################################################
#
# Figure S2. Supplementing feed with GOS in early life affects changes 
# in the cecal microbiota caused by *Salmonella* challenge
#
#######################################################################

# packages
library("readr")
library("FSA")
library("dplyr")
library("tidyr")
library("ggplot2")
library("forcats")
library("purrr")
library("cowplot")

set.seed(1216)

# input
alpha_path <- "results/mothur/sal.trim.contigs.sort.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.opti_mcc.0.03.pick.groups.ave-std.summary"

# output
figpath <- "submission/figS2_shannon.pdf"

# formatting
group_meta_id <- c("Trial", "Feed", "Challenge", "Age", "Replicate")
cohort_levels <- c("ctl unc","ctl sal","gos unc", "gos sal")
text_size = 9

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

# are shannon data normally distributed?
print("# are data normally distributed?")
alpha %>% 
      select(Feed, Challenge, Age, shannon) %>% 
      group_by(Feed, Challenge, Age) %>% 
      summarise(shapiro = shapiro.test(shannon)$p.value)

alpha_nest <- alpha %>% 
                      group_nest(Age) 
# get stats
shannon_kw <- alpha_nest %>% 
                  mutate(kruskal_raw = map(data, ~ kruskal.test(.x$shannon, .x$Cohort)),
                  kruskal = map(kruskal_raw, broom::tidy)) %>% 
                  unnest(kruskal)

shannon_dunn <- alpha_nest %>% 
                    filter(Age == "35 (15)") %>% 
                    mutate(dunntest = map(data, ~.x %>% dunnTest(shannon ~ Cohort, 
                                               data=., method="bh")))


# plot alpha
shannon.p <- alpha %>% 
                    ggplot(aes(x = Cohort, y = shannon, fill = feed)) +
                        geom_boxplot(outlier.shape = NA,lwd=0.3) +
                        geom_point(aes(shape = Challenge), size = 1,
                                       position = position_jitterdodge(0.2)) +
                        scale_y_continuous(limits = c(0, 5), 
                                           breaks = seq(0, 5, by = 1), 
                                           expand = c(0,0)) +
                        scale_x_discrete(labels = c("Mock","SE","Mock","SE")) +
                        facet_grid(cols = vars(Age)) +
                        theme_bw() +
                        theme(legend.position = "None",
                              axis.text.y = element_text(colour = "black"),
                              axis.text.x = element_text(colour = "black"),
                              axis.title = element_text(colour = "black"),
                              panel.grid.minor = element_blank(), 
                              panel.grid.major = element_blank()) +
                        scale_fill_manual(values = c("#f5e4ae","#ffffff")) +
                        xlab("") + ylab("Shannon entropy") 


# print plot
shannon_annotate.p <- ggdraw(shannon.p) +
  #draw_line(x = c(0.112, 0.165), y = 0.12, color = "black", size = 0.35) +
  draw_line(x = c(0.0877, 0.143), y = 0.055, color = "black", size = 0.35) +
  draw_label("ctl", x = 0.11535, y = 0.03, size = text_size) +
  draw_line(x = c(0.195, 0.2503), y = 0.055, color = "black", size = 0.35) +
  draw_label("jGOS", x = 0.22265, y = 0.03, size = text_size) +
  
  draw_line(x = c(0.3235, 0.3788), y = 0.055, color = "black", size = 0.35) +
  draw_label("ctl", x = 0.35115, y = 0.03, size = text_size) +
  draw_line(x = c(0.4308, 0.4861), y = 0.055, color = "black", size = 0.35) +
  draw_label("jGOS", x = 0.45845, y = 0.03, size = text_size) +
  
  draw_line(x = c(0.5593, 0.6146), y = 0.055, color = "black", size = 0.35) +
  draw_label("ctl", x = 0.58695, y = 0.03, size = text_size) +
  draw_line(x = c(0.6666, 0.7219), y = 0.055, color = "black", size = 0.35) +
  draw_label("jGOS", x = 0.69425, y = 0.03, size = text_size) +
  
  draw_line(x = c(0.7951, 0.8504), y = 0.055, color = "black", size = 0.35) +
  draw_label("ctl", x = 0.82275, y = 0.03, size = text_size) +
  draw_line(x = c(0.9024, 0.9577), y = 0.055, color = "black", size = 0.35) +
  draw_label("jGOS", x = 0.93005, y = 0.03, size = text_size) +
  
  # ctl sal - gos sal
  draw_line(x = c(0.8504, 0.9577), y = 0.79, color = "black", size = 0.35) +
  draw_label(round(shannon_dunn[[1,3]][[1]]$res[2,4],3), x = 0.90405, y = 0.81, size = text_size) +
  # gos sal - gos unc
  draw_line(x = c(0.9024, 0.9577), y = 0.72, color = "black", size = 0.35) +
  draw_label(round(shannon_dunn[[1,3]][[1]]$res[6,4],3), x = 0.93005, y = 0.74, size = text_size) +
  # ctl unc - gos sal
  draw_line(x = c(0.7951, 0.9577), y = 0.86, color = "black", size = 0.35) +
  draw_label(round(shannon_dunn[[1,3]][[1]]$res[3,4],3), x = 0.8764, y = 0.88, size = text_size)


# export plot to pdf
pdf(file = figpath, paper = "a4", height = 4)
shannon_annotate.p
dev.off()


