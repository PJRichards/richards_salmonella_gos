################################################################################
#
# Figure S3. ELISA detection of S. Enteritidis specific serum antibodies from 
# chickens infected with S. Enteritidis or mock-infected and raised on either a 
# control or GOS-supplemented diet
#
################################################################################

# packages
library("ggplot2")
library("dplyr")

set.seed(1216)

# in path
abs_path <- "resources/zootechnical/late_serum_Ig.csv"

# output
figpath <- "submission/figS3_serum.pdf"

# read in data
abs_raw <- read.csv(file = abs_path, header = TRUE) 
days <- c(20,22,24,28,35)
labeldays <- c("20","22 (2)","24 (4)","28 (8)","35 (15)")

# format data
# add Cohort grouping variable
abs <- abs_raw %>% 
                mutate(Cohort = paste(Diet, "x", Challenge)) %>% 
                rename(Absorbance = abs)

# plot data
p.raw <- abs %>% 
          ggplot(aes(x=Age, y=Absorbance, group=Cohort, color=Cohort)) +
                geom_line(aes(linetype=Cohort)) +
                geom_point(aes(shape=Cohort)) +
                geom_errorbar(aes(ymin=Absorbance-SEM, ymax=Absorbance+SEM, 
                                  linetype=Cohort),
                              width=0.2,
                              position=position_dodge(width=0.2)) +
                scale_linetype_manual(values = c(1,1,2,2)) +
                scale_shape_manual(values = c(16,16,17,17)) +
                scale_color_manual(values = rep(c("#000000","#995c5c"),2)) +
                scale_x_continuous(breaks = days, labels = labeldays) +
                facet_wrap(vars(antibody), scales = "free_y", nrow = 2) +
                geom_vline(xintercept = 20)

# add theme formatting to plot
p <- p.raw +
      theme_bw() +
        theme(axis.text.y = element_text(colour = "black"),
              axis.text.x = element_text(colour = "black"),
              axis.title = element_text(colour = "black"),
              panel.grid.minor = element_blank(), 
              panel.grid.major = element_blank()
              )

# export plots to pdf
pdf(file = figpath, paper = "a4")
p
dev.off()

