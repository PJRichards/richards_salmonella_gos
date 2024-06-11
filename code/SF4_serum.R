################################################################################
#
# Figure S4. ELISA detection of S. Enteritidis specific serum antibodies from 
# chickens infected with S. Enteritidis or mock-infected and raised on either a 
# control or GOS-supplemented diet
#
################################################################################

# packages
library("ggplot2")
library("dplyr")
library("cowplot")

set.seed(1216)

# p-values
IgY_28_ctl <- 0.049
IgY_28_jGOS <- 0.042
IgY_35_ctl <- 0.024
IgY_35_jGOS <- 0.004

IgA_35_ctl <- 0.008
IgA_35_jGOS <- 0.006

# in path
abs_path <- "resources/zootechnical/late_serum_Ig.csv"

# output
figpath <- "submission/figS4_serum.pdf"

# read in data
abs_raw <- read.csv(file = abs_path, header = TRUE) 

days <- c(20,22,24,28,35)
labeldays <- c("20","22 (2)","24 (4)","28 (8)","35 (15)")
line_size <- 0.35
header_size <- 8

# format data
# add Cohort grouping variable
abs <- abs_raw %>% 
                mutate(Cohort = paste(Diet, "x", Challenge),
                       feed = if_else(Diet == "ctl", "ctl",
                                      if_else(Diet == "jGOS" & Age < 28,"jGOS","ctl")),
                       Age = paste0(Age," (",as.numeric(Age)-20,")")
                       ) %>% 
                rename(Absorbance = abs)

# subset data
iga <- abs %>% filter(antibody=="IgA")
igy <- abs %>% filter(antibody=="IgY")

# add theme formatting to plot
plot_theme <- list(facet_wrap(facets = vars(Age), nrow = 1),
                   geom_errorbar(aes(ymin=Absorbance-SEM, ymax=Absorbance+SEM), 
                                 width=0.1),
                   scale_x_discrete(labels = c("Mock","SE","Mock","SE")),
                   scale_fill_manual(values = c("#ffffff","#f5e4ae")),
                   theme_bw(),
                   theme(legend.position = "None",
                         axis.text.y = element_text(colour = "black"),
                         axis.text.x = element_text(colour = "black"),
                         axis.title = element_text(colour = "black"),
                         panel.grid.minor = element_blank(), 
                         panel.grid.major = element_blank()),
                   xlab("")
                   )

# plot data
iga.p <- iga %>% 
              ggplot(aes(fill=feed, y=Absorbance, x=Cohort)) + 
                    geom_bar(colour="black", position="dodge", stat="identity") + 
                    scale_y_continuous(name = "Absorbance (450 nm)", 
                                       limits = c(0,2), expand = c(0,0)) + 
                     
                    plot_theme

igy.p <- igy %>% 
              ggplot(aes(fill=feed, y=Absorbance, x=Cohort)) + 
                    geom_bar(colour="black", position="dodge", stat="identity") + 
                    scale_y_continuous(name = "Absorbance (450 nm)",
                    limits = c(0,0.8), expand = c(0,0)) + 
                plot_theme

# format figure                
p <- plot_grid(igy.p,
               iga.p,
               ncol = 1, labels = c("A","B"))

annotate.p <- ggdraw(p) +
                    # A. IgY
                    draw_line(x = c(0.103, 0.155), y = 0.539, size = line_size) +
                      draw_label("ctl", x=0.129, y= 0.526, size = header_size) +
                    draw_line(x = c(0.2075, 0.261), y = 0.539, size = line_size) +
                      draw_label("jGOS", x=0.23425, y=0.526, size = header_size) +
                    
                    draw_line(x = c(0.3345, 0.3875), y = 0.539, size = line_size) +
                      draw_label("ctl", x=0.361, y= 0.526, size = header_size) +
                    draw_line(x = c(0.44, 0.4935), y = 0.539, size = line_size) +
                      draw_label("jGOS", x=0.46675, y=0.526, size = header_size) +
    
                    draw_line(x = c(0.567, 0.62), y = 0.539, size = line_size) +
                      draw_label("ctl", x=0.5935, y= 0.526, size = header_size) +
                    draw_line(x = c(0.672, 0.725), y = 0.539, size = line_size) +
                      draw_label("jGOS", x=0.6985, y=0.526, size = header_size) +
  
                    draw_line(x = c(0.7985, 0.8515), y = 0.539, size = line_size) +
                      draw_label("ctl", x=0.825, y= 0.526, size = header_size) +
                    draw_line(x = c(0.9035, 0.9565), y = 0.539, size = line_size) +
                      draw_label("jGOS", x=0.93, y=0.526, size = header_size) +
                    
                    ## stat labels
                    # IgY_28_ctl
                    draw_line(x = c(0.567, 0.62), y = 0.9, size = line_size) +
                      draw_label(IgY_28_ctl, x=0.5935, y= 0.913, size = header_size) +
                    
                    # IgY_28_jGOS
                    draw_line(x = c(0.672, 0.725), y = 0.9, size = line_size) +
                      draw_label(IgY_28_jGOS, x=0.6985, y= 0.913, size = header_size) +
  
                    # IgY_35_ctl
                    draw_line(x = c(0.7985, 0.8515), y = 0.9, size = line_size) +
                      draw_label(IgY_35_ctl, x=0.825, y= 0.913, size = header_size) +
                    
                    # IgY_35_jGOS
                    draw_line(x = c(0.9035, 0.9565), y = 0.9, size = line_size) +
                      draw_label(IgY_35_jGOS, x=0.93, y=0.913, size = header_size) +

                    # B
                    draw_line(x = c(0.103, 0.155), y = 0.04, size = line_size) +
                      draw_label("ctl", x=0.129, y= 0.027, size = header_size) +
                    draw_line(x = c(0.2075, 0.261), y = 0.04, size = line_size) +
                      draw_label("jGOS", x=0.23425, y=0.027, size = header_size) +
        
                    draw_line(x = c(0.3345, 0.3875), y = 0.04, size = line_size) +
                      draw_label("ctl", x=0.361, y= 0.027, size = header_size) +
                    draw_line(x = c(0.44, 0.4935), y = 0.04, size = line_size) +
                      draw_label("jGOS", x=0.46675, y=0.027, size = header_size) +
  
                    draw_line(x = c(0.567, 0.62), y = 0.04, size = line_size) +
                      draw_label("ctl", x=0.5935, y= 0.027, size = header_size) +
                    draw_line(x = c(0.672, 0.725), y = 0.04, size = line_size) +
                      draw_label("jGOS", x=0.6985, y=0.027, size = header_size) +
  
                    draw_line(x = c(0.7985, 0.8515), y = 0.04, size = line_size) +
                      draw_label("ctl", x=0.825, y= 0.027, size = header_size) +
                    draw_line(x = c(0.9035, 0.9565), y = 0.04, size = line_size) +
                      draw_label("jGOS", x=0.93, y=0.027, size = header_size) +
    
                    ## stat labels
                    # IgA_35_ctl 
                    draw_line(x = c(0.7985, 0.8515), y = 0.4, size = line_size) +
                      draw_label(IgA_35_ctl, x=0.825, y= 0.413, size = header_size) +
                      
                    # IgA_35_jGOS 
                    draw_line(x = c(0.9035, 0.9565), y = 0.4, size = line_size) +
                      draw_label(IgA_35_jGOS, x=0.93, y=0.413, size = header_size)
  
# export plots to pdf
pdf(file = figpath, paper = "a4", height = 6)
annotate.p
dev.off()

