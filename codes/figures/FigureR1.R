#### load packages ####
targetPackages <- c('tidyverse','arrow')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### Figure R1 (Appendix Figure 1 in Peer Review File) ####
##### load dataset #####
df_spider_single <- read_parquet("../data/6_spiderACI/Experiment1_single/parquet/df_spider_single.parquet")

##### make plot #####
g_hist_spider <- 
  ggplot(df_spider_single, aes(x = Spider_speed)) +
  geom_histogram(binwidth = 1) +
  coord_cartesian(xlim = c(0, 50)) +
  xlab("Moving speed of spiders (mm/s)") +
  ylab("Count") +
  theme_bw()

ggsave("../figures/FigureR1.pdf", g_hist_spider, w = 3, h = 3)
