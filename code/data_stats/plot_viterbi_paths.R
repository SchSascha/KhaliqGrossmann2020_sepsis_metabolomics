library(ggplot2)
library(data.table)
library(Hmisc)

#Import central functions
source("../function_definitions.R")

#Set path
out_dir <- "../../results/data_stats/"

#Read viterbi paths from file
viterbi_paths <- fread(input = '../HMM/tram_viterbi_paths.csv', sep = "\t", header = TRUE, data.table = FALSE)

#Transform paths to long format table
viterbi_paths_long <- melt(viterbi_paths, id.vars = c("Status", "State", "Case"))

#Plot viterbi paths
viterbi_p <- ggplot(viterbi_paths_long, aes(x = State, y = value, group = Status, color = Status)) +
  facet_wrap(facets = ~ variable, ncol = 5, nrow = ceiling(length(unique(viterbi_paths_long$variable))/5)) +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", geom = "line") +
  stat_summary(fun.ymin = "min", fun.ymax = "max", fun.y = "mean", size = 0.5, geom = "errorbar") +
  #stat_summary(fun.data = "mean_sdl", geom = "errorbar") +
  theme_bw()
plot(viterbi_p)
