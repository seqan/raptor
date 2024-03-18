# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

library(ggplot2)
library(RColorBrewer)
library(here)
library(patchwork)

data <- read.table(here("fig_4_many_user_bins.tsv"), header = T, sep = "\t")

data$search_time <- data$search_time / 3600
data$index_size <- data$index_size / (1024 * 1024 * 1024)
data$build_time <- data$build_time / 60

my.colors <- c( "darkgreen", "#8fbf5f", "orange")

layout <- "
AAABBCC
AAABBCC
AAADDEE
###DDEE
"

theme <- theme_classic(base_size=7) +
         theme(axis.line = element_line(colour = "black", linewidth=1/.pt),
               axis.ticks = element_line(colour = "black", linewidth = 1/.pt),
               axis.text.x = element_text(angle=45, hjust=1, colour="black"),
               axis.text.y = element_text(colour = "black"),
               text=element_text(colour="black"),
               plot.title = element_text(face = "bold"),
               axis.title.y = element_text(face="bold"))

xscale <- scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11), labels = data$Number.of.bins[1:11])
xscale2 <- scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11), labels=c("1"="1024", "2"="", "3"="4096", "4"="", "5"="16384", "6"="", "7"="65536", "8"="", "9"="262144", "10"="", "11"="1048576"))

gline <- geom_line(aes(color = Who), linewidth=1/.pt)
gpoint <- geom_point(aes(color = Who), size =1/.pt)
xlabbins <- xlab("Number of user bins")

g1 <- ggplot(data, aes(x = count, y = search_time)) +
  gline + gpoint + xlabbins + theme + xscale +
  ylab("Query time in minutes") +
  scale_color_manual(values = my.colors) +
  labs(title = "a")

g2 <- ggplot(data, aes(x = count, y = search_ram)) +
  gline + gpoint + xlabbins + theme + xscale2 +
  ylab("Query RAM in GiB") +
  scale_color_manual(values = my.colors) +
  scale_y_continuous(limits = c(0, 16)) +
  labs(title = "b")

g3 <- ggplot(data, aes(x = count, y = index_size)) +
  gline + gpoint + xlabbins + theme + xscale2 +
  ylab("Index size in GiB") +
  scale_color_manual(values = my.colors) +
  scale_y_continuous(limits = c(0, 16)) +
  labs(title = "c")

g4 <- ggplot(data, aes(x = count, y = build_time)) +
  gline + gpoint + xlabbins + theme + xscale2 +
  ylab("Build time in minutes") +
  scale_color_manual(values = my.colors) +
  scale_y_log10() +
  labs(title = "d")

g5 <- ggplot(data, aes(x = count, y = build_ram)) +
  gline + gpoint + xlabbins + theme + xscale2 +
  ylab("Build RAM in GiB") +
  scale_color_manual(values = my.colors) +
  labs(title = "e")

gg <- g1 + g2 + g3 + g4 + g5 + plot_layout(design = layout, guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(), legend.text=element_text(colour="black", size=7))

pdf(here("fig_4_many_user_bins_tmp.pdf"), width = 5, height = 3.5)
print(gg)
invisible(dev.off())

