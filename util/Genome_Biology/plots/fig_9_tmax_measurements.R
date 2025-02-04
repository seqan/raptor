# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

library(ggplot2)
library(reshape2)
library(here)

table2 <- data.frame(tmax = as.factor(c(64,128, 192,256, 512,1024,2048,4096,8192)),
                     exp_time  = c(1.00, 0.84, 0.78, 0.83, 1.06, 1.23, 1.51, 2.26, 3.59),
                     exp_mem   = c(1.00, 0.79, 0.68, 0.68, 0.72, 0.80, 0.90, 0.93, 0.84),
                     real_mem = c(1.00, 0.790698, 0.639535, 0.732558, 0.686047, 0.813953, 0.895349, 0.965116, 0.755814),
                     real_time = c(0.647888, 0.682035, 0.440796, 0.486141, 0.43423, 0.530725, 0.728071, 1.18781, 1.71542))

table2 <- data.frame(tmax = as.factor(c(64,128, 192,256, 512,1024,2048,4096,8192)),
                     exp_time  = c(1.00, 0.84, 0.78, 0.83, 1.06, 1.23, 1.51, 2.26, 3.59),
                     exp_mem   = c(1.00, 0.79, 0.68, 0.68, 0.72, 0.80, 0.90, 0.93, 0.84),
                     real_mem = c(1.00, 0.790698, 0.639535, 0.732558, 0.686047, 0.813953, 0.895349, 0.965116, 0.755814),
                     real_time = c(0.756273, 0.644432, 0.469482, 0.437686, 0.444667, 0.54607, 0.713417, 1.11696, 1.80044))

table2$real_time <- table2$real_time / table2$real_time[1] # make runtime relative to 64 HIBF

table2 <- cbind(table2, exp_total = table2$exp_time * table2$exp_mem)
table2 <- cbind(table2, real_total = table2$real_time * table2$real_mem)

cor(table2$exp_total, table2$real_total)
cor(table2$exp_mem, table2$real_mem)
cor(table2$exp_time, table2$real_time)

ttable <- t(table2)

melted <- melt(table2)
melted <- cbind(melted, bench = c(rep("Expected cost", 18), rep("Real cost", 18), rep("Expected cost",9), rep("Real cost",9)))
melted <- cbind(melted, versus = c(rep("Query time",9), rep("Space",18), rep("Query time",9), rep("Space * Query time",18)))

colnames(melted) <- c("tmax", "variable", "Ratio", "bench", "versus")

my.colors <- c( "darkblue", "darkred")

gg <- ggplot(melted, aes(x=tmax, y=Ratio, group=variable)) +
  geom_line(linetype=3, aes(col=bench), linewidth=1/.pt) +
  geom_point(aes(col=bench), size=1/.pt) +
  scale_color_manual(values = my.colors) +
  facet_grid(~ versus) +
  coord_cartesian(ylim = c(0, 4)) +
  theme_bw(base_size=7) +
  xlab(expression(t[m][a][x])) +
  theme(axis.line = element_line(colour = "black", linewidth=1/.pt),
        axis.ticks = element_line(colour = "black", linewidth = 1/.pt),
        axis.text.x = element_text(angle=45, hjust=1, colour="black"),
        axis.text.y = element_text(colour = "black"),
        text=element_text(colour="black"),
        plot.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())

pdf(here("fig_9_tmax_measurements_tmp.pdf"), width = 5, height = 1.8)
print(gg)
invisible(dev.off())
