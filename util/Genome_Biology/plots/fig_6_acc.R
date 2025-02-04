# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

library(lubridate)
library(ggplot2)
library(ggnewscale)
library(here)
library(patchwork)

refseq_data <- read.table(here("fig_6_acc.csv"), sep = ",", quote = "\"", header = T,strip.white = T, row.names = NULL, stringsAsFactors = F)

refseq_data$pre_time <- as.numeric(as.period(ms(refseq_data$pre_time), unit = "sec"))
refseq_data$build_time <- as.numeric(as.period(ms(refseq_data$build_time), unit = "sec"))
refseq_data$query_time <- as.numeric(as.period(ms(refseq_data$query_time), unit = "sec"))
refseq_data$rep_kmers <- factor(refseq_data$rep_kmers, levels = unique(refseq_data$rep_kmers))
refseq_data$k_or_m <- as.factor(refseq_data$k_or_m)

my.colors <- c("#e29b03", "#386db9", "#8fbf5f", "#c14040", "#8f5fbf","darkgreen", "#7497c9")

common_theme<-  theme_classic(base_size=7) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.8/.pt),
        axis.ticks = element_line(colour = "black", linewidth = 1/.pt),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black", face="bold"),
        plot.title = element_text(colour = "black", face="bold"),
        text=element_text(colour="black"))

wk.colors<-c("#3871c1ff","#00adefff", "#37ab9cff", "#2c3792ff", "#8e143dff", "#ef2e32ff", "#f14624ff", "#f78712ff", "#f9a72bff", "#f2ef1dff")
wk.colors<-setNames(wk.colors, levels(refseq_data$rep_kmers))
labels <- levels(refseq_data$rep_kmers)
labels <- setNames(labels, labels)
refseq_data_without_others <- refseq_data[c(1,2,3,4,5,6,7,8,9,10,11,12),]

plot_refseq <- function(my_aes){
  ggplot() +
    geom_bar(data=refseq_data_without_others, my_aes, stat="identity",position=position_dodge2(width = 1.2/.pt, preserve = "single"), linewidth=0.5/.pt, color="black") +
    common_theme +
    scale_fill_manual(values = wk.colors, labels=labels[1:4], breaks = names(wk.colors)[1:4], name = "k",guide = guide_legend(title.position = "top", order = 1, keywidth=1/.pt, keyheight=0.8/.pt)) +
    new_scale_fill() +
    labs(x=NULL) +
    geom_bar(data=refseq_data_without_others, my_aes, stat="identity",position=position_dodge2(width = 1.2/.pt, preserve = "single"), linewidth=0.5/.pt, color="black") +
    scale_fill_manual(values = wk.colors, labels=labels[5:10], breaks = names(wk.colors)[5:10], name = "Minimizers", guide = guide_legend(title.position = "top", order = 0, keywidth=1/.pt, keyheight=0.8/.pt)) +
    theme(legend.position = "right",
          legend.direction = "vertical")
}

g1 <- ggplot(refseq_data, aes(x=index.size, y=Accuracy)) +
  geom_point(aes(shape=k_or_m, color=k_or_m), size = 3/.pt, show.legend=FALSE)+
  geom_text(label=paste(refseq_data$method, refseq_data$rep_kmers), size=1/.pt) +
  ggtitle("a") +
  ylab("Accuracy") +
  xlab("Index size in GiB") +
  common_theme
g2 <- plot_refseq(aes(x=method, y=query_time, fill=rep_kmers)) + ylab("Query time in seconds") + ggtitle("b")
g3 <- plot_refseq(aes(x=method, y=index.size, fill=rep_kmers)) + ylab("Index size in GiB") + ggtitle("c")

layout <- "
AAABC
"

gg<-g1 + g2 + g3 +
  plot_layout(design = layout, guides = "collect") &
  theme(legend.position = "right", legend.direction = "vertical", legend.margin=margin(r=0, l=0), legend.box.margin=margin(r=0, l=0))

pdf(here("fig_6_acc_tmp.pdf"), width = 5, height = 3)
print(gg)
invisible(dev.off())
