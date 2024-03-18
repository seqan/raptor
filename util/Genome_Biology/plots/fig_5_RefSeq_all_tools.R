# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

library(lubridate)
library(here)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

timings <- data.frame(
  c("COBS", "00:01:15", 1),
  c("Bifrost", "01:17:43", 1),
  c("Metagraph", "00:03:46", 1),
  c("SeqOthello", "00:29:34", 1),
  c("Mantis", "00:07:27", 1),
  c("IBF", "00:02:31", 1),
  c("HIBF", "00:01:01", 1),

  c("COBS", "00:11:34", 2),
  c("Bifrost", "01:21:55", 2),
  c("Metagraph", "00:13:24", 2),
  c("SeqOthello", "01:43:48", 2),
  c("Mantis", "00:18:47", 2),
  c("IBF", "00:03:59", 2),
  c("HIBF", "00:01:15", 2),

  c("COBS", "00:57:13", 3),
  c("Bifrost", "01:39:00", 3),
  c("Metagraph", "00:52:09", 3),
  c("SeqOthello", NA, 3),
  c("Mantis", "01:05:47", 3),
  c("IBF", "00:08:28", 3),
  c("HIBF", "00:02:29", 3),

  c("COBS", "01:50:02", 4),
  c("Bifrost", "02:02:33", 4),
  c("Metagraph", "01:37:50", 4),
  c("SeqOthello", NA, 4),
  c("Mantis", "02:18:33", 4),
  c("IBF", "00:13:08", 4),
  c("HIBF", "00:04:11", 4)
)

timings <- data.frame(t(timings))
colnames(timings) <- c("Method", "Time (HH:)MM:SS", "Point")

timings$`Time (HH:)MM:SS` <- period_to_seconds(hms(as.character(timings$`Time (HH:)MM:SS`)))
timings$Point <- as.numeric(timings$Point)

memory <- data.frame(
  c("COBS", 106.02, 1),
  c("Bifrost", 282.89, 1),
  c("Metagraph", 231.42, 1),
  c("SeqOthello", 30.70, 1),
  c("Mantis", 468.55, 1),
  c("IBF", 374.04, 1),
  c("HIBF", 133.05, 1),

  c("COBS", 127.47, 2),
  c("Bifrost",  282.89, 2),
  c("Metagraph", 249.02, 2),
  c("SeqOthello", 287.71, 2),
  c("Mantis", 470.67, 2),
  c("IBF",  374.38, 2),
  c("HIBF", 132.70, 2),

  c("COBS", 222.91, 3),
  c("Bifrost", 282.89, 3),
  c("Metagraph", 249.56, 3),
  c("SeqOthello", NA, 3),
  c("Mantis", 472.93, 3),
  c("IBF", 376.83, 3),
  c("HIBF", 134.57, 3),

  c("COBS",  342.51, 4),
  c("Bifrost", 282.89, 4),
  c("Metagraph", 254.97, 4),
  c("SeqOthello", NA, 4),
  c("Mantis", 496.40, 4),
  c("IBF", 379.83, 4),
  c("HIBF", 136.47, 4)
)

memory <- data.frame(t(memory))
colnames(memory) <- c("Method", "Memory", "Point")
memory$Memory <- as.numeric(as.character(memory$Memory))
memory$Point <- as.numeric(memory$Point)

index.size <- data.frame(
  c("COBS", 104),
  c("Bifrost", 68),
  c("Metagraph", 216),
  c("SeqOthello", 315),
  c("Mantis", 496),
  c("IBF", 374),
  c("HIBF", 133)
)

index.size <- data.frame(t(index.size))
colnames(index.size) <- c("Method", "Size")
index.size$Size <- as.numeric(as.character(index.size$Size))

index.time <- data.frame(
  c("COBS",	"00:34:04"),
  c("Bifrost", "17:46:25")
  c("Metagraph",	"56:36:22"),
  c("SeqOthello",	"49:26:04"),
  c("Mantis", "45:19:00"),
  c("IBF", "00:32:37"),
  c("HIBF", "00:38:28")
)

index.time <- data.frame(t(index.time))
colnames(index.time) <- c("Method", "Time")
index.time$Time <- period_to_seconds(hms(as.character(index.time$Time)))

#               bifrost    COBS       HIBF       IBF        mantis     metagraph    SeqOthello
my.colors <- c("#e29b03", "#3017b0", "darkgreen", "#8fbf5f", "#826546", "#9e4dc9", "#43a5f0")

common_theme<-  theme_classic(base_size=7) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.8/.pt),
        axis.ticks = element_line(colour = "black", linewidth = 1/.pt),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black", face="bold"),
        plot.title = element_text(colour = "black", face="bold"),
        text=element_text(colour="black"))

g1<- ggplot(timings, aes(x = Point, y = `Time (HH:)MM:SS`)) +
  geom_line(aes(color = Method), linewidth=1/.pt) +
  geom_point(aes(color = Method), size =2/.pt) +
  scale_color_manual(values = my.colors) +
  scale_y_continuous(breaks=c(0, 1800, 3600, 5400, 7200), labels=c("0"="0.0", "1800"="0.5", "3600"="1.0", "5400"="1.5", "7200"="2.0")) +
  scale_x_continuous(labels=c("1" = "1K", "2" = "1Mio", "3" = "5Mio", "4" = "10Mio")) +
  ggtitle("a") +
  ylab("Query time in hours") +
  xlab("Number of transcripts") +
  common_theme +
  theme(axis.title.x = element_text(margin = margin(t = -20)))

g2<- ggplot(index.time, aes(x=Method, y=Time, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge2(width = 2/.pt, preserve = "single"), linewidth=0.5/.pt, color="black", show.legend = FALSE) +
  ggtitle("b") +
  ylab("Build time in hours") +
  labs (x = NULL) +
  scale_y_continuous(breaks=c(0, 72000, 144000), labels=c("0"="0", "72000"="20", "144000"="40")) +
  scale_fill_manual(values = my.colors) +
  common_theme +
  theme(axis.text.x = element_text(angle=45, hjust=1))

g3<- ggplot(memory, aes(x = Point, y = Memory)) +
  geom_line(aes(color = Method), linewidth=1/.pt, show.legend = FALSE) +
  geom_point(aes(color = Method), size =2/.pt, show.legend = FALSE) +
  scale_color_manual(values = my.colors) +
  scale_x_continuous(labels=c("1" = "1K", "2" = "1Mio", "3" = "5Mio", "4" = "10Mio")) +
  ggtitle("c") +
  expand_limits(y = 0) +
  ylab("Query RAM in GiB") +
  xlab("Number of transcripts") +
  common_theme

g4<- ggplot(index.size, aes(x=Method, y=Size, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge2(width = 2/.pt, preserve = "single"), color="black", linewidth=0.5/.pt, show.legend = FALSE) +
  theme_classic(base_size=7) +
  ggtitle("d") +
  labs (x = NULL) +
  ylab("Index size in GiB") +
  scale_fill_manual(values = my.colors) +
  common_theme +
  theme(axis.text.x = element_text(angle=45, hjust=1))

layout <- "
AAAB
AAAC
AAAD
"

gg<- g1  +
    g2 + g3 + g4 +
  plot_layout(design = layout, guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(), legend.text=element_text(colour="black", size=7))

pdf(here("fig_5_RefSeq_all_tools_tmp.pdf"), width = 5, height = 4.5)
print(gg)
invisible(dev.off())
