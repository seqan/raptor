library(lubridate)

timings <- data.frame(
  c("COBS", "00:01:15", "1", "1T"),
  c("metagraph", "00:09:30", "1", "1T"),
  c("SeqOthello", "00:29:34", "1", "1T"),
  c("raptor", "00:01:01", "1", "1T"),
  c("mantis", "00:07:27", "1", "1T"),
  c("COBS", "00:11:34", "2", "1M"),
  c("metagraph", "00:55:30", "2", "1M"),
  c("SeqOthello", "01:43:48", "2", "1M"),
  c("mantis", "00:18:47", "2", "1M"),
  c("raptor", "00:01:15", "2", "1M"),
  c("COBS", "00:57:13", "3", "5Mio"),
  c("metagraph", "04:27:43", "3", "5Mio"),
  c("SeqOthello", "99:99:99", "3", "5Mio"),
  c("raptor", "00:02:29", "3", "5Mio"),
  c("mantis", "01:05:47", "3", "5Mio"),
  c("COBS", "01:50:02", "4", "10Mio"),
  c("metagraph", "08:00:44", "4", "10Mio"),
  c("SeqOthello", "99:99:99", "4", "10Mio"),
  c("mantis", "02:18:33", "4", "10Mio"),
  c("raptor", "00:04:11", "4", "10Mio"),
  c("ibf", "00:02:31", "1", "1T"),
  c("ibf", "00:03:59", "2", "1Mio"),
  c("ibf", "00:08:28", "3", "5Mio"),
  c("ibf", "00:13:08", "4", "10Mio"),
  c("bifrost", "01:17:43", "1", "1T"),
  c("bifrost", "01:21:55", "2", "1Mio"),
  c("bifrost", "01:39:00", "3", "5Mio"),
  c("bifrost", "02:02:33", "4", "10Mio")
)

timings <- data.frame(t(timings))
colnames(timings) <- c("Method", "Time (HH:)MM:SS", "Point", "Label")

timings$`Time (HH:)MM:SS` <- period_to_seconds(hms(as.character(timings$`Time (HH:)MM:SS`)))
timings$Point <- as.numeric(timings$Point )

memory <- data.frame(
  c("COBS", 127.47, 2),
  c("metagraph", 336.31, 2),
  c("SeqOthello", 287.71, 2),
  c("mantis", 470.67, 2),
  c("raptor", 132.70, 2),
  c("COBS", 106.02, 1),
  c("metagraph", 321.57, 1),
  c("SeqOthello", 30.70, 1),
  c("raptor", 133.05, 1),
  c("mantis", 468.55, 1),
  c("COBS", 222.91, 3),
  c("metagraph", 334.96, 3),
  c("SeqOthello", NA, 3),
  c("raptor", 134.57, 3),
  c("mantis", 472.93, 3),
  c("SeqOthello", NA, 4),
  c("metagraph", 342.08, 4),
  c("COBS",  342.51, 4),
  c("raptor", 136.47, 4),
  c("mantis", 496.40, 4),
  c("ibf", 374.04, 1),
  c("ibf",  374.38, 2),
  c("ibf", 376.83, 3),
  c("ibf", 379.83, 4),
  c("bifrost", 282.89, 1),
  c("bifrost",  282.89, 2),
  c("bifrost", 282.89, 3),
  c("bifrost", 282.89, 4)
)

memory <- data.frame(t(memory))
colnames(memory) <- c("Method", "Memory", "Point")
memory$Memory <- as.numeric(as.character(memory$Memory))
memory$Point <- as.numeric(memory$Point)

index.size <- data.frame(
  c("metagraph", 301),
  c("mantis", 496),
  c("COBS", 104),
  c("SeqOthello", 315),
  c("raptor", 133),
  c("ibf", 374),
  c("bifrost", 68)
)

index.size <- data.frame(t(index.size))
colnames(index.size) <- c("Method", "Size")
index.size$Size <- as.numeric(as.character(index.size$Size))

index.time <- data.frame(
  c("metagraph",	"11:58:30"),
  c("COBS",	"00:34:04"),
  c("SeqOthello",	"49:26:04"),
  c("mantis", "45:19:00"),
  c("raptor", "00:38:28"),
  c("ibf", "00:32:37"),
  c("bifrost", "17:46:25")
)

index.time <- data.frame(t(index.time))
colnames(index.time) <- c("Method", "Time")
index.time$Time <- period_to_seconds(hms(as.character(index.time$Time)))

library(ggplot2)
library(RColorBrewer)

#yellow  blue green  red purple darkgreen lightblue
#bifrost COBS  IBF  mantis metagraph raptor SeqOthello
my.colors <- c("#e29b03", "#386db9", "#8fbf5f", "#c14040", "#8f5fbf","darkgreen", "#7497c9")

g1<- ggplot(timings, aes(x = Point, y = `Time (HH:)MM:SS`)) +
  geom_line(aes(color = Method), size=1) +
  geom_point(aes(color = Method), size =2) +
  scale_color_manual(values = my.colors) +
  scale_y_continuous(limits = c(0, 30000)) +
  theme_classic() +
  ggtitle("Query Time") +
  theme(axis.line = element_line(colour = "black", size = 0.8),
        axis.ticks = element_line(colour = "black", size = 1))

g1.2<- ggplot(timings, aes(x = Point, y = `Time (HH:)MM:SS`)) +
  geom_line(aes(color = Method), size=1) +
  geom_point(aes(color = Method), size =2) +
  scale_color_manual(values = my.colors) +
  scale_y_continuous(trans='log10') +
  theme_classic() +
  ggtitle("Log(Query Time)") +
  theme(axis.line = element_line(colour = "black", size = 0.8),
        axis.ticks = element_line(colour = "black", size = 1))

g1.3<- ggplot(timings, aes(x = Point, y = `Time (HH:)MM:SS`)) +
  geom_line(aes(color = Method), size=1) +
  geom_point(aes(color = Method), size =2) +
  scale_color_manual(values = my.colors) +
  scale_y_continuous(limits=c(0, 7000)) +
  scale_x_continuous(limits=c(0, 2)) +
  theme_classic() +
  ggtitle("Log(Query Time)") +
  theme(axis.line = element_line(colour = "black", size = 0.8),
        axis.ticks = element_line(colour = "black", size = 1))

g2<- ggplot(memory, aes(x = Point, y = Memory)) +
  geom_line(aes(color = Method), size=1) +
  geom_point(aes(color = Method), size =2) +
  scale_color_manual(values = my.colors) +
  theme_classic() +
  ggtitle("Query RAM") +
  theme(axis.line = element_line(colour = "black", size = 0.8),
        axis.ticks = element_line(colour = "black", size = 1))

g3<- ggplot(index.size, aes(x=Method, y=Size, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge2(width = 2, preserve = "single"), color="black") +
  theme_classic() +
  ggtitle("Final Index Size") +
  scale_fill_manual(values = my.colors) +
  theme(axis.line = element_line(colour = "black", size = 0.8),
        axis.ticks = element_line(colour = "black", size = 1))

g4<- ggplot(index.time, aes(x=Method, y=Time, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge2(width = 2, preserve = "single"), color="black") +
  theme_classic() +
  ggtitle("Build Time") +
  scale_fill_manual(values = my.colors) +
  theme(axis.line = element_line(colour = "black", size = 0.8),
        axis.ticks = element_line(colour = "black", size = 1))

library(patchwork)

layout <- "
AAAB
AAAC
AAAD
"
g1  +
  inset_element(
    g1.2,
    left = 0.1,
    bottom = 0.6,
    right = 0.5,
    top = 0.95,
    align_to = "full"
  ) +
    g4 + g2 + g3 +
  plot_layout(design = layout, guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")

