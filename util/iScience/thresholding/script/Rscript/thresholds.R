library(ggplot2)
library(reshape2)

#######################################################################################################################
# p = 100
#######################################################################################################################

counts <- c(13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31)
heuristic  <- c(1, 1, 2, 2, 1, 2, 1, 3, 3, 2, 2, 0, 0, 1, 1, 2, 0, 0, 2)
simple  <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
indirect  <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
overlap    <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
both  <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

df <- data.frame(counts, heuristic, simple, indirect, overlap, both)
colnames(df) <- c("counts", "Heuristic", "Simple", "Indirect", "Overlap", "Both")
df2 <- melt(df, id = "counts", variable.name="Method")

mini1 <- ggplot(df2, aes(x = counts, y = value, color = Method)) + 
  geom_line(size = 1) + 
  geom_point(size = 3) + 
  labs(title=expression(bold("Threshold vs minimizer counts, p = 100")),
       x = "Number of minimizers", 
       y = "Threshold") + 
  scale_colour_discrete(name = "Method", labels =c("Heuristic", 
                                                   expression("("*bar(italic("I"))*","*bar(italic("O"))*")"), 
                                                   expression("("*italic("I")*","*bar(italic("O"))*")"), 
                                                   expression("("*bar(italic("I"))*","*italic("O")*")"), 
                                                   expression("("*italic("I")*","*italic("O")*")"))) +
  theme_bw() +
  theme(text = element_text(size=30),
        plot.title = element_text(size=30, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        legend.text.align = 0)

#######################################################################################################################
# p = 150
#######################################################################################################################

counts <- c(23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46)
heuristic  <- c(9, 10, 9, 9, 11, 12, 11, 12, 12, 13, 14, 10, 15, 12, 17, 17, 12, 14, 14, 14, 18, 13, 18, 20)
simple  <- c(5, 6, 6, 7, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16)
indirect  <- c(4, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15)
overlap    <- c(5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 17)
both  <- c(4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16)

df <- data.frame(counts, heuristic, simple, indirect, overlap, both)
colnames(df) <- c("counts", "Heuristic", "Simple", "Indirect", "Overlap", "Both")
df2 <- melt(df, id = "counts", variable.name="Method")

mini2 <- ggplot(df2, aes(x = counts, y = value, color = Method)) + 
  geom_line(size = 1) + 
  geom_point(size = 3) + 
  labs(title=expression(bold("Threshold vs minimizer counts, p = 150")),
       x = "Number of minimizers", 
       y = "Threshold") + 
  scale_colour_discrete(name = "Method", labels =c("Heuristic", 
                                                   expression("("*bar(italic("I"))*","*bar(italic("O"))*")"), 
                                                   expression("("*italic("I")*","*bar(italic("O"))*")"), 
                                                   expression("("*bar(italic("I"))*","*italic("O")*")"), 
                                                   expression("("*italic("I")*","*italic("O")*")"))) +
  theme_bw() +
  theme(text = element_text(size=30),
        plot.title = element_text(size=30, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        legend.text.align = 0)

#######################################################################################################################
# p = 250
#######################################################################################################################

counts <- c(44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74)
heuristic  <- c(27, 29, 30, 30, 31, 31, 31, 32, 32, 33, 34, 33, 34, 35, 36, 34, 36, 37, 39, 38, 39, 40, 39, 41, 41, 44, 39, 42, 46, 47, 38)
simple  <- c(25, 26, 26, 27, 28, 29, 29, 30, 31, 31, 32, 33, 33, 34, 35, 36, 36, 37, 38, 38, 39, 40, 40, 41, 42, 43, 43, 44, 45, 45, 46)
indirect  <- c(24, 25, 25, 26, 27, 28, 28, 29, 30, 30, 31, 32, 32, 33, 34, 35, 35, 36, 37, 37, 38, 39, 39, 40, 41, 42, 42, 43, 44, 44, 45)
overlap    <- c(25, 26, 27, 27, 28, 29, 29, 30, 31, 32, 32, 33, 34, 34, 35, 36, 36, 37, 38, 39, 39, 40, 41, 41, 42, 43, 44, 44, 45, 46, 46)
both  <- c(24, 25, 26, 26, 27, 28, 28, 29, 30, 31, 31, 32, 33, 33, 34, 35, 35, 36, 37, 38, 38, 39, 40, 40, 41, 42, 43, 43, 44, 45, 45)

df <- data.frame(counts, heuristic, simple, indirect, overlap, both)
colnames(df) <- c("counts", "Heuristic", "Simple", "Indirect", "Overlap", "Both")
df2 <- melt(df, id = "counts", variable.name="Method")

mini3 <- ggplot(df2, aes(x = counts, y = value, color = Method)) + 
  geom_line(size = 1) + 
  geom_point(size = 3) + 
  labs(title=expression(bold("Threshold vs minimizer counts, p = 250")),
       x = "Number of minimizers", 
       y = "Threshold") + 
  scale_colour_discrete(name = "Method", labels =c("Heuristic", 
                                                   expression("("*bar(italic("I"))*","*bar(italic("O"))*")"), 
                                                   expression("("*italic("I")*","*bar(italic("O"))*")"), 
                                                   expression("("*bar(italic("I"))*","*italic("O")*")"), 
                                                   expression("("*italic("I")*","*italic("O")*")"))) +
  theme_bw() +
  theme(text = element_text(size=30),
        plot.title = element_text(size=30, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        legend.text.align = 0)

#######################################################################################################################
# Save plots
#######################################################################################################################

ggsave(file="<DIR>/mini_100.svg", plot=mini1, width=16, height=9)
ggsave(file="<DIR>/mini_150.svg", plot=mini2, width=16, height=9)
ggsave(file="<DIR>/mini_250.svg", plot=mini3, width=16, height=9)