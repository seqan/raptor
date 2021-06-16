library(ggplot2)
library(reshape2)

tau <- c(0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0)

#######################################################################################################################
# p = 100
#######################################################################################################################

simple <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.99999, 0.99999, 0.99999, 0.99984, 0.99982)
indirect <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.99999, 0.99999, 0.99999, 0.99996, 0.99982, 0.99982)
overlapping <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.99999, 0.99999, 0.99992, 0.99982)
indirect_overlapping <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.99999, 0.99999, 0.99998, 0.99986, 0.99982)

df <- data.frame(tau, simple, indirect, overlapping, indirect_overlapping)
colnames(df) <- c("tau", "Simple", "Indirect", "Overlap", "Both")
df2 <- melt(df, id = "tau", variable.name="Method")

tnr1 <- ggplot(df2, aes(x = tau, y = value, color = Method)) + 
  geom_line(size = 1) + 
  geom_point(size = 3) + 
  ylim(0.999, 1) +
  labs(title=expression(bold("TNR vs "*tau*", p = 100")),
       x =expression(tau), 
       y ="TNR") + 
  scale_colour_discrete(name = "Method", labels =c(expression("("*bar(italic("I"))*","*bar(italic("O"))*")"), 
                                                   expression("("*italic("I")*","*bar(italic("O"))*")"), 
                                                   expression("("*bar(italic("I"))*","*italic("O")*")"), 
                                                   expression("("*italic("I")*","*italic("O")*")"))) +
  theme_bw() +
  theme(text = element_text(size=30),
        plot.title = element_text(size=30, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=40, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        legend.text.align = 0)

#######################################################################################################################
# p = 150
#######################################################################################################################

simple <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.99967)
indirect <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.99969)
overlapping <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.99966)
indirect_overlapping <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.99966)

df <- data.frame(tau, simple, indirect, overlapping, indirect_overlapping)
colnames(df) <- c("tau", "Simple", "Indirect", "Overlap", "Both")
df2 <- melt(df, id = "tau", variable.name="Method")

tnr2 <- ggplot(df2, aes(x = tau, y = value, color = Method)) + 
  geom_line(size = 1) + 
  geom_point(size = 3) + 
  ylim(0.999, 1) +
  labs(title=expression(bold("TNR vs "*tau*", p = 150")),
       x =expression(tau), 
       y ="TNR") + 
  scale_colour_discrete(name = "Method", labels =c(expression("("*bar(italic("I"))*","*bar(italic("O"))*")"), 
                                                   expression("("*italic("I")*","*bar(italic("O"))*")"), 
                                                   expression("("*bar(italic("I"))*","*italic("O")*")"), 
                                                   expression("("*italic("I")*","*italic("O")*")"))) +
  theme_bw() +
  theme(text = element_text(size=30),
        plot.title = element_text(size=30, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=40, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        legend.text.align = 0)

#######################################################################################################################
# p = 250
#######################################################################################################################

simple <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
indirect <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
overlapping <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
indirect_overlapping <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

df <- data.frame(tau, simple, indirect, overlapping, indirect_overlapping)
colnames(df) <- c("tau", "Simple", "Indirect", "Overlap", "Both")
df2 <- melt(df, id = "tau", variable.name="Method")

tnr3 <- ggplot(df2, aes(x = tau, y = value, color = Method)) + 
  geom_line(size = 1) + 
  geom_point(size = 3) + 
  ylim(0.999, 1) +
  labs(title=expression(bold("TNR vs "*tau*", p = 250")),
       x =expression(tau), 
       y ="TNR") + 
  scale_colour_discrete(name = "Method", labels =c(expression("("*bar(italic("I"))*","*bar(italic("O"))*")"), 
                                                   expression("("*italic("I")*","*bar(italic("O"))*")"), 
                                                   expression("("*bar(italic("I"))*","*italic("O")*")"), 
                                                   expression("("*italic("I")*","*italic("O")*")"))) +
  theme_bw() +
  theme(text = element_text(size=30),
        plot.title = element_text(size=30, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=40, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        legend.text.align = 0)

#######################################################################################################################
# Save plots
#######################################################################################################################

ggsave(file="<DIR>/tnr_100.svg", plot=tnr1, width=16, height=9)
ggsave(file="<DIR>/tnr_150.svg", plot=tnr2, width=16, height=9)
ggsave(file="<DIR>/tnr_250.svg", plot=tnr3, width=16, height=9)
