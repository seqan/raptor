library(ggplot2)
library(reshape2)

tau <- c(0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0)

#######################################################################################################################
# p = 100
#######################################################################################################################

simple <- c(0.00073, 0.453, 0.5908, 0.68955, 0.76454, 0.82772, 0.86499, 0.89532, 0.92886, 0.94558, 0.96399, 0.97541, 0.98477, 0.99145, 0.99459, 0.99793, 0.99917, 0.99965, 0.99998, 1.0, 1.0)
indirect <- c(0.00073, 0.49942, 0.65422, 0.76251, 0.80768, 0.86499, 0.89651, 0.92886, 0.95377, 0.96445, 0.97541, 0.98477, 0.99173, 0.99459, 0.99793, 0.99905, 0.99952, 0.99993, 0.99999, 1.0, 1.0)
overlapping <- c(0.00073, 0.40691, 0.54404, 0.63588, 0.72183, 0.79746, 0.8367, 0.87545, 0.91107, 0.92886, 0.95386, 0.96399, 0.97798, 0.98477, 0.99173, 0.99578, 0.99863, 0.99948, 0.99992, 1.0, 1.0)
indirect_overlapping <- c(0.00073, 0.45214, 0.59084, 0.72158, 0.77742, 0.83325, 0.87545, 0.91107, 0.92886, 0.95386, 0.96445, 0.97798, 0.98482, 0.99173, 0.99553, 0.99794, 0.99929, 0.99982, 0.99998, 1.0, 1.0)

df <- data.frame(tau, simple, indirect, overlapping, indirect_overlapping)
colnames(df) <- c("tau", "Simple", "Indirect", "Overlap", "Both")
df2 <- melt(df, id = "tau", variable.name="Method")

tpr1 <- ggplot(df2, aes(x = tau, y = value, color = Method)) + 
  geom_line(size = 1) + 
  geom_point(size = 3) + 
  labs(title=expression(bold("TPR vs "*tau*", p = 100")),
       x =expression(tau), 
       y ="TPR") + 
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

simple <- c(0.00036, 0.29619, 0.41984, 0.52279, 0.60512, 0.66181, 0.72736, 0.77462, 0.81752, 0.86014, 0.88979, 0.91792, 0.93439, 0.95534, 0.97285, 0.98427, 0.99134, 0.99798, 0.99903, 0.99996, 1.0)
indirect <- c(0.00036, 0.34301, 0.49576, 0.60005, 0.6777, 0.74428, 0.78867, 0.82893, 0.86795, 0.90259, 0.93095, 0.95452, 0.97154, 0.9827, 0.98938, 0.99542, 0.99798, 0.99903, 0.99986, 1.0, 1.0)
overlapping <- c(0.00036, 0.25031, 0.36626, 0.46799, 0.55167, 0.62782, 0.69667, 0.74434, 0.7882, 0.82794, 0.86095, 0.8973, 0.92558, 0.9528, 0.96921, 0.9827, 0.99127, 0.99561, 0.99903, 0.99995, 1.0)
indirect_overlapping <- c(0.00036, 0.29841, 0.4423, 0.54899, 0.62782, 0.69962, 0.75958, 0.80379, 0.84256, 0.87809, 0.90844, 0.93434, 0.95601, 0.97271, 0.98423, 0.99134, 0.99565, 0.99892, 0.99983, 0.99998, 1.0)

df <- data.frame(tau, simple, indirect, overlapping, indirect_overlapping)
colnames(df) <- c("tau", "Simple", "Indirect", "Overlap", "Both")
df2 <- melt(df, id = "tau", variable.name="Method")

tpr2 <- ggplot(df2, aes(x = tau, y = value, color = Method)) + 
  geom_line(size = 1) + 
  geom_point(size = 3) + 
  labs(title=expression(bold("TPR vs "*tau*", p = 150")),
       x =expression(tau), 
       y ="TPR") + 
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

simple <- c(0.0001, 0.1566, 0.26191, 0.34652, 0.41111, 0.49367, 0.5342, 0.5998, 0.66763, 0.72825, 0.77037, 0.81721, 0.86526, 0.89641, 0.92671, 0.95544, 0.97447, 0.98821, 0.99603, 0.9994, 1.0)
indirect <- c(0.0001, 0.20377, 0.31854, 0.41111, 0.50114, 0.5671, 0.63298, 0.69975, 0.76045, 0.80187, 0.84307, 0.88704, 0.9218, 0.94837, 0.96607, 0.98064, 0.9893, 0.99582, 0.99882, 0.99987, 1.0)
overlapping <- c(0.0001, 0.12663, 0.2231, 0.30371, 0.37675, 0.43801, 0.50115, 0.5671, 0.63298, 0.66763, 0.73032, 0.78888, 0.83178, 0.87406, 0.90883, 0.94232, 0.96829, 0.98385, 0.99429, 0.99904, 1.0)
indirect_overlapping <- c(0.0001, 0.16935, 0.27954, 0.36991, 0.45064, 0.53154, 0.59917, 0.66544, 0.69975, 0.76046, 0.81666, 0.86589, 0.89641, 0.92671, 0.95542, 0.97388, 0.98558, 0.99402, 0.9982, 0.99978, 1.0)

df <- data.frame(tau, simple, indirect, overlapping, indirect_overlapping)
colnames(df) <- c("tau", "Simple", "Indirect", "Overlap", "Both")
df2 <- melt(df, id = "tau", variable.name="Method")

tpr3 <- ggplot(df2, aes(x = tau, y = value, color = Method)) + 
  geom_line(size = 1) + 
  geom_point(size = 3) + 
  labs(title=expression(bold("TPR vs "*tau*", p = 250")),
       x =expression(tau), 
       y ="TPR") + 
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

ggsave(file="<DIR>/tpr_100.svg", plot=tpr1, width=16, height=9)
ggsave(file="<DIR>/tpr_150.svg", plot=tpr2, width=16, height=9)
ggsave(file="<DIR>/tpr_250.svg", plot=tpr3, width=16, height=9)