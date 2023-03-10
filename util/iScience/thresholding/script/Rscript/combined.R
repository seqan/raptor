library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

grid_arrange_shared_legend <-
  function(...,
           ncol = length(list(...)),
           nrow = 1,
           position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(text = element_text(size=20), legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }

source("<DIR>/tpr.R")
source("<DIR>/tnr.R")
source("<DIR>/fscore.R")

tpr1 <- tpr1 + geom_line(size = 0.25) + labs(title="p = 100", x="", y = "TPR") + theme(text = element_text(size=10), legend.position = "none",  axis.title.x=element_blank(), axis.title.y = element_text(size=20, face="bold"))
tpr2 <- tpr2 + geom_line(size = 0.25) + labs(title="p = 150", x="", y = "") + theme(text = element_text(size=10), legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank())
tpr3 <- tpr3 + geom_line(size = 0.25) + labs(title="p = 250", x="", y = "") + theme(text = element_text(size=10), legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank())
tnr1 <- tnr1 + geom_line(size = 0.25) + labs(title="", x="", y = "TNR") + theme(text = element_text(size=10), legend.position = "none",  axis.title.x=element_blank(), axis.title.y = element_text(size=20, face="bold"))
tnr2 <- tnr2 + geom_line(size = 0.25) + labs(title="", x="", y = "") + theme(text = element_text(size=10), legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank())
tnr3 <- tnr3 + geom_line(size = 0.25) + labs(title="", x="", y = "") + theme(text = element_text(size=10), legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank())
f1 <- f1 + geom_line(size = 0.25) + labs(title="", x="", y = "F-score") + theme(text = element_text(size=10), legend.position = "none",  axis.title.x=element_blank(), axis.title.y = element_text(size=20, face="bold"))
f2 <- f1 + geom_line(size = 0.25) + labs(title="", x="", y = "") + theme(text = element_text(size=10), legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank())
f3 <- f1 + geom_line(size = 0.25) + labs(title="", x="", y = "") + theme(text = element_text(size=10), legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank())
tpr1 <- delete_layers(tpr1, "GeomPoint")
tpr2 <- delete_layers(tpr2, "GeomPoint")
tpr3 <- delete_layers(tpr3, "GeomPoint")
tnr1 <- delete_layers(tnr1, "GeomPoint")
tnr2 <- delete_layers(tnr2, "GeomPoint")
tnr3 <- delete_layers(tnr3, "GeomPoint")
f1 <- delete_layers(f1, "GeomPoint")
f2 <- delete_layers(f2, "GeomPoint")
f3 <- delete_layers(f3, "GeomPoint")

the_plot <- grid_arrange_shared_legend(tpr1, tpr2, tpr3, tnr1, tnr2, tnr3, f1, f2, f3, ncol=3, nrow=3, position="bottom")

ggsave(file="<DIR>/complete.svg", plot=the_plot, width=16, height=9)
