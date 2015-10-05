################################################################################
#                                                                              #
# Functions for calculating metrics of diversity, evenness, rarity, etc.       #
# These are not included in other diversity packages, e.g., Vegan              #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Ken Locey                                                        #
# Obtained from:                                                               #
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/      #
#                                                                              #
################################################################################
#                                                                              #
# Recent Changes:                                                              #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1.                                                                   #
#         2. Add warnings                                                      #
#                                                                              #
################################################################################

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


theme_black <- function(base_size = 12, base_family = "Helvetica") {
  theme(
    line =               element_line(colour = "black", size = 1.5, linetype = 1,
                                      lineend = "butt"),
    rect =               element_rect(fill = "white", colour = "black", size = 1.5, linetype = 1),
    text =               element_text(family = base_family, face = "plain",
                                      colour = "black", size = base_size,
                                      hjust = 0.5, vjust = 1.2, angle = 0, lineheight = 0.9),
    axis.text =          element_text(size = rel(0.8), colour = "white"),
    strip.text =         element_text(size = rel(0.8), colour = "white"),
    
    axis.line =          element_blank(),
    axis.text.x =        element_text(vjust = 1),
    axis.text.y =        element_text(hjust = 1),
    axis.ticks =         element_line(colour = "white", size = 0.8),
    axis.title =         element_text(colour = "white"),
    axis.title.x =       element_text(vjust = -0.5, size = rel(1.5)),
    axis.title.y =       element_text(angle = 90, size = rel(1.5)),
    #axis.ticks.length =  unit(0.3, "lines"),
    #axis.ticks.margin =  unit(0.5, "lines"),
    
    legend.background =  element_rect(colour = "black"),
    #legend.margin =      unit(0.2, "cm"),
    legend.key =         element_rect(fill = "black", colour = "white"),
    #legend.key.size =    unit(1.2, "lines"),
    legend.key.height =  NULL,
    legend.key.width =   NULL,
    legend.text =        element_text(size = rel(0.8), colour = "white"),
    legend.text.align =  NULL,
    legend.title =       element_text(size = rel(0.8), face = "bold", hjust = 0, colour = "white"),
    legend.title.align = NULL,
    legend.position =    "none",
    legend.direction =   "vertical",
    legend.justification = "center",
    legend.box =         NULL,
    #opts(legend.position = "none"),
    
    panel.background =   element_rect(fill = "black", colour = NA),
    panel.border =       element_rect(fill = NA, colour = "white"),
    panel.grid.major =   element_line(colour = "grey20", size = 0.2),
    panel.grid.minor =   element_line(colour = "grey5", size = 0.5),
    #panel.margin =       unit(0.25, "lines"),
    
    strip.background =   element_rect(fill = "grey30", colour = "grey10"),
    strip.text.x =       element_text(),
    strip.text.y =       element_text(angle = -90),
    
    plot.background =    element_rect(colour = "black", fill = "black"),
    plot.title =         element_text(colour = "white", size = rel(1.8)),
    #plot.margin =        unit(c(1, 1, 0.5, 0.5), "lines"),
    
    complete = TRUE
  )
}
# Check that it is a complete theme
attr(theme_black(), "complete")



addlinetoplot <- function(varx, vary) { 
  list(
    geom_line(data=NULL, aes(x=varx, y=vary)), 
    geom_point(data=NULL, aes(x=varx, y=vary))
  )
}