# cite: http://stat.ethz.ch/R-manual/R-patched/library/stats/html/reshape.html
# cite: http://www.cs.grinnell.edu/~rebelsky/Courses/MAT115/2008S/R/stacked-bar-graphs.html
# cite: http://www.harding.edu/fmccown/r/#barcharts

library(RColorBrewer)
sequential <- brewer.pal(6, "BuGn")

Loblolly[1:10,]

wide <- reshape(Loblolly,
                v.names = "height",
                timevar = "age",
                idvar = "Seed",
                direction = "wide")

wide[1:5,]

wide$h0.3 <- wide$height.3
wide$h3.5 <- wide$height.5 - wide$height.3
wide$h5.10 <- wide$height.10 - wide$height.5
wide$h10.15 <- wide$height.15 - wide$height.10
wide$h15.20 <- wide$height.20 - wide$height.15
wide$h20.25 <- wide$height.25 - wide$height.20

barplot(t(wide[,8:13]),
        names.arg = wide$Seed, # x-axis labels
        cex.names = 0.7, # makes x-axis labels small enough to show all
        col = sequential, # colors
        xlab = "Seed Source",
        ylab = "Height, Feet",
        xlim = c(0,20), # these two lines allow space for the legend
        width = 1) # these two lines allow space for the legend
legend("bottomright",
       legend = c("20-25", "15-20", "10-15", "5-10", "3-5", "0-3"), #in order from top to bottom
       fill = sequential[6:1], # 6:1 reorders so legend order matches graph
       title = "Years")
