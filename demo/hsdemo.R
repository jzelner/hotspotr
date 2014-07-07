require(hotspotr)
require(ggplot2)
require(pracma)

sensplot <- 20
resisplot <- 3
pointsize <- 3

x <- runif(1000)
y <- runif(1000)

hs <- random_hotspot(x,y, 0.30, 0.8, 0.2)

cr <- hs[["coord"]]
hs <- data.frame(x = x, y = y, z = as.factor(hs[["z"]]))


x <- hotspot_map(hs, dbm_score_rr, p = 0.03, color_samples = 100, pbar = FALSE)

plot(x)

# retreat_color_df <- retreat_color_d[["data"]]
# retreat_colors <- retreat_color_d[["colors"]]

# retreat_point <- retreat_point + geom_tile(aes(x =x, y = y, fill = score, alpha = 0.01), data = retreat_color_df) + scale_fill_manual(values = retreat_colors) + coord_cartesian(xlim = c(min(retreat_color_df$x), max(retreat_color_df$x)), ylim = c(min(retreat_color_df$y), max(retreat_color_df$y))) + xlab("Longitude") + ylab("Latitude")  
# print(retreat_point)
