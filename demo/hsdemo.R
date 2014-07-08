require(hotspotr)
require(ggplot2)

sensplot <- 20
resisplot <- 3
pointsize <- 3

x <- runif(1000)
y <- runif(1000)

hs <- random_hotspot(x,y, 0.30, 0.8, 0.2)

cr <- hs[["coord"]]
hs <- data.frame(x = x, y = y, z = as.factor(hs[["z"]]))


x <- hotspot_map(hs, dbm_score_rr, p = 0.03, color_samples = 100, pbar = TRUE)
p <- plot(x)
ggsave("test.png")