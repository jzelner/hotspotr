#' hsmap
#'
#' Create a hotspot map object
#' @export
hsmap <- function(levels, colors, points) {
	x <- structure(list(levels = levels, colors = colors, points = points), class = "hsmap")
}

#' plot.hsmap
#'
#' @param x object of type hsmap
#' @export
plot.hsmap <- function(x, point = TRUE, pointtype = c(20, 3), xlab = "Longitude", ylab = "Latitude", map = NULL) {

	if (!is.null(map)) {
		g <- map 
	} else {
		g <- ggplot() 
	}

	if (point == TRUE) {
		g <- g + geom_point(aes(x = x$points$x, y = x$points$y, alpha = 0.01, shape = factor(x$points$z), size = 1)) + scale_shape_manual(values = c(pointtype[1],pointtype[2])) + scale_size_identity()
	}


	g <- g + geom_tile(aes(x =x, y = y, fill = score, alpha = 0.01), data = x$levels) + scale_fill_manual(values = x$colors) + coord_cartesian(xlim = c(min(x$levels$x), max(x$levels$x)), ylim = c(min(x$levels$y), max(x$levels$y))) + xlab(xlab) + ylab(xlab)  


	g <- g + theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank()) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + guides(fill = FALSE, alpha = FALSE) 

	print(g)
}

#' Generate a hotspot map
#'
#' This function takes a function for computing some kind of spatial risk,
#' calculates spatial densities and generates a color scale indicating which areas
#' have a statistically significant concentration of cases.
#' @param in_df data frame consisting of x and y coordinates of input data, and any 
#' other information used by the mapping function.
#' @param fn mapping function that takes a data frame as input and returns a data frame #' with a grid of points of size dim x dim.  
#' @param bounds gives the coordinates of a bounding box to restrict the map to (otherwise uses the extreme values of points in in_df) 
#' @param pbar show a progress bar while calculating the color scale. 
#' @export
hotspot_map <- function(in_df, fn, color_samples = 100, p = 0.005, dim = in_dim, bounds = NULL, pbar = TRUE) {

	data_df <- fn(in_df, p = p)

	#Now repeat with the same dataset with case and control statuses randomized
	random_df_in <- in_df
	max_scores <- c()
	min_scores <- c()
	all_scores <- c()

	message(sprintf("Sampling color scale %d times",c(color_samples)))
	if (pbar == TRUE) {
		pb <- txtProgressBar()
	}
	for (i in 1:color_samples) {
		#Randomly permute case/control labels
		random_df_in$z <- sample(random_df_in$z)
		random_df <- fn(random_df_in, p = p)

		min_scores <- append(min_scores, min(random_df$score))
		max_scores <- append(max_scores, max(random_df$score))
		if (pbar == TRUE) {
			setTxtProgressBar(pb, i/color_samples)
		}
	}
	if (pbar == TRUE) {
		close(pb)
	}

	all_scores <- c(min_scores, max_scores)
	# print(all_scores)
	color_df <- data_df
	sq <- quantile(all_scores, probs = seq(0.0,1.0,0.025))

	included_colors <- c()
	cr <- colorRampPalette(c("turquoise", "orange"))

	colors <- c("blue4", "blue3")
	colors <- append(colors, cr(length(sq)-3)) 
	colors <- append(colors, c("red3", "red4"))

	pvals <- append(-1,append(seq(0.0,1.0,0.025),1.25))
	color_df$pval <- rep("",nrow(color_df))

	current_index <- 1
	if (sum(data_df$score < sq[1]) > 0) {
		color_df$score[data_df$score < sq[1]] <- current_index
		color_df$pval[data_df$score < sq[1]] <- pvals[1]

		included_colors <- append(included_colors, colors[1])
		current_index <- current_index + 1
	}


	for (i in 1:(length(sq)-1)) {
		score_indices <- (data_df$score >= sq[i]) & (data_df$score < sq[i+1])
		if (sum(score_indices) > 0) {
			color_df$score[score_indices] <- current_index
			color_df$pval[score_indices] <- pvals[i+1]
			included_colors <- append(included_colors, colors[i+1])
			current_index <- current_index + 1
		}
	}

	if (sum(data_df$score > sq[length(sq)]) > 0) {
		color_df$score[data_df$score > sq[length(sq)]] <- current_index
		color_df$pval[data_df$score > sq[length(sq)]] <- pvals[length(pvals)]

		included_colors <- append(included_colors, colors[length(colors)])
	}	



	color_df$pval <- as.numeric(color_df$pval)
	color_df$score <- as.factor(color_df$score)
	return(hsmap(color_df, included_colors, in_df))
}	
#' Generate a random hotspot.
#'
#' Given a set of points, generate a square hotspot with width h in the middle. 
#' Useful for testing whether a map based on a given set of points will work. 
#' Specify the within-spot risk (in_p) and outside-spot risk (out_p), so that the 
#' relative risk (RR) is equal to in_p / out_p.

#' @param x vector of x coordinates for input points
#' @param y vector of y coordinates for input points
#' @param h2 width of hotspot window
#' @param in_p risk of being in case group inside hotspot
#' @param out_p risk of being in case group outside hotspot
#' @export
random_hotspot <- function(x,y, h2, in_p, out_p) {
	
	h <- h2/2
	mid_x <- (max(x) + min(x))/2
	mid_y <- (max(y) + min(y))/2

	hp_x <- c(mid_x-h, mid_x+h)
	hp_y <- c(mid_y-h, mid_y+h)

	in_spot <- (x > hp_x[1] & x < hp_x[2]) & (y > hp_y[1] & y < hp_y[2])

	z <- rbinom(length(x), 1, out_p)
	z[in_spot] <- rbinom(sum(in_spot), 1, in_p)

	vals <- list(z = z, coord = c(hp_x, hp_y) )
	return(vals)

}

#' Given a map with a single hotspot, extract the squares with an approximate p-value
#' greater than or equal to a given threshhold
#'
#' @param df color dataframe output from make_hotspot
#' @param datad data frame containing data points
#' @param pv threshhold p-value
#' @export
hotspot_area <- function(df, datad, pv) {
	case_indices <- c()
	all_indices <- 1:nrow(datad)
	hs_rows <- subset(df, pval >= pv)
	#Scan across unique y values
	y_vals <- sort(unique(hs_rows$y))
	for (i in 1:(length(y_vals)-1)) {
		x_df <- subset(hs_rows, y >= y_vals[i] & y <= y_vals[i+1])
		x_vals <- sort(unique(x_df$x))
		for (j in 1:(length(x_vals)-1)) {
			case_indices <- append(case_indices, all_indices[datad$x >= x_vals[j] & datad$x <= x_vals[j+1] & datad$y >= y_vals[i] & datad$y <= y_vals[i+1]])
		}
	}
	rdf <- datad
	rdf$in_cluster <- rep(0, nrow(datad))
	rdf$in_cluster[case_indices] <- 1
	return(rdf)
}

