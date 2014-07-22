#' hsmap
#'
#' Create a hotspot map object
#' @export
hsmap <- function(levels, colors, points) {
	x <- structure(list(levels = levels, colors = colors, points = points), class = "hsmap", fname = NULL)
}

#' plot.hsmap
#'
#' @param x object of type hsmap
#' @export
plot.hsmap <- function(x, point = TRUE, caseonly = FALSE, pointtype = c(20, 3), xlab = "Longitude", ylab = "Latitude", map = NULL, fname = NULL) {

	if (!is.null(map)) {
		g <- map
	} else {
		.e <- environment()
		g <- ggplot(environment = .e) + coord_cartesian(xlim = c(min(x$levels$x), max(x$levels$x)), ylim = c(min(x$levels$y), max(x$levels$y))) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
	}

	if (point == TRUE) {

		point_df <- data.frame(x = x$points$x, y = x$points$y, z = as.factor(x$points$z), size = 1)

		if (caseonly == TRUE) {
			point_df <- subset(point_df, z == 1)
		}
		g <- g + geom_point(aes(x = x, y = y, shape = z, size = 1), data = point_df) + scale_shape_manual(values = c(pointtype[1],pointtype[2])) + scale_size_identity()
	}


	g <- g + geom_tile(aes(x =x, y = y, fill = score, alpha = 0.01), data = x$levels) + scale_fill_manual(values = x$colors) + xlab(xlab) + ylab(xlab)


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
      plot.background=element_blank()) + guides(fill = FALSE, alpha = FALSE)

	if (!is.null(fname)) {
		ggsave(fname)
	} else {
		return(g)
	}
}

get_colorscale <- function(low_colors = c("blue4", "blue3"), high_colors = c("red3", "red4"), inner_range = c("turquoise", "orange")) {
	cr <- colorRampPalette(inner_range)
	inner_colors <- cr(90)
	all_colors <- append(low_colors[1], append(rep(low_colors[2],4), inner_colors))

	all_colors <- append(append(all_colors, rep(high_colors[1],4)), high_colors[2])

	return(all_colors)
}

get_score_percentile <- function(scoredist, scores) {
	d <- ecdf(scoredist)
	p <- round(d(scores),2)
	p[p == 0] <- 0.01
	return(p)
}

make_hsmap <- function(data_df, scores, in_df) {

		cs <- get_colorscale()

		pval <- get_score_percentile(scores, data_df$score)

		data_df$pval <- pval
		data_df$score <- as.factor(pval)

		included_colors <- cs[sort(unique(pval*100))]

		return(hsmap(data_df, included_colors, in_df))
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

	data_df <- fn(in_df, p = p, bounds = bounds)

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

	scores <- c(min_scores, max_scores)
	hm <- make_hsmap(data_df, scores, in_df)

	return(hm)
}

#' Make a hotspot map using julia
#' @export
hotspot_map_julia <- function(in_df, r = 1.0, p = 0.005, dim = 100) {
	require("rjson")

	max_x <- max(in_df$x)
	min_x <- min(in_df$x)
	max_y <- max(in_df$y)
	min_y <- min(in_df$y)

	out_d <- list()
	out_d[["x"]] = (in_df$x - min(in_df$x)) / (max(in_df$x) - min(in_df$x))
	out_d[["y"]] = (in_df$y - min(in_df$y)) / (max(in_df$y) - min(in_df$y))
	out_d[["z"]] = as.numeric(as.character(in_df$z))
	out_d[["r"]] = r
	out_d[["p"]] = p
	out_d[["dim"]] = dim

	write(toJSON(out_d), "test.json")
	
	cmd <- "julia -e \'using Hotspot; z = Hotspot.hsmap_from_json(\"test.json\"); f = open(\"out.json\",\"w\"); write(f, z); close(f);\'"

	system(cmd)

	map_d <- fromJSON(file="out.json")

	data_df <- data.frame(x = ((max_x-min_x)*map_d[["x"]])+min_x, y = ((max_y-min_y)*map_d[["y"]])+min_y, score = map_d[["score"]])

	scores <- c(map_d[["min"]], map_d[["max"]])

	hm <- make_hsmap(data_df, scores, in_df)

	return(hm)
}

#' Make multiple hotspot maps using julia
#' @export
batch_hotspot_map_julia <- function(in_maps, in_df, r = 1.0, p = 0.005, dim = 100) {
	require("rjson")
	max_x_vals <- c()
	min_x_vals <- c()
	max_y_vals <- c()
	min_y_vals <- c()

	processed_maps <- list()
	for (i in 1:length(in_maps)) {
		in_df <- in_maps[[i]]
		max_x <- max(in_df$x)
		max_x_vals <- append(max_x_vals, max_x)
		
		min_x <- min(in_df$x)
		min_x_vals <- append(min_x_vals, min_x)

		max_y <- max(in_df$y)
		max_y_vals <- append(max_y_vals, max_y)

		min_y <- min(in_df$y)
		min_y_vals <- append(min_y_vals, min_y)

		out_d <- list()
		out_d[["x"]] = (in_df$x - min(in_df$x)) / (max(in_df$x) - min(in_df$x))
		out_d[["y"]] = (in_df$y - min(in_df$y)) / (max(in_df$y) - min(in_df$y))
		out_d[["z"]] = as.numeric(as.character(in_df$z))
		out_d[["r"]] = r
		out_d[["p"]] = p
		out_d[["dim"]] = dim

		processed_maps[[i]] <- out_d
	}

	cmd <- "julia -e \'using Hotspot; z = Hotspot.hsmap_from_json(\"test.json\"); f = open(\"out.json\",\"w\"); write(f, z); close(f);\'"

	write(toJSON(processed_maps), "test.json")
	system(cmd)

	map_list <- fromJSON(file="out.json")

	return_maps <- list()
	for (i in 1:length(map_list)) {
		map_d <- map_list[[i]]

		data_df <- data.frame(x = ((max_x-min_x)*map_d[["x"]])+min_x, y = ((max_y-min_y)*map_d[["y"]])+min_y, score = map_d[["score"]])

		scores <- c(map_d[["min"]], map_d[["max"]])

		return_maps[[i]] <- make_hsmap(data_df, scores, in_maps[[i]])
	}

	return(return_maps)
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
