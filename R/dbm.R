require(pracma)

#' Get distance of every point in data from fixed point (cx, cy)
#' @export
all_dist <- function(x, y, cx,cy) {
	return(sqrt((x-cx)**2 + (y-cy)**2))
}

#Create a data frame with (x,y) coordinates of N uniformly spaced cells
risk_grid <- function(minX, maxX, minY, maxY, dim) {
	xvals <- seq(from = minX, to = maxX, by = ((maxX-minX)/(dim-1)))
	yvals <- seq(from = minY, to = maxY, by = ((maxY-minY)/(dim-1)))
	z <- meshgrid(xvals, yvals)
	df <- data.frame(x = as.vector(z$X), y = as.vector(z$Y))
	return(df)
}

#' Create a square grid of points of the pre-specified dimension
#'
#' @param x vector of x coordinates of data points
#' @param y vector of y coordinates of data points
#' @param dim dimension of resulting square grid
#' @export
grid_from_points <- function(x,y, dim) {
	rg <- risk_grid(min(x), max(x), min(y), max(y), dim)
	return(rg)
}

#' Calculate an empirical cdf of point distances relative to a fixed point.
#'
#' @param x vector of x coordinates of data points
#' @param y vector of y coordinates of data points
#' @param cx fixed point x coordinate
#' @param cy fixed point y coordinate
#' @export
fixed_point_cdf <- function(x, y, cx, cy) {
	d <- all_dist(x,y,cx,cy)
	dcdf <- ecdf(d)
	return(dcdf)
}

corner_points <- function(x, y, h = 0.1) {
	x <- na.omit(x)
	y <- na.omit(y)

	xout <- c(min(x)-h, min(x)-h, max(x)+h, max(x)+h)
	yout <- c(min(y)-h, max(y)+h, min(y)-h, max(y)+h)
	return(data.frame(x = xout, y = yout))
}

#' Surround a set of points with a circle of points at a pre-specified
#' radius.
#'
#' @param x vector of x coordinates for points to surround
#' @param y vector of y coordinates for points to surround
#' @param np number of external points
#' @param r radius of circle, measured in distance from center of the 
#' set of points.
#' @export
circle_points <- function(x, y, np = 6, r = 0.5) {
	x <- na.omit(x)
	y <- na.omit(y)

	mid_x <- (max(x) + min(x))/2
	mid_y <- (max(y) + min(y))/2

	offsets <- seq(from = 0, to = (2*pi)-(2*pi/20),  by = (2*pi/20))
	xout <- c()
	yout <- c()
	for (i in offsets) {
		xout <- append(xout, mid_x + r*sin(i))
		yout <- append(yout, mid_y + (r*cos(i)))	
	}
	
	return(data.frame(x = xout, y = yout))
}

adaptive_dbm_score <- function(x, y, dim = 100, p0 = 0.1) {

	#Make the grid of points we're going to use
	map_points <- grid_from_points(x,y,dim)

	#Get a set of points outside the region
	# ref_points <- corner_points(x,y,h)
	ref_points <- circle_points(x,y,1.0)

	for (i in 1:nrow(ref_points)) {

		fpx <- ref_points$x[i]
		fpy <- ref_points$y[i]

		#Now get the CDFs corresponding to these points
		scorefunc <- fixed_point_cdf(x, y, fpx, fpy)

		#Now translate the grid coordinates into distances from the fixed 
		#point
		distances <- all_dist(map_points$x,map_points$y,fpx,fpy)

		dist_quantile <- scorefunc(distances) 

		#Get low and high quantiles, assuming that d-(p0/2) > 0 and d+(p0/2) < 1
		low_q <- dist_quantile - (p0/2)
		high_q <- dist_quantile + (p0/2)

		#Now correct quantiles
		high_q[low_q < 0] <- high_q[low_q < 0] - low_q[low_q < 0]
		low_q[low_q < 0] <- 0

		low_q[high_q > 1] <- low_q[high_q > 1] - (high_q[high_q > 1]-1)
		high_q[high_q > 1] <- 1

		s0_diff <- as.vector(quantile(scorefunc, high_q) -quantile(scorefunc, low_q))
		

		if (i == 1) {
			rd <- data.frame(x = map_points$x, y = map_points$y, score = s0_diff)
		} else {
			rd$score <- rd$score + s0_diff
		}
	}

	return(rd)

}

#' Generate a scored map of the local proportion of cases.
#'
#' Returns a data frame with x and y coordinates of grid squares of 
#' pre-specified dimension and the dbm score for each square.
#'
#' @param df data frame containing (x,y) coordinates of cases and case control status (z: 0,1)
#' @export
dbm_score_rr <- function(df, dim = 100, h = 0.005, p = 0.005, bounds = NULL) {

	x <- df$x
	y <- df$y

	#Subset out the cases
	c_df <- subset(df, z == 1)
	cx <- c_df$x
	cy <- c_df$y

	#Get the set of coordinates for the grid for the map
	if (is.null(bounds)) {
		map_points <- grid_from_points(x,y,dim)
	} else {
		map_points <- grid_from_points(bounds[1:2],bounds[3:4],dim)		
	}

	#Get a set of points outside the region
	ref_points <- circle_points(x,y,h)

	#For each outside reference point, evaluate the score function and store
	for (i in 1:nrow(ref_points)) {

		fpx <- ref_points$x[i]
		fpy <- ref_points$y[i]

		#Now get the CDFs corresponding to these points
		scorefunc_0 <- fixed_point_cdf(x, y, fpx, fpy)
		scorefunc_1 <- fixed_point_cdf(cx, cy, fpx, fpy)

		#Now translate the grid coordinates into distances from the fixed 
		#point
		distances <- all_dist(map_points$x,map_points$y,fpx,fpy)

		s0_diff <- scorefunc_0(distances+p)-scorefunc_0(distances-p)
		
		s1_diff <- scorefunc_1(distances+p)-scorefunc_1(distances-p)


		if (i == 1) {
			rd <- data.frame(x = map_points$x, y = map_points$y, score = s1_diff - s0_diff)
		} else {
			rd$score <- rd$score + (s1_diff - s0_diff)
		}
	}
	
	return(rd)

}

adaptive_dbm_score_rr <- function(df, dim = 100, h = 2.0, p0 = 0.1) {

	x <- df$x
	y <- df$y

	#Subset out the cases
	c_df <- subset(df, z == 1)
	cx <- c_df$x
	cy <- c_df$y

	#Get the set of coordinates for the grid for the map
	map_points <- grid_from_points(x,y,dim)

	#Get a set of points outside the region
	ref_points <- circle_points(x,y,h)

	#For each outside reference point, evaluate the score function and store
	for (i in 1:nrow(ref_points)) {

		fpx <- ref_points$x[i]
		fpy <- ref_points$y[i]

		#Now get the CDFs corresponding to these points
		scorefunc_0 <- fixed_point_cdf(x, y, fpx, fpy)
		scorefunc_1 <- fixed_point_cdf(cx, cy, fpx, fpy)

		#Translate the grid coordinates into distances from the fixed 
		#point
		distances <- all_dist(map_points$x,map_points$y,fpx,fpy)

		dist_quantile <- scorefunc_0(distances) 

		#Get low and high quantiles, assuming that d-(p0/2) > 0 and d+(p0/2) < 1
		low_q <- dist_quantile - (p0/2)
		high_q <- dist_quantile + (p0/2)

		#Now correct quantiles
		high_q[low_q < 0] <- high_q[low_q < 0] - low_q[low_q < 0]
		low_q[low_q < 0] <- 0

		low_q[high_q > 1] <- low_q[high_q > 1] - (high_q[high_q > 1]-1)
		high_q[high_q > 1] <- 1


		# s0_diff <- as.vector(quantile(scorefunc_0, high_q) -quantile(scorefunc_0, low_q))

		s0_diff <- p0

		s1_diff <- as.vector(scorefunc_1(quantile(scorefunc_0, high_q)) -scorefunc_1(quantile(scorefunc_0, low_q)))


		if (i == 1) {
			rd <- data.frame(x = map_points$x, y = map_points$y, score = s1_diff - s0_diff)
		} else {
			rd$score <- rd$score + (s1_diff - s0_diff)
		}
	}

	rd$score <- rd$score / nrow(ref_points)
	
	return(rd)

}