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
hotspot_map <- function(in_df, fn, color_samples = 100, p = 0.1, dim = in_dim, bounds = NULL, pbar = TRUE) {

	data_df <- fn(in_df)

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
		random_df <- fn(random_df_in)

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

	current_index <- 1
	if (sum(data_df$score < sq[1]) > 0) {
		color_df$score[data_df$score < sq[1]] <- current_index
		included_colors <- append(included_colors, colors[1])
		current_index <- current_index + 1
	}
	
	
	for (i in 1:(length(sq)-1)) {
		score_indices <- (data_df$score >= sq[i]) & (data_df$score < sq[i+1])
		if (sum(score_indices) > 0) {
			color_df$score[score_indices] <- current_index
			included_colors <- append(included_colors, colors[i+1])
			current_index <- current_index + 1
		}
	}

	if (sum(data_df$score > sq[length(sq)]) > 0) {
		color_df$score[data_df$score > sq[length(sq)]] <- current_index
		included_colors <- append(included_colors, colors[length(colors)])
	}	


	
	color_df$score <- as.factor(color_df$score)
	return(list("data" = color_df, "colors" = included_colors, "rawscore" = data_df))
}	

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
			print(datad$x >= x_vals[j])
			case_indices <- append(case_indices, all_indices[datad$x >= x_vals[j] & datad$x <= x_vals[j+1] & datad$y >= y_vals[i] & datad$y <= y_vals[i+1]])
		}
	}
	rdf <- datad
	rdf$in_cluster <- rep(0, nrow(datad))
	rdf$in_cluster[case_indices] <- 1
	return(rdf)
}

