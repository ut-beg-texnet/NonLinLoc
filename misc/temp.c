	/* incorporate station weights */
	
	if (iSetStationDistanceWeights) {

		for (nrow = 0; nrow < num_arrivals; nrow++) {
			for (ncol = 0; ncol < num_arrivals; ncol++) {
				 wt_matrix[nrow][ncol] *= arrival[nrow].station_weight * arrival[ncol].station_weight;
				 edt_matrix[nrow][ncol] *= arrival[nrow].station_weight * arrival[ncol].station_weight;
			}
		}
	
		if (message_flag >= 5)
			DisplayDMatrix("Weight (after Station Distance Weighting)", wt_matrix, num_arrivals, num_arrivals);
				
	
	}
	

