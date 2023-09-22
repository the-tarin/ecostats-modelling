set.seed(123)
library(lubridate)
library(circular)

# x and y in square survey site with dimensions 10km x 10km with markers every 2000m
# mics
x_mic = seq(from = 2000, to = 8000, by = 2000)
y_mic = seq(from = 2000, to = 8000, by = 2000)

mic_coords = expand.grid(x = x_mic, y = y_mic)
colnames(mic_coords) = c('x_mic', 'y_mic')
mic_coords = as.matrix(mic_coords)

# expectation of gibbon groups in survey area
intensity = 0.000001 # modelling parameter
x_perimeter = 10000
y_perimeter = 10000
survey_area = x_perimeter * y_perimeter
expected_gibbon_groups = intensity * survey_area

# homogeneous Poisson point process for getting gibbon coordinates
x_gibbon_group <- runif(expected_gibbon_groups, 0, x_perimeter)
y_gibbon_group <- runif(expected_gibbon_groups, 0, y_perimeter)

gibbon_group_coords = cbind(x_gibbon_group, y_gibbon_group)
gibbon_group_coords = as.matrix(gibbon_group_coords)

# plot of mic grid and Poisson point process distributed gibbon groups
{
plot(x_gibbon_group, y_gibbon_group, type = "p", col = "blue", pch = 1, xlab = "X-coordinate (m)", ylab = "Y-coordinate (m)")
points(mic_coords[,1], mic_coords[,2], type = "p", col = "black", pch = 15)
}

# Calculate the bearing from mic to gibbon group
calc_bearing = function(x_point_1, y_point_1, x_point_2, y_point_2) {
  bearing = atan2(
    y_point_2 - y_point_1,
    x_point_2 - x_point_1
  )
  
  # convert to range of radians [0,2pi] where bearing is taken clockwise from north (0 radians)
  if (bearing <= 0) {
    bearing = -1*bearing + 0.5*pi 
  } else if ((bearing > 0) && (bearing <= 0.5*pi)) {
    bearing = (0.5*pi - bearing)
  } else {
    bearing = (pi - bearing) + 1.5*pi
  }
  return(bearing)
}


# euclidean distance for detector to gibbon observation
dist_mic_gibbon <- matrix(NA, nrow = nrow(mic_coords), ncol = nrow(gibbon_group_coords))
bearing_mic_gibbon <- matrix(NA, nrow = nrow(mic_coords), ncol = nrow(gibbon_group_coords))

for (i in 1:nrow(dist_mic_gibbon)) {
  for (j in 1:ncol(dist_mic_gibbon)) {
    dist_mic_gibbon[i, j] <- sqrt(sum((mic_coords[i, ] - gibbon_group_coords[j, ])^2))
    bearing_mic_gibbon[i, j] <- calc_bearing(mic_coords[i, 1], mic_coords[i, 2], gibbon_group_coords[j, 1], gibbon_group_coords[j, 2])
  }
}

# half normal hazard distribution for probability matrix
# assuming each mic has same detection function
lambda = 50 # modelling parameter
sigma = 500 # modelling parameter

hazard_half_normal = function(lambda, sigma, dist_matrix) {
  prob_matrix = 1 - exp(-lambda * exp((-1 * dist_matrix^2) / (2 * sigma^2)))
  return (prob_matrix)
}

detection_prob_matrix = hazard_half_normal(lambda = lambda, sigma = sigma, dist_matrix = dist_mic_gibbon)

# weighted Bernoulli trails for each gibbon group to mic for detected/undetected
detection_matrix <- matrix(NA, nrow = nrow(mic_coords), ncol = nrow(gibbon_group_coords))

for (i in 1:nrow(dist_mic_gibbon)) {
  for (j in 1:ncol(dist_mic_gibbon)) {
    detection_matrix[i, j] <- rbinom(size = 1, n = 1, prob = detection_prob_matrix[i,j])
  }
}

### plotting results

# plot detected gibbon groups for chosen mic
plot_select = 12

for (i in plot_select) {
  detection_gibbon_group_idx = which(detection_matrix[i,] == 1)
  plot(x_gibbon_group, y_gibbon_group, type = "p", col = "blue", pch = 1, xlab = "X-coordinate (m)", ylab = "Y-coordinate (m)")
  points(mic_coords[,1], mic_coords[,2], type = "p", col = "black", pch = 15)
  points(mic_coords[i,1], mic_coords[i,2], type = "p", col = "red", pch = 15)
  points(x_gibbon_group[detection_gibbon_group_idx], y_gibbon_group[detection_gibbon_group_idx], type = "p", col = "green", pch = 15)
}

# plot mics which detected a certain gibbon group
plot_select = 36

for (i in plot_select) {
  detection_mic_idx = which(detection_matrix[,i] == 1)
  plot(x_gibbon_group, y_gibbon_group, type = "p", col = "blue", pch = 1, xlab = "X-coordinate (m)", ylab = "Y-coordinate (m)")
  points(mic_coords[,1], mic_coords[,2], type = "p", col = "black", pch = 15)
  points(x_gibbon_group[i], y_gibbon_group[i], type = "p", col = "red", pch = 15)
  points(mic_coords[detection_mic_idx,1], mic_coords[detection_mic_idx,2], type = "p", col = "green", pch = 15)
}

### create dataframes for analysis

# convert meter coordinates to lat, lon for RShiny App
convert_coords = function(x_coord, y_coord) {
  ref_lat = 14.238574
  ref_lon = 70
  earth_radius = 6371000
  
  # calculate lat, lon
  # issue with this calculation
  lat = ref_lat + (y_coord / earth_radius) * (180 / pi)
  lon = ref_lon + (x_coord / (earth_radius * cos(lat * pi / 180))) * (180 / pi)
  
  lat_lon = cbind(lat, lon)
  
  return(lat_lon)
}

mic_lat_lon = data.frame()
mic_lat_lon = convert_coords(mic_coords[,1], mic_coords[,2])
mic_df = cbind(1:nrow(mic_lat_lon), mic_lat_lon)
colnames(mic_df)[1] = "mic_id"

# these are ground truth locations of each gibbon group and exact call times
gibbon_group_df = cbind(1:nrow(gibbon_group_coords), gibbon_group_coords)
colnames(gibbon_group_df)[1] = "gibbon_group_id"

# generate random ground truth call times
dawn_chorus_start_time = as.POSIXct("2023-01-01 04:00:00", tz = "UTC")
dawn_chorus_end_time = as.POSIXct("2023-01-01 05:00:00", tz = "UTC")

call_datetimes = as.POSIXct(runif(nrow(gibbon_group_df), dawn_chorus_start_time, dawn_chorus_end_time), origin = "1970-01-01", tz = "UTC")
call_timestamps = as.integer(call_datetimes)

# gibbon_group_df = cbind(gibbon_group_df, call_datetimes)
# colnames(gibbon_group_df)[4] = "call_datetime"

### generate recording dataframe

# speed of sound (m/s)
speed_of_sound = 331
von_mises_kappa = 50 # modelling parameter

recording_df = data.frame()

for (i in 1:nrow(mic_df)) {
  ground_truth_animal_ID = which(detection_matrix[i,] == 1)
  mic_ID = rep(i, times = length(ground_truth_animal_ID))
  
  # compute measured call time which considers speed of sound
  ground_truth_call_timestamp = call_timestamps[as.integer(gibbon_group_df[ground_truth_animal_ID])]
  ground_truth_call_datetime = as.POSIXct(ground_truth_call_timestamp, origin = "1970-01-01", tz = "UTC")
  measured_call_timestamp = ground_truth_call_timestamp + dist_mic_gibbon[i, ground_truth_animal_ID] / speed_of_sound
  measured_call_datetime = as.POSIXct(measured_call_timestamp, origin = "1970-01-01", tz = "UTC")
  
  # compute measured bearing
  ground_truth_bearing = bearing_mic_gibbon[i, ground_truth_animal_ID]
  # von mises distribution for bearing measurement
  measured_bearing = data.frame()
  for (j in 1:length(ground_truth_animal_ID)) {
    bearing = rvonmises(n = 1, mu = ground_truth_bearing[j], kappa = von_mises_kappa)
    measured_bearing = rbind(measured_bearing, bearing)
  }
  colnames(measured_bearing) = "measured_bearing"
  
  # binding dataframes for each mic
  recording_temp = cbind(mic_ID, ground_truth_animal_ID, ground_truth_call_datetime, measured_call_datetime, ground_truth_bearing, measured_bearing)
  recording_df = rbind(recording_df, recording_temp)
}

write.csv(mic_df, "../output/mic.csv", row.names=FALSE)
write.csv(gibbon_group_df, "../output/gibbon_group.csv", row.names=FALSE)
write.csv(recording_df, "../output/recording.csv", row.names=FALSE)

### map plots for testing lat/lon
plot(mic_df[,2], mic_df[,3], type = "p", pch = 19, col = "blue", xlab = "Longitude", ylab = "Latitude", main = "Scatterplot of Coordinates")
