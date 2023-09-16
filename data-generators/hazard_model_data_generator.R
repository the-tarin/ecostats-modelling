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
x_gibbon_group <- round(runif(expected_gibbon_groups, 0, x_perimeter))
y_gibbon_group <- round(runif(expected_gibbon_groups, 0, y_perimeter))

gibbon_group_coords = cbind(x_gibbon_group, y_gibbon_group)
gibbon_group_coords = as.matrix(gibbon_group_coords)

# plot of mic grid and Poisson point process distributed gibbon groups
{
plot(x_gibbon_group, y_gibbon_group, type = "p", col = "blue", pch = 1, xlab = "X-coordinate (m)", ylab = "Y-coordinate (m)")
points(mic_coords[,1], mic_coords[,2], type = "p", col = "black", pch = 15)
}

# euclidean distance for detector to gibbon observation
dist_mic_gibbon <- matrix(NA, nrow = nrow(mic_coords), ncol = nrow(gibbon_group_coords))
bearing_mic_gibbon <- matrix(NA, nrow = nrow(mic_coords), ncol = nrow(gibbon_group_coords))

for (i in 1:nrow(dist_mic_gibbon)) {
  for (j in 1:ncol(dist_mic_gibbon)) {
    dist_mic_gibbon[i, j] <- sqrt(sum((mic_coords[i, ] - gibbon_group_coords[j, ])^2))
    bearing_mic_gibbon[i, j] <- atan((mic_coords[i, 2] - gibbon_group_coords[j, 2]) / (mic_coords[i, 1] - gibbon_group_coords[j, 1]))
  }
}

# half normal hazard distribution for probability matrix
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
plot_select = 20

for (i in plot_select) {
  detection_mic_idx = which(detection_matrix[,i] == 1)
  print(detection_mic_idx)
  plot(x_gibbon_group, y_gibbon_group, type = "p", col = "blue", pch = 1, xlab = "X-coordinate (m)", ylab = "Y-coordinate (m)")
  points(mic_coords[,1], mic_coords[,2], type = "p", col = "black", pch = 15)
  points(x_gibbon_group[i], y_gibbon_group[i], type = "p", col = "red", pch = 15)
  points(mic_coords[detection_mic_idx,1], mic_coords[detection_mic_idx,2], type = "p", col = "green", pch = 15)
}

### create dataframes for analysis
mic_df = cbind(1:nrow(mic_coords), mic_coords)
colnames(mic_df)[1] = "mic_id"

# these are ground truth locations of each gibbon group and exact call times
gibbon_group_df = cbind(1:nrow(gibbon_group_coords), gibbon_group_coords)
colnames(gibbon_group_df)[1] = "gibbon_group_id"

# generate random ground truth call times
# todo: generate calls sequentially as events are dependent
# todo: sort out datetime zones
dawn_chorus_start_time = as.POSIXct("2023-01-01 04:00:00", tz = "UTC")
dawn_chorus_end_time = as.POSIXct("2023-01-01 05:00:00", tz = "UTC")

call_datetimes = as.POSIXct(runif(nrow(gibbon_group_df), dawn_chorus_start_time, dawn_chorus_end_time), origin = "1970-01-01", tz = "UTC")
call_timestamps = as.integer(call_datetimes)

gibbon_group_df = cbind(gibbon_group_df, call_datetimes)
colnames(gibbon_group_df)[4] = "call_datetime"

# speed of sound (m/s)
speed_of_sound = 331
von_mises_kappa = 4 # modelling parameter

recording_df = data.frame()

for (i in 1:nrow(mic_df)) {
  # construct detection dataframe for each mic
  recording_temp = data.frame()
  ground_truth_animal_ID = which(detection_matrix[i,] == 1)
  
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
  recording_temp = cbind(rep(i, times = length(ground_truth_animal_ID)), ground_truth_animal_ID, ground_truth_call_datetime, measured_call_datetime, ground_truth_bearing, measured_bearing)
  recording_df = rbind(recording_df, recording_temp)
}

