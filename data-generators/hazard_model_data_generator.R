set.seed(123)

# x and y in square survey site with dimensions 10km x 10km with markers every 2000m
# mics
x_mic = seq(from = 2000, to = 8000, by = 2000)
y_mic = seq(from = 2000, to = 8000, by = 2000)

mic_coords = expand.grid(x = x_mic, y = y_mic)
colnames(mic_coords) = c('x_mic', 'y_mic')
mic_coords = as.matrix(mic_coords)

# expectation of gibbon groups in survey area
intensity = 0.000001
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
points(mic_coords[,1], mic_coords[,2], type = "p", col = "red", pch = 15)
}

# euclidean distance for detector to gibbon observation
dist_mic_gibbon <- matrix(NA, nrow = nrow(mic_coords), ncol = nrow(gibbon_group_coords))

for (i in 1:nrow(dist_mic_gibbon)) {
  for (j in 1:ncol(dist_mic_gibbon)) {
    dist_mic_gibbon[i, j] <- sqrt(sum((mic_coords[i, ] - gibbon_group_coords[j, ])^2))
  }
}

# half normal hazard distribution for probability matrix
lambda = 50
sigma = 500

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

# plot detected gibbon groups for chosen mic
plot_select = 9

for (i in plot_select) {
  detection_gibbon_group_idx = which(detection_matrix[i,] == 1)
  print(detection_gibbon_group_idx)
  plot(x_gibbon_group, y_gibbon_group, type = "p", col = "blue", pch = 1, xlab = "X-coordinate (m)", ylab = "Y-coordinate (m)")
  points(mic_coords[,1], mic_coords[,2], type = "p", col = "red", pch = 15)
  points(x_gibbon_group[detection_gibbon_group_idx], y_gibbon_group[detection_gibbon_group_idx], type = "p", col = "green", pch = 15)
}

