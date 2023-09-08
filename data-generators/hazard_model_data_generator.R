set.seed(123)

# x and y in square survey site with dimensions 10kmx10km with markers every 2000m
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

# half normal hazard distribution
lambda = 500
sigma = 50









