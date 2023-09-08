set.seed(123)

# x and y in square survey site with dimensions 10kmx10km with markers every 2000m
# mics
x_mic = seq(from = 2000, to = 8000, by = 2000)
y_mic = seq(from = 2000, to = 8000, by = 2000)

mic_coords = expand.grid(x = x_mic, y = y_mic)
colnames(mic_coords) = c('x_mic', 'y_mic')

# expectation of gibbon groups in survey area
intensity = 0.000001
x_perimeter = 10000
y_perimeter = 10000
survey_area = x_perimeter * y_perimeter
expected_gibbon_groups = intensity * survey_area

# homogeneous Poisson point process for getting gibbon coordinates
x_gibbon_group <- runif(expected_gibbon_groups, 0, x_perimeter)
y_gibbon_group <- runif(expected_gibbon_groups, 0, y_perimeter)

gibbon_group_coords = cbind(x_gibbon_group, y_gibbon_group)

# TODO: this plot looks like ass, please fix
plot(x_gibbon_group, y_gibbon_group)



