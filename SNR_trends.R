# Load necessary libraries
library(ggplot2)
library(tidyr)
library(stringr)
library(RColorBrewer)

# Define the data directly from the CSV file
folder_path <- "/home/hb0358/Disk2/Disk2/Data/MTT_mbs/Statistics_GW21/"
data <- as.data.frame(read.csv(paste0(folder_path, "SNR_pathways_raw.csv"), check.names = FALSE))

# Reshape the data into long format
data_long <- gather(data, key = "Time", value = "Value", -Tech)

# Convert the Time column to numeric hours using str_replace_all from stringr
data_long$Time <- as.numeric(str_replace_all(data_long$Time, c(" Hours" = "", " Hour" = "")))

# Filter the data for the desired technologies
plot_data <- data_long[data_long$Tech %in% c("18S", "16S", "NTS", "MTT_Pr_PW", "MTT_PrS_PW", "MTT_Pr_KO", "MTT_PrS_KO", "MTT_PrS_KO_NOC1_S2"), ]

# Create a new column for line size
plot_data$LineSize <- ifelse(plot_data$Tech %in% c("x"), 2, 1)

# Define a colorblind-friendly palette
tech_colors <- brewer.pal(n = length(unique(plot_data$Tech)), name = "Set2")
names(tech_colors) <- unique(plot_data$Tech)

# Add the Color column to the plot_data data frame
plot_data$Color <- tech_colors[match(plot_data$Tech, names(tech_colors))]

# Plot the data with proportional spacing on the x-axis
ggplot(plot_data, aes(x = Time, y = Value, color = Tech, group = Tech)) +
  geom_line(aes(size = LineSize)) +
  #geom_point(aes(shape = Tech)) +
  scale_x_continuous(breaks = c(1, 12, 24, 48, 96, 168, 336), labels = c("1 Hour", "12 Hours", "24 Hours", "48 Hours", "96 Hours", "168 Hours", "336 Hours")) +
  scale_size_identity() +
  scale_color_manual(values = tech_colors) +
  theme_minimal() +
  labs(title = "Trends of Different Technologies Over Time",
       x = "Time",
       y = "Value",
       color = "Technology")
