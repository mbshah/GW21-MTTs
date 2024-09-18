# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Read the Excel file
file_path <- "/home/hb0358/Disk2/Disk2/Data/MTT_mbs/Statistics_GW21/SNR_Analysis_2_all.xlsx"
data <- read_excel(file_path)

# Separate the significance columns
significance_cols <- data %>%
  dplyr::select(Time, Pair_n, ends_with("Significance"))

# Reshape significance columns to long format
long_significance <- significance_cols %>%
  pivot_longer(
    cols = -c(Time, Pair_n),
    names_to = "Measure",
    names_pattern = "(.*)_Significance",
    values_to = "Significance"
  )

# Reshape the main data to long format
long_data <- data %>%
  dplyr::select(Time, Pair_n, ends_with("mean"), ends_with("SD"), ends_with("meanD"), ends_with("Factor")) %>%
  pivot_longer(
    cols = -c(Time, Pair_n),
    names_to = c("Measure", "Type"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  pivot_wider(
    names_from = Type,
    values_from = value
  )

# Merge the significance data back into the main long data
long_data <- long_data %>%
  left_join(long_significance, by = c("Time", "Pair_n", "Measure"))

long_data <- long_data %>%
  mutate(Technology = case_when(
    Measure %in% c("MTTO", "MTTP") ~ "Metatranscriptomics",
    Measure %in% c("18S", "16S") ~ "Metabarcoding",
    Measure == "NTS" ~ "Non Target Analysis",
    TRUE ~ NA_character_  # This will handle any other cases, if needed
  ))

# Remove skewed time points
long_data <- long_data %>%
  mutate(mean = ifelse(Measure == "NTS" & Time == 96, NA, mean)) %>%
  mutate(meanD = ifelse(Measure == "NTS" & Time == 96, NA, meanD)) %>%
  mutate(SD = ifelse(Measure == "NTS" & Time == 96, NA, SD)) %>%
  mutate(Factor = ifelse(Measure == "NTS" & Time == 96, NA, Factor)) %>%
  mutate(Significance = ifelse(Measure == "NTS" & Time == 96, NA, Significance)) %>%
  mutate(Technology = ifelse(Measure == "NTS" & Time == 96, NA, Technology))

long_data <- long_data %>%
  mutate(mean = ifelse(Measure == "MTTO" & Time == 12, NA, mean)) %>%
  mutate(meanD = ifelse(Measure == "MTTO" & Time == 12, NA, meanD)) %>%
  mutate(SD = ifelse(Measure == "MTTO" & Time == 12, NA, SD)) %>%
  mutate(Factor = ifelse(Measure == "MTTO" & Time == 12, NA, Factor)) %>%
  mutate(Significance = ifelse(Measure == "MTTO" & Time == 12, NA, Significance)) %>%
  mutate(Technology = ifelse(Measure == "MTTO" & Time == 12, NA, Technology))

#long_data <- long_data %>% mutate(adjusted_time = case_when(Time == 1 ~ '18S.1', Time == 12 ~ '18S.12', Time == 24 ~ '18S.24', Time == 48 ~ '18S.48', Time == 96 ~ '18S.96', Time == 168 ~ '18S.168', Time == 336 ~ '18S.336'))

# Create the plot with significance markers
box_data <- data.frame(
  xmin = c(0.5, 5.5, 10.5, 15.5, 20.5, 25.5, 30.5),
  xmax = c(5.5, 10.5, 15.5, 20.5, 25.5, 30.5, 35.5),
  ymin = -Inf,
  ymax = Inf
)
time_label_data <- data.frame(
  x = c(3, 8, 13, 18, 23, 28, 33),
  y = 0.8,
  labels = c("1 hour", "12 hours", "24 hours", "48 hours", "96 hours", "168 hours", "240 hours")
)
bar_plot <- ggplot(long_data, aes(x = interaction(Measure, Time), y = mean, fill = Pair_n, group = Measure)) +
  geom_bar(stat = "identity", color = "white") +
  geom_errorbar(aes(ymin = meanD - SD, ymax = meanD + SD, group = Pair_n),
                position = position_dodge(),
                stat = "identity", width = 0.5, size = 0.5) +
  geom_text(aes(label = paste0(ifelse(is.na(Significance), "", Significance),"\n", ifelse(is.na(Factor), "", Factor)), y = 0),
            position = position_dodge(0.9), size = 3, vjust = 0.5) +
  geom_text(data = time_label_data, aes(x = x, y = y, label = labels), color = "black", size = 5, inherit.aes = FALSE) +
  geom_rect(data = box_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", inherit.aes = FALSE) +
  labs(x = "Time Points", y = "Bray-Curtis Dissimilarity", fill = "Pair") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "B: Plot with all SNR values significance and Dissimilairty")

bar_plot2 <- ggplot(long_data, aes(x = as.factor(Time), y = mean, fill = Pair_n)) +
  geom_bar(stat = "identity", color = "white",position = position_stack(reverse = TRUE)) +
  geom_errorbar(aes(ymin = meanD - SD, ymax = meanD + SD),
                position = position_dodge(),
                width = 0.5, size = 0.5) +
  geom_text(aes(label = paste0(ifelse(is.na(Significance), "", Significance),"\n", ifelse(is.na(Factor), "", Factor)), y = 0),
            position = position_dodge(0.9), size = 3, vjust = 0.5) +
  facet_wrap(~ Measure, scales = "free_y", nrow = 1) +
  labs(x = "Time Points", y = "Bray-Curtis Dissimilarity", fill = "Pair") +theme_minimal()+
  theme(
    strip.background = element_rect(color = "black", fill = "white", size = 1, linetype = "solid"),
    strip.text = element_text(face = "bold"),axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  ) +
  labs(title = "B: Plot with all SNR values significance and Dissimilarity")

# Create the line plot without the "Pair" legend
line_plot <- ggplot(long_data, aes(x = Time, y = mean, fill = Pair_n, group = interaction(Measure, Pair_n))) +
  geom_line(data = long_data[!is.na(long_data$Factor),], aes(x = Time, y = Factor, color = Measure, group = Measure), size = 1) +
  geom_point(data = long_data[!is.na(long_data$Factor),], aes(x = Time, y = Factor,  shape= Measure), size = 3) +
  scale_x_continuous(breaks = unique(long_data$Time),
                     labels = as.character(unique(long_data$Time))) +
  labs(x = "Time Points", y = "Signal to Noise Ratios") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "A: Plot to visualise the SNR trends") +
  guides(fill = "none")

line_plot2 <- ggplot(long_data, aes(x = Time, y = mean, fill = Pair_n, group = interaction(Measure, Pair_n))) +
  #geom_line(data = long_data[!is.na(long_data$Factor),], aes(x = Time, y = Factor, color = Measure), size = 1, na.rm = FALSE) +
  #geom_line(data = long_data %>% filter(Measure == "MTTO" & (Time == 1 | Time == 24)),
  #          aes(x = Time, y = Factor), color = "red", size = 1, linetype = "dotted")+
  geom_point(data = long_data[!is.na(long_data$Factor),], aes(x = Time, y = Factor,  shape= Technology, colour = Measure), size = 3) +
  scale_x_continuous(breaks = unique(long_data$Time),
                     labels = as.character(unique(long_data$Time))) +
  labs(x = "Time Points", y = "Signal to Noise Ratios") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "A: Plot to visualise the SNR trends") +
  guides(fill = "none")


# Define custom colors and shapes
custom_colors <- c("18S" = "blue", "16S" = "darkblue",
                   "MTTP" = "lightcoral", "MTTO" = "darkred",
                   "NTS" = "black")

custom_shapes <- c("Metabarcoding" = 19, #"Metabarcoding" = 1,  # Empty circle for 18S, filled circle for 16S
                   "Metatranscriptomics" = 17, #"Metatranscriptomics" = 2, # Empty triangle for MTTP, filled triangle for MTTO
                   "Non Target Analysis" = 15)              # Filled square for NTS

# Create the line plot with custom colors and shapes
line_plot3 <- ggplot(long_data, aes(x = Time, y = mean, fill = Pair_n, group = interaction(Measure, Pair_n))) +
  # Solid lines for all continuous data
  geom_line(data = long_data[(long_data$Measure == "MTTO") & !(long_data$Time %in% c(12,1)), ],
            aes(x = Time, y = Factor, color = Measure), size = 0.4) +
  geom_point(data = long_data[(long_data$Measure == "MTTO"), ],
             aes(x = Time, y = Factor, shape = Technology, color = Measure), size = 3, shape = 17) +
  geom_line(data = long_data[(long_data$Measure == "MTTP"), ],
            aes(x = Time, y = Factor, color = Measure), size = 0.4) +
  geom_point(data = long_data[(long_data$Measure == "MTTP"), ],
             aes(x = Time, y = Factor, shape = Technology, color = Measure), size = 3, shape = 2) +
  geom_line(data = long_data[(long_data$Measure == "18S"), ],
            aes(x = Time, y = Factor, color = Measure), size = 0.4) +
  geom_point(data = long_data[(long_data$Measure == "18S"), ],
             aes(x = Time, y = Factor, shape = Technology, color = Measure), size = 3, shape = 1) +
  geom_line(data = long_data[(long_data$Measure == "16S") , ],
            aes(x = Time, y = Factor, color = Measure), size = 0.4) +
  geom_point(data = long_data[(long_data$Measure == "16S"), ],
             aes(x = Time, y = Factor, shape = Technology, color = Measure), size = 3, shape = 19) +
  geom_line(data = long_data[(long_data$Measure == "NTS") & (long_data$Time >= 168), ],
            aes(x = Time, y = Factor, color = Measure), size = 0.4) +
  geom_line(data = long_data[(long_data$Measure == "NTS") & (long_data$Time <= 48), ],
            aes(x = Time, y = Factor, color = Measure), size = 0.4) +
  geom_point(data = long_data[(long_data$Measure == "NTS"), ],
             aes(x = Time, y = Factor, shape = Technology, color = Measure), size = 3, shape = 15) +
  # Dotted lines
  geom_line(data = long_data %>% filter(Measure == "MTTO" & (Time == 1 | Time == 24)),
            aes(x = Time, y = Factor, color = Measure), size = 0.4, linetype = "dotted") +
  geom_line(data = long_data %>% filter(Measure == "NTS" & (Time == 48 | Time == 168)),
            aes(x = Time, y = Factor, color = Measure), size = 0.4, linetype = "dotted") +
  #geom_point(data = long_data[!is.na(long_data$Factor),],
  #           aes(x = Time, y = Factor, shape = Technology, color = Measure), size = 3) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = custom_shapes) +
  scale_x_continuous(breaks = unique(long_data$Time),
                     labels = as.character(unique(long_data$Time))) +
  labs(x = "Time Points", y = "Signal to Noise Ratios") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title = "A: Plot to visualise the SNR trends") +
  guides(fill = "none")

# Combine the plots using patchwork
combined_plot1 <- line_plot3 / bar_plot2

# Display the combined plot
print(combined_plot1)
