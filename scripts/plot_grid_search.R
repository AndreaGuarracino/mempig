library(tidyverse)

# Read the data
data <- read.table("/home/guarracino/Desktop/C4.test.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char = '?')

results <- data %>%
  mutate(status = map2_lgl(sample.id, 
                           str_c(haplotype.1, "|", haplotype.2),
                           ~{
                             # Check if sample.id is in both haplotypes
                             has_sample_id <- str_detect(str_c(.y), .x)
                             
                             # If sample_id present, check #1/#2 pattern
                             if(has_sample_id) {
                               haps <- str_split(.y, "\\|")[[1]]
                               # Check if haplotype.1 has #1 and haplotype.2 has #2
                               # OR haplotype.1 has #2 and haplotype.2 has #1
                               (str_detect(haps[1], "#1") && str_detect(haps[2], "#2")) ||
                                 (str_detect(haps[1], "#2") && str_detect(haps[2], "#1"))
                             } else {
                               FALSE
                             }
                           })) %>%
  mutate(status = if_else(status, "okay", "not_okay"))

# Then group and summarize
summary <- results %>%
  filter (sample.id == 'HG02080') %>%
  group_by(l, p) %>%
  summarize(
    okay_count = sum(status == "okay"),
    not_okay_count = sum(status == "not_okay"),
    .groups = "drop"
  )

# Calculate ratios
plot_data <- summary %>%
  mutate(
    total = okay_count + not_okay_count,
    ratio = not_okay_count / total
  )

# Get the actual maximum ratio for scaling
max_ratio <- max(plot_data$ratio)
min_ratio <- min(plot_data$ratio)

# Create the heatmap
p <- ggplot(plot_data, aes(x = factor(l), y = factor(p))) +
  geom_tile(aes(fill = ratio)) +
  # Second layer: green tiles for zero not_okay counts
  geom_tile(data = plot_data %>% filter(not_okay_count == 0),
            fill = "#99FF99") +
  # Add text layer for not_okay counts
  geom_text(aes(label = not_okay_count, 
                # Make text white if background is dark, black if background is light
                color = if_else(ratio > mean(c(min_ratio, max_ratio)), "white", "black")),
            size = 2.5) +
  scale_fill_gradient(
    high = "black",    # minimum ratio in data
    low = "white",   # maximum ratio in data
    limits = c(min_ratio, max_ratio)
  ) +
  scale_color_identity() +  # Use the colors directly as specified in the color aesthetic
  labs(
    #title = paste("mempig - C4 grid search"),
    x = "min-MEM-length",
    y = "max-num-pangenome-MEMs",
    fill = "nok / (nok +ok)"
  ) +
  theme_minimal() +
  coord_equal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 6),
    legend.position = "bottom",
    plot.title = element_text(margin = margin(b = 0), hjust = 0.5),
  )
p
# Save the plot
ggsave(
  "C4.grid-search.heatmap.HG02080.png",
  plot = p,
  width = 17,
  height = 8,
  dpi = 300,
  bg = "white"  # ensures white background
)



p2 <- ggplot(results %>% filter(l > 9 & l < 15 & p < 25), 
       aes(x = sample.id, y = 1, fill = status, width = 0.8)) + # Added width parameter
  geom_tile() +
  facet_grid(fct_rev(factor(p)) ~ l, scales = "free_y") + # Added fct_rev here
  scale_fill_manual(values = c("okay" = "darkgrey", "not_okay" = "firebrick")) +
  scale_y_continuous(breaks = c(1), labels = c(""), limits = c(0.5, 1.5)) +
  labs(
    x = "Sample ID",
    y = "",
    fill = "Status"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    strip.text = element_text(size = 10),
    panel.spacing = unit(0.25, "lines"),
    panel.grid = element_blank(),
    legend.position = "bottom",
  )

ggsave(
  "C4.grid-search.by_sample.subset.png",
  plot = p2,
  width = 17,
  height = 8,
  dpi = 300,
  bg = "white"  # ensures white background
)


# Create summary per sample
sample_summary <- results %>%
  group_by(sample.id) %>%
  summarize(
    okay_count = sum(status == "okay"),
    not_okay_count = sum(status == "not_okay"),
    total = n(),
    okay_ratio = okay_count / total
  ) %>%
  arrange(okay_ratio)  # Sort by ratio of okay results

# Create the third plot
p3 <- ggplot(sample_summary, 
             aes(x = reorder(sample.id, okay_ratio))) +
  geom_col(aes(y = not_okay_count), fill = "firebrick", alpha = 0.8) +
  geom_col(aes(y = okay_count), fill = "forestgreen", alpha = 0.8) +
  labs(
    x = "Sample ID",
    y = "Count of parameter pairs (l,p)",
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 11),
    plot.title = element_text(margin = margin(b = 0), hjust = 0.5),
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = function(x) pretty(x, n = 20)  # Increase number of breaks
  )
p3

ggsave(
  "C4.grid-search.by_sample.ok-nok.png",
  plot = p3,
  width = 17,
  height = 8,
  dpi = 300,
  bg = "white"  # ensures white background
)

