setwd("C:/Users/uqqzhao/UQ/1projects/IMB_group_Palpant/202306_NatureProtocols/202409_BriefingsInBioinformatics/revision/Scalability") 
library(ggplot2)
library(dplyr)
library(readr)

memory_data <- read.table("TRIAGEgene_memory_results.txt", header = TRUE, sep = "\t")
runtime_data <- read.table("TRIAGEgene_runtime_results.txt", header = TRUE, sep = "\t")

# MB -> GB
memory_data <- memory_data %>%
  mutate(Max_memory = case_when(
    grepl("MB", Max_memory) ~ as.numeric(sub(" MB", "", Max_memory)) / 1024,
    grepl("GB", Max_memory) ~ as.numeric(sub(" GB", "", Max_memory)),
    TRUE ~ as.numeric(Max_memory)
  ))

#unique(memory_data$Max_memory)
#memory_data

# mean and sd
memory_summary <- memory_data %>%
  group_by(Sample_size) %>%
  summarise(
    Mean_memory = mean(Max_memory, na.rm = TRUE),
    SD_memory = sd(Max_memory, na.rm = TRUE)
  )

runtime_summary <- runtime_data %>%
  group_by(Sample_size) %>%
  summarise(
    Mean_runtime = mean(Runtime..mins., na.rm = TRUE),
    SD_runtime = sd(Runtime..mins., na.rm = TRUE)
  )

combined_data <- merge(memory_summary, runtime_summary, by = "Sample_size")

ggplot(combined_data, aes(x = Sample_size)) +
  # Y-axis leftside
  geom_line(aes(y = Mean_runtime), color = "#1f78b4", size = 1, linetype = "solid") +
  geom_point(aes(y = Mean_runtime), color = "#1f78b4", size = 2) +
  geom_errorbar(aes(ymin = Mean_runtime - SD_runtime, ymax = Mean_runtime + SD_runtime), 
                width = 0.5, color = "#1f78b4") +
  
  geom_ribbon(aes(ymin = Mean_runtime - SD_runtime, ymax = Mean_runtime + SD_runtime), 
              fill = "#1f78b4", alpha = 0.2) +
  
  # Y-axis rightside
  geom_line(aes(y = Mean_memory * 6.67), color = "#ff7f00", size = 1, linetype = "dashed") +
  geom_point(aes(y = Mean_memory * 6.67), color = "#ff7f00", size = 2) +
  geom_errorbar(aes(ymin = (Mean_memory - SD_memory) * 6.67, ymax = (Mean_memory + SD_memory) * 6.67), 
                width = 0.5, color = "#ff7f00") +
  
  geom_ribbon(aes(ymin = (Mean_memory - SD_memory) * 6.67, ymax = (Mean_memory + SD_memory) * 6.67), 
              fill = "#ff7f00", alpha = 0.2) +
  
  scale_y_continuous(
    name = "Runtime (mins)", limits = c(0, 80),
    sec.axis = sec_axis(~ . / 6.67, name = "Max memory (GB)",
                        breaks = c(0, 3, 6, 9, 12), labels = c(0, 3, 6, 9, 12))
  ) +
  
  labs(x = "Sample size") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(color = "#1f78b4"),
    axis.title.y.right = element_text(color = "#ff7f00")
  )

ggsave("TRIAGEgene_scalability.pdf", width = 6, height = 6)

