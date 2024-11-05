setwd("C:/Users/uqqzhao/UQ/1projects/IMB_group_Palpant/202306_NatureProtocols/202409_BriefingsInBioinformatics/revision/Scalability") 

library(dplyr)
library(ggplot2)

memory_data <- read.table("TRIAGEcluster_memory_results.txt", header = TRUE, sep = "\t")
memory_data <- memory_data %>%
  mutate(Max_memory = case_when(
    grepl("MB", Max_memory) ~ as.numeric(sub(" MB", "", Max_memory)) / 1024,
    grepl("GB", Max_memory) ~ as.numeric(sub(" GB", "", Max_memory)),
    TRUE ~ as.numeric(Max_memory)
  ))

runtime_data <- read.table("TRIAGEcluster_runtime_results.txt", header = TRUE, sep = "\t")


memory_summary <- memory_data %>%
  group_by(Cell_number) %>%
  summarise(
    Mean_memory = mean(Max_memory, na.rm = TRUE),
    SD_memory = sd(Max_memory, na.rm = TRUE)
  )

runtime_summary <- runtime_data %>%
  group_by(Cell_number) %>%
  summarise(
    Mean_runtime = mean(Runtime..mins., na.rm = TRUE),
    SD_runtime = sd(Runtime..mins., na.rm = TRUE)
  )


combined_data <- merge(memory_summary, runtime_summary, by = "Cell_number")


ggplot(combined_data, aes(x = Cell_number)) +
  # Runtime
  geom_line(aes(y = Mean_runtime), color = "#1f78b4", size = 1, linetype = "solid") +
  geom_point(aes(y = Mean_runtime), color = "#1f78b4", size = 2) +
  geom_errorbar(aes(ymin = Mean_runtime - SD_runtime, ymax = Mean_runtime + SD_runtime), 
                width = 0.5, color = "#1f78b4") +
  
  geom_ribbon(aes(ymin = Mean_runtime - SD_runtime, ymax = Mean_runtime + SD_runtime), 
              fill = "#1f78b4", alpha = 0.2) +
  
  # Memory
  geom_line(aes(y = Mean_memory), color = "#ff7f00", size = 1, linetype = "dashed") +
  geom_point(aes(y = Mean_memory), color = "#ff7f00", size = 2) +
  geom_errorbar(aes(ymin = Mean_memory - SD_memory, ymax = Mean_memory + SD_memory), 
                width = 0.5, color = "#ff7f00") +
  
  geom_ribbon(aes(ymin = Mean_memory - SD_memory, ymax = Mean_memory + SD_memory), 
              fill = "#ff7f00", alpha = 0.2) +
  
  scale_y_continuous(
    name = "Runtime (mins)", limits = c(0, 60),
    sec.axis = sec_axis(~ ., name = "Max memory (GB)", breaks = seq(0, 60, by = 10))
  ) +
  
  scale_x_continuous(breaks = c(20000, 40000, 60000, 80000, 100000)) +
  labs(x = "The number of cells") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(color = "#1f78b4"),
    axis.title.y.right = element_text(color = "#ff7f00")
  )

ggsave("TRIAGEcluster_scalability.pdf", width = 6, height = 6)

