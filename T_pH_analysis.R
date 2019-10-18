# libraries
library(tidyverse)
library(readxl)
library(ggpubr)

# load file
species = read_xlsx('species_list_final.xlsx', sheet = 'Sheet1') %>%
    mutate(Mutation = as.factor(Mutation))

# get summary statistics
species %>% group_by(Mutation) %>%
	summarise(Mean = mean(Temperature, na.rm = T),
			  SD = sd(Temperature, na.rm = T))

# plot Temperature vs mutation with ANOVA test
species %>%
	filter(Temperature != is.na(Temperature)) %>%
	ggplot(aes(x = Mutation, y = Temperature)) +
	geom_boxplot(aes(fill = Mutation)) +
	geom_jitter(shape = 16, position = position_jitter(0.2)) +
	theme_classic() +
	stat_compare_means(method = 'aov', label.x = 1.5, label.y = 110)


quartz.save(file = 'temperature_boxplot.pdf',
    type = 'pdf', dpi = 300, height = 6, width = 6)


# plot pH vs mutation with ANOVA test
species %>%
	filter(pH != is.na(pH)) %>%
	ggplot(aes(x = Mutation, y = pH)) +
	geom_boxplot(aes(fill = Mutation)) +
	geom_jitter(shape = 16, position = position_jitter(0.2)) +
	theme_classic() +
	stat_compare_means(method = 'aov', label.x = 1.5)


quartz.save(file = 'pH_boxplot.pdf',
    type = 'pdf', dpi = 300, height = 6, width = 6)






