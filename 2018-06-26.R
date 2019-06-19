library(vegan)
library(lme4)
library(tidyverse)
library(gridExtra)

options(tibble.print_max = Inf)

befta.all <- read_csv("/Users/grahamprescott/Documents/1. Research/BEFTA birds/befta.csv")

# remove unidentified records, calculate total abundance, and number of species for each point count
befta.all %>% 
  filter(common != "unidentified") %>% 
  group_by(plot, stage, trial) %>% 
  summarise(species = n_distinct(common),
             abundance = length(common)) %>%
  mutate(treatment = ifelse(plot == "C10"|plot=="C19"|plot=="D29"|plot=="F09"|plot=="G09"|plot=="G14", 
                            "Reduced", ifelse(plot == "C11"|plot=="C17"|plot=="D28"|plot=="F04"|plot=="G07"|plot=="G16", "Enhanced", "Normal"))) -> befta.summary

# average across point counts

befta.summary %>%
  group_by(plot, stage) %>%
  summarise(mean.species = mean(species),
            mean.abundance = mean(abundance),
            sd.species = sd(species),
            sd.abundance = sd(abundance)) %>%
  mutate(treatment = ifelse(plot == "C10"|plot=="C19"|plot=="D29"|plot=="F09"|plot=="G09"|plot=="G14", 
                            "Reduced", ifelse(plot == "C11"|plot=="C17"|plot=="D28"|plot=="F04"|plot=="G07"|plot=="G16", "Enhanced", "Normal"))) -> befta.plot.means

ggplot(subset(befta.summary, stage == "Before")) + geom_boxplot(aes(x = treatment, y = species)) +
  theme_classic(18) + labs(title = "Before", y = "Species richness", x = "") + 
  theme(plot.title = element_text(hjust = 0.5)) -> figure.1a

ggplot(subset(befta.summary, stage == "After")) + geom_boxplot(aes(x = treatment, y = species)) +
  theme_classic(18) + labs(title = "After", y = "", x = "") + 
  theme(plot.title = element_text(hjust = 0.5)) -> figure.1b

ggplot(subset(befta.summary, stage == "Before")) + geom_boxplot(aes(x = treatment, y = abundance)) +
  theme_classic(18) + labs(title = "Before", y = "Abundance", x = "") + 
  theme(plot.title = element_text(hjust = 0.5)) -> figure.2a

ggplot(subset(befta.summary, stage == "After")) + geom_boxplot(aes(x = treatment, y = abundance)) +
  theme_classic(18) + labs(title = "After", y = "", x = "") + 
  theme(plot.title = element_text(hjust = 0.5)) -> figure.2b

grid.arrange(figure.1a, figure.1b, ncol = 2)
grid.arrange(figure.2a, figure.2b, ncol = 2)

richness.before.model <- aov(data = subset(befta.summary, stage == "Before"), species ~ treatment)
summary(richness.before.model) # F = 0.755, p = 0.474
plot(richness.before.model) # diagnostic plots look ok

richness.after.model <- aov(data = subset(befta.summary, stage == "After"), species ~ treatment)
summary(richness.after.model) # F = 2.288, p = 0.118
plot(richness.after.model) # diagnositc plots look ok

abundance.before.model <- aov(data = subset(befta.summary, stage == "Before"), abundance ~ treatment)
summary(abundance.before.model) # F = 1.736, p = 0.184
plot(abundance.before.model) # normal Q-Q plot a bit off

abundance.after.model <- aov(data = subset(befta.summary, stage == "After"), species ~ treatment)
summary(abundance.after.model) # F = 2.28, p = 0.118
plot(abundance.after.model) # diagnostics are fine
