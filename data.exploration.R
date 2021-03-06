# test commit to see if GitHub connectio working
# for ggplot2, dplyr, magritr etc
library(tidyverse)
# for multivariate analysis and species accumulation curves
library(vegan)
# for mixed effects models
libary(lme4)

library(gridExtra)

# load data
birds.full <- read_csv("./BEFTA birds 2013-2019 final 12.06.2020.csv")

# identify structure
str(birds.full)
# 1,547 records

# filter full dataset 

birds.full %>%
  mutate(common = str_replace_all(common, c("Cinereous Tit" = "Japanese Tit",
                                            "Plaintive cuckoo" = "Plaintive Cuckoo"))) %>%
  mutate(latin = str_replace_all(latin, c("Parus cinereus" = "Parus minor"))) %>%
  subset(uncertain == 0) %>%
  group_by(family, latin, common) %>%
  count(sort=TRUE) %>%
  #head()
  #write.csv("./TableS1.csv")
  view()

birds.full %>%
  subset(uncertain != 0) %>%
  view()
# filter dataset
birds <- birds.full %>%
  # remove flyovers and observations > 70 m [to only include records in core plot]
  subset(distance < 30) %>%
  # remove records of uncertain ID [doubtful records, those marked 'Unidentified', Sunda frogmouth etc had value of 1 in 'uncertain' column]
  subset(uncertain == 0) %>%
  # remove red junglefowl / domestic chicken
  subset(common != "Red Junglefowl") %>%
  #remove Eurasian tree sparrow
  subset(common != "Eurasian Tree Sparrow") %>%
  # remove Oriental Mapgie Robin (poaching pressure biases results)
  subset(common != "Oriental Magpie-Robin") %>%
  # remove Sunda Frogmout (nocturnal birds not adequately sampled)
  subset(common != "Sunda Frogmouth")
  # merge the bulbuls together into one species; rename Cinereous to Japanese Tit (following recent paper)
  #mutate(common = str_replace_all(common, c("Olive-winged Bulbul" = "Bulbul","Yellow-vented Bulbul" = "Bulbul", "Sooty-headed Bulbul" = "Bulbul","Cinereous Tit" = "Japanese Tit"))) %>%
  # adjust Latin names too 
  #mutate(latin = str_replace_all(latin, c("Pycnonotus plumosus" = "Pycnonotus sp.","Pycnonotus goiavier" = "Pycnonotus sp.","Pycnonotus aurigaster" = "Pycnonotus sp.","Parus cinereus" = "Parus minor"))) 
  

str(birds)
# 240 observations

# species by no. records
birds %>%
  group_by(latin, common) %>%
  count(sort=TRUE) %>%
  #write.csv("./Table1.csv")
  view()

# add mean counts per stage/treatment
birds %>%
  group_by(common, treatment) %>%
  count() %>%
  view()


# species by no. records, per treatment, before
birds %>%
  subset(stage == "Before") %>%
  group_by(common, treatment) %>%
  count() %>%
  view()

# species by no. records, per treatment, before
birds %>%
  subset(stage == "After") %>%
  group_by(common, treatment) %>%
  count() %>%
  view()

# at a glance some striking differences in abundance for ashy tailorbird, bar-winged prinia, yellow-vented bulbul

# graph this:
birds %>%
  mutate(stage = fct_relevel(stage, "Before", "After")) %>%
  group_by(latin, treatment, stage) %>%
  count() %>%
  ggplot() +
  geom_point(aes(x = latin, y = n, color = treatment, shape = treatment)) +
  scale_color_manual(values = c("dark green", "blue", "red")) +
  theme_bw(12) +
  coord_flip() + 
  facet_wrap(~stage) +
  labs(x = "", y = "total abundance", title = "Comparison of raw counts of individual species by treatment and stage")

# graph this:
birds %>%
  mutate(stage = fct_relevel(stage, "Before", "After")) %>%
  group_by(common, treatment, stage, round) %>%
  count() %>%
  ggplot() +
  geom_point(aes(x = common, y = n, color = treatment, shape = stage)) +
  coord_flip() + 
  facet_wrap(~round) +
  labs(x = "", y = "total abundance", title = "Comparison of raw counts of individual species by treatment and stage")

# calculate and display no individuals recorded per point count
birds %>%
  group_by(round, plot) %>%
  mutate(n = length(common)) %>%
  mutate(sr = n_distinct(common)) %>%
  ggplot() +
  geom_point(aes(x = plot, y = n, shape = stage, color = treatment)) +
  coord_flip() +
  facet_wrap(~round)  +
  labs( y = "number of individuals", x = "plot", title = "plot level variation in number of individuals per point count")

# calculate and display species richness recorded per point count
birds %>%
  group_by(round, plot) %>%
  mutate(n = length(common)) %>%
  mutate(sr = n_distinct(common)) %>%
  ggplot() +
  geom_point(aes(x = plot, y = sr, shape = stage, color = treatment)) +
  coord_flip() +
  facet_wrap(~round)  +
  labs( y = "number of species recorded", x = "plot", title = "plot level variation in observed species richness per point count")

# calculate and display overall differences in abundance and species richness by stage and treatment
birds %>%
  mutate(stage = fct_relevel(stage, "Before", "After")) %>%
  group_by(round, plot) %>%
  mutate(n = length(common)) %>%
  mutate(sr = n_distinct(common)) %>%
  ggplot() +
  geom_boxplot(aes(x = stage, y = sr, fill = stage)) +
  scale_fill_manual(values = c("white", "dark grey")) +
  theme_bw(12) +
  ylim(0,5) +
  facet_wrap(~treatment)  +
  labs( y = "number of species recorded", x = "", title = "")

# same but coloured by treatment
birds %>%
  mutate(stage = fct_relevel(stage, "Before", "After")) %>%
  group_by(round, plot) %>%
  mutate(n = length(common)) %>%
  mutate(sr = n_distinct(common)) %>%
  ggplot() +
  geom_boxplot(aes(x = stage, y = sr, fill = treatment)) +
  scale_fill_manual(values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme_bw(12) +
  ylim(0,5) +
  facet_wrap(~treatment)  +
  labs( y = "number of species recorded", x = "", title = "")

# calculate and display overall differences in abundance and species richness by stage and treatment
birds %>%
  mutate(stage = fct_relevel(stage, "Before", "After")) %>%
  group_by(round, plot) %>%
  mutate(n = length(common)) %>%
  mutate(sr = n_distinct(common)) %>%
  ggplot() +
  geom_boxplot(aes(x = stage, y = n, fill = stage)) +
  scale_fill_manual(values = c("white", "dark grey")) +
  theme_bw(12) +
  ylim(0,8) +
  facet_wrap(~treatment)  +
  labs( y = "number of individuals recorded", x = "", title = "")

# same but coloured by treatment
birds %>%
  mutate(stage = fct_relevel(stage, "Before", "After")) %>%
  group_by(round, plot) %>%
  mutate(n = length(common)) %>%
  mutate(sr = n_distinct(common)) %>%
  ggplot() +
  geom_boxplot(aes(x = stage, y = n, fill = treatment)) +
  scale_fill_manual(values = c("#1b9e77", "#7570b3", "#d95f02")) +
  theme_bw(12) +
  ylim(0,8) +
  facet_wrap(~treatment)  +
  labs( y = "number of individuals recorded", x = "", title = "")

# together (note, not sure how to integrate Shannon diversity)

birds %>%
  mutate(stage = fct_relevel(stage, "Before", "After")) %>%
  group_by(round, plot) %>%
  mutate(Abundance = length(common)) %>%
  mutate(Richness = n_distinct(common)) %>%
  #mutate(H = diversity(index = "Shannon")) %>%
  select(stage, treatment, plot, round, Abundance, Richness) %>%
  pivot_longer(-c(stage, treatment, plot, round), names_to = "measure", values_to = "value") %>%
  ggplot() +
  geom_boxplot(aes(x = treatment, y = value, fill = treatment)) +
  scale_fill_manual(values = c("dark green", "blue", "red")) +
  theme_bw(12) +
  ylim(0,8) +
  facet_grid(rows = vars(measure), cols = vars(stage))  +
  labs( y = "", x = "", title = "")

#ggsave(filename = paste0("richness.abundance.before.after.",format(Sys.time(), "%d-%b-%Y %H.%M"), ".png"), dpi = "print", units = "mm", device="png")

birds %>%
  select(common, plot, round) %>%
  unite("plot.round", plot:round) %>%
  group_by(plot.round, common) %>%
  summarise(count = length(common)) %>%
  head()
#  pivot_wider(names_from = common, values_from = count, values_fill = list(count = 0)) %>%
#  view()

# create new variable
wide.birds <- birds %>%
  select(common, plot, round) %>%
  unite("plot.round", plot:round) %>%
  group_by(plot.round, common) %>%
  summarise(count = length(common)) %>%
  pivot_wider(names_from = common, values_from = count, values_fill = list(count = 0))

# calculate species accumulation curves
plot(specaccum(wide.birds[,-1]))
specpool(wide.birds[,-1])

# prep for joining with plot data [handy for ordination graphs etc]
wide.birds <- wide.birds %>%
  separate(plot.round, c("plot", "round")) 

# prep plot data
plot.summary <- birds %>%
  select(round, plot, treatment, stage) %>%
  group_by(round, plot, treatment, stage) %>%
  tally()

# join together to create 'wider.birds' dataset for NMDS and species richness estimates by treatment, stage, etc. 
wider.birds <- merge(plot.summary, wide.birds)
#view(wider.birds)

# expected species richness after treatment - enhanced plots
wider.birds %>%
  subset(stage == "After") %>%
  subset(treatment == "enhanced") %>%
  select(-plot, -round, -n, -treatment) %>%
  specpool()

# expexted species richness after treatment - normal plots
wider.birds %>%
  subset(stage == "After") %>%
  subset(treatment == "normal") %>%
  select(-plot, -round, -n, -treatment) %>%
  specpool()

# expexted species richness after treatment - reduced plots
wider.birds %>%
  subset(stage == "After") %>%
  subset(treatment == "reduced") %>%
  select(-plot, -round, -n, -treatment) %>%
  specpool()

# my preliminary feeling is that we will see differences in community composition but not overall richness

# note that this is NMDS for everything together..
# NMDS on the data (excluding first 5 columns)
nmds.all <- metaMDS(wider.birds[,-c(1:5)], k = 3)
scores(nmds.all)

nmds.graph <- data.frame(scores(nmds.all), wider.birds[,c(1:5)])
#view(nmds.graph)

# display for all sites, before and after
nmds.graph %>%
  ggplot() +
  theme_classic(16) +
  geom_point(aes(x = NMDS1, y = NMDS2, shape = treatment, color = stage))

# display for all sites, after
nmds1b <- nmds.graph %>%
  subset(stage == "After") %>%
  ggplot() +
  theme_classic() +
  geom_point(aes(x = NMDS1, y = NMDS2, shape = treatment, color = treatment)) +
  labs( x = "NMDS1", y = "NMDS2", title = "b) 2018-9, after treatment")

# display for all sites, before 
nmds1a <- nmds.graph %>%
  subset(stage == "Before") %>%
  ggplot() +
  theme_classic() +
  theme(legend.position = "none") +
  geom_point(aes(x = NMDS1, y = NMDS2, shape = treatment, color = treatment)) +
  labs( x = "NMDS1", y = "NMDS2", title = "a) 2013, before treatment")

grid.arrange(nmds1a, nmds1b, ncol =2)

# but I can see that being an issue. What if we did the NMDS calculations separately?

head(wider.birds)

wider.birds.before <- wider.birds %>%
  subset(stage == "Before")

nmds.before <- metaMDS(vegdist(wider.birds.before[,-c(1:5)]), k = 2)
str(wider.birds.before[,c(1:5)])
str(scores(nmds.before))

figure.nmds.before <- data.frame(wider.birds.before[,c(1:5)], scores(nmds.before))
ggplot(figure.nmds.before) +
  theme_classic() +
  geom_point(aes(x = NMDS1, y = NMDS2, shape = treatment, color = treatment)) +
  scale_color_manual(values = c("dark green", "blue", "red")) +
  

#ggsave(filename = paste0("NMDS.before.only.", 
#                        format(Sys.time(), "%d-%b-%Y %H.%M"), ".png"), 
#       dpi = "print", units = "mm", device="png")

ggplot(figure.nmds.before) +
  theme_classic() +
  geom_text(aes(x = NMDS1, y = NMDS2, label = plot, color = treatment)) +
  scale_color_manual(values = c("dark green", "blue", "red")) 
  
#ggsave(filename = paste0("NMDS.before.only.w.plot.labels.", 
#                           format(Sys.time(), "%d-%b-%Y %H.%M"), ".png"), 
#         dpi = "print", units = "mm", device="png")

# now for the after plot:

head(wider.birds)

wider.birds.after <- wider.birds %>%
  subset(stage == "After")

nmds.after <- metaMDS(vegdist(wider.birds.after[,-c(1:5)]), k = 2)
# k=2, Run 20 stress 0.0001820022 
str(wider.birds.after[,c(1:5)])
str(scores(nmds.after))
# k=2, Run 20 stress 0.1019276 
figure.nmds.after <- data.frame(wider.birds.after[,c(1:5)], scores(nmds.after))
ggplot(figure.nmds.after) +
  theme_classic() +
  geom_point(aes(x = NMDS1, y = NMDS2, shape = treatment, color = treatment)) +
  scale_color_manual(values = c("dark green", "blue", "red")) +
  
  
#ggsave(filename = paste0("NMDS.after.only.", 
#                           format(Sys.time(), "%d-%b-%Y %H.%M"), ".png"), 
#         dpi = "print", units = "mm", device="png")

ggplot(figure.nmds.after) +
  theme_classic() +
  geom_text(aes(x = NMDS1, y = NMDS2, label = plot, color = treatment)) +
  scale_color_manual(values = c("dark green", "blue", "red")) 

#ggsave(filename = paste0("NMDS.after.only.w.plot.labels.", 
#                         format(Sys.time(), "%d-%b-%Y %H.%M"), ".png"), 
#       dpi = "print", units = "mm", device="png")

# PERMANOVA:

# for stage x treatment interactions:

data.sp <- wider.birds[,-c(1:5)]
data.env <- wider.birds[,c(1:5)]

data.sp.before <- wider.birds.before[,-c(1:5)]
data.env.before <- wider.birds.before[,c(1:5)]

data.sp.after <- wider.birds.after[,-c(1:5)]
data.env.after <- wider.birds.after[,c(1:5)]

#adonis(data.sp~stage*treatment, data=data.env, permutations = 999)

#adonis(data.sp.before~treatment, data=data.env.before, permutations = 999)

#adonis(data.sp.after~treatment, data=data.env.after, permutations = 999)

# note that the NMDS graph looks weird, and it's possible the PERMANOVA is being skwewd by two outliers - G09 and D28
# consider redoing with total species community observed at plot

