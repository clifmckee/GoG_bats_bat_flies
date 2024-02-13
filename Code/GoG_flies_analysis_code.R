### Code for Gulf of Guinea bat flies study
### Works: 2024-01-31

# Load packages
library(vegan)
library(ape)
library(phytools)
library(tidyverse)
library(cowplot)
library(ggtree)
library(ggthemes)
library(reshape2)
library(sp)
library(maps)
library(mapdata)
library(binom)
library(MuMIn)
library(readxl)
library(broom)
library(here)

###########################
### Map of study system ###
###########################

# Read in IUCN species distributions
load(file='Data/mammterr.RData')

# Create a list of species names to filter mammterr
allnames = c('Eidolon helvum', 'Rousettus aegyptiacus')
# Filter mammterr down to just binomial names
all.binomial = mammterr$binomial

# Which rows of the data are the species I care about
keep = list()
for(i in 1:length(allnames)){
  keep[[i]] = which(all.binomial == allnames[i])
}
x = keep[[1]]
for(i in 2:length(allnames)){
  x = c(x, keep[[i]])
}
keep = x
myspecies.distr = mammterr[keep,]

# Separate polygons by species, simplify the resolution, and convert to data frames
E_helvum <- fortify(myspecies.distr[myspecies.distr$binomial == 'Eidolon helvum', ]) %>%
  mutate(sp = 'Eidolon helvum')
R_aegyptiacus <- fortify(myspecies.distr[myspecies.distr$binomial == 'Rousettus aegyptiacus', ]) %>%
  mutate(sp = 'Rousettus aegyptiacus')

# Combine species data frames
speciesDF <- bind_rows(E_helvum, R_aegyptiacus)

# Make world map
world <- map_data('worldHires')
map_A <- ggplot() +
    geom_polygon(data=speciesDF, aes(x=long, y=lat, group=group, fill=sp), alpha=0.5) +
    geom_map(data=world, map=world, aes(map_id=region), fill=NA, colour='grey50', linewidth=0.1) +
    geom_rect(aes(xmin=0, xmax=9.5, ymin=-1.8, ymax=6), colour='black', fill=NA) +
    scale_fill_manual(name='Species', values=c('#009E73', '#E69F00')) +
    scale_x_continuous(name='Longitude', limits=c(-20, 60), breaks=seq(-20, 60, 20)) +
    scale_y_continuous(name='Latitude', limits=c(-40, 40), breaks=seq(-40, 40, 20)) +
    theme_cowplot(font_size=12) +
    theme(legend.position='bottom')

# Make Gulf of Guinea map
map_B <- ggplot() +
    geom_polygon(data=speciesDF, aes(x=long, y=lat, group=group, fill=sp), alpha=0.5) +
    geom_map(data=world, map=world, aes(map_id=region), fill=NA, colour='grey50', linewidth=0.1) +
    scale_fill_manual(name='Species', values=c('#009E73', '#E69F00')) +
    scale_x_continuous(name='Longitude', breaks=seq(0, 10, 1)) +
    scale_y_continuous(name='Latitude', breaks=seq(-2, 6, 1)) +
    annotate('text', x=0.8, y=5.4, label='Ghana', size=4) +
    annotate('text', x=7.7, y=3.5, label='Bioko', size=4) +
    annotate('text', x=6.3, y=1.7, label='Príncipe', size=4) +
    annotate('text', x=5.3, y=0.3, label='São Tomé', size=4) +
    annotate('text', x=4.6, y=-1.4, label='Annobón', size=4) +
    annotate('segment', x=8.92, xend=9.11, y=3.8, yend=3.96, linetype=3) +
    annotate('segment', x=7.46, xend=8.53, y=1.71, yend=3.185, linetype=3) +
    annotate('segment', x=7.52, xend=9.27, y=1.61, yend=1.24, linetype=3) +
    annotate('segment', x=6.68, xend=7.35, y=0.44, yend=1.5, linetype=3) +
    annotate('segment', x=6.82, xend=8.62, y=0.27, yend=-0.62, linetype=3) +
    annotate('segment', x=5.655, xend=6.5, y=-1.35, yend=-0.013, linetype=3) +
    annotate('segment', x=5.7, xend=8.96, y=-1.42, yend=-1.36, linetype=3) +
    theme_cowplot(font_size=12) +
    theme(legend.position = "none") +
    coord_cartesian(xlim = c(0, 9.5), ylim = c(-1.8, 6))

# Combine maps
combo_map <- plot_grid(map_A, map_B, labels=c('A', 'B'), label_size=16, ncol=2, axis = "lb", align='hv')
ggsave('Results/GoG_maps.pdf', plot = combo_map, height=5, width=10, units='in')

###############################
### Ectoparasite prevalence ###
###############################

# Read in data and filter to Eidolon helvum
Eh_flies <- read_excel("Data/bat_fly_prevalence.xlsx", sheet = "Islands raw spreadsheet") %>%
  filter(Species == "Eidolon helvum") %>%
  mutate(Island = factor(Island,
                         levels = c("Bioko", "Príncipe", "São Tomé", "Annobón")),
         age_recode = factor(case_when(Age == "A" ~ "adult",
                                Age == "SI" ~ "sexually immature",
                                Age == "J" ~ "juvenile",
                                Age == "Free-flying J/SI" ~ "neonate", 
                                TRUE ~ Age),
                             levels = c("neonate", "juvenile", "sexually immature", "adult")),
         sex_recode = factor(case_when(Gender == "F" ~ "female",
                                       Gender == "M" ~ "male",
                                       TRUE ~ Gender),
                             levels = c("female", "male")))

# Bat fly prevalence across islands for Eidolon helvum
flies_prev_islands <- Eh_flies %>%
  group_by(Island) %>%
  summarize(n_bats = n(),
            bats_w_flies = sum(Nycteribiid_present),
            bats_w_fly_count = sum(!is.na(Nycteribiid_count)),
            avg_flies_per_bat = round(mean(Nycteribiid_count, na.rm = TRUE), 2),
            lower_IQR_flies = quantile(Nycteribiid_count, probs = 0.25, na.rm = TRUE),
            upper_IQR_flies = quantile(Nycteribiid_count, probs = 0.75, na.rm = TRUE)) %>%
  mutate(binom.wilson(x = bats_w_flies, n = n_bats)) %>%
  select(!c(x, n)) %>%
  rename(prevalence = mean, prev_lower = lower, prev_upper = upper) %>%
  mutate(prev_CI = paste0(round(prevalence, 2), " (",
                          round(prev_lower, 2), "–",
                          round(prev_upper, 2), ")"))

# Calculate difference in prevalence between islands for Eidolon helvum
prop.test(x = flies_prev_islands$bats_w_flies,
          n = flies_prev_islands$n_bats,
          alternative = "two.sided") %>%
  tidy()

# Calculate difference in fly counts between islands for Eidolon helvum
kruskal.test(Nycteribiid_count ~ Island, data = Eh_flies) %>%
  tidy()
glm(Nycteribiid_count ~ Island, data = Eh_flies, family = poisson(link = "log")) %>%
  tidy()
ggplot(Eh_flies, aes(x = Island, y = Nycteribiid_count)) +
  geom_boxplot(outlier.color = NA, color = "black", fill = NA, width = 0.5) +
  geom_jitter(shape = 1, size = 0.5, width = 0.25) +
  labs(x = "Location",
       y = "Number of bat flies") +
  theme_cowplot(font_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# Bat fly prevalence across age groups for Eidolon helvum
flies_prev_age <- Eh_flies %>%
  group_by(age_recode) %>%
  summarize(n_bats = n(),
            bats_w_flies = sum(Nycteribiid_present),
            bats_w_fly_count = sum(!is.na(Nycteribiid_count)),
            avg_flies_per_bat = round(mean(Nycteribiid_count, na.rm = TRUE), 2),
            lower_IQR_flies = quantile(Nycteribiid_count, probs = 0.25, na.rm = TRUE),
            upper_IQR_flies = quantile(Nycteribiid_count, probs = 0.75, na.rm = TRUE)) %>%
  mutate(binom.wilson(x = bats_w_flies, n = n_bats)) %>%
  select(!c(x, n)) %>%
  rename(prevalence = mean, prev_lower = lower, prev_upper = upper) %>%
  mutate(prev_CI = paste0(round(prevalence, 2), " (",
                          round(prev_lower, 2), "–",
                          round(prev_upper, 2), ")"))

# Calculate difference in prevalence between ages for Eidolon helvum
prop.test(x = flies_prev_age$bats_w_flies,
          n = flies_prev_age$n_bats,
          alternative = "two.sided") %>%
  tidy()

# Calculate difference in fly counts between islands for Eidolon helvum
kruskal.test(Nycteribiid_count ~ age_recode, data = Eh_flies) %>%
  tidy()
glm(Nycteribiid_count ~ age_recode, data = Eh_flies, family = poisson(link = "log")) %>%
  tidy()
ggplot(Eh_flies, aes(x = age_recode, y = Nycteribiid_count)) +
  geom_boxplot(outlier.color = NA, color = "black", fill = NA, width = 0.5) +
  geom_jitter(shape = 1, size = 0.5, width = 0.25) +
  labs(x = "Age",
       y = "Number of bat flies") +
  theme_cowplot(font_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# Bat fly prevalence across sexes for Eidolon helvum
flies_prev_sex <- Eh_flies %>%
  group_by(sex_recode) %>%
  summarize(n_bats = n(),
            bats_w_flies = sum(Nycteribiid_present),
            bats_w_fly_count = sum(!is.na(Nycteribiid_count)),
            avg_flies_per_bat = round(mean(Nycteribiid_count, na.rm = TRUE), 2),
            lower_IQR_flies = quantile(Nycteribiid_count, probs = 0.25, na.rm = TRUE),
            upper_IQR_flies = quantile(Nycteribiid_count, probs = 0.75, na.rm = TRUE)) %>%
  mutate(binom.wilson(x = bats_w_flies, n = n_bats)) %>%
  select(!c(x, n)) %>%
  rename(prevalence = mean, prev_lower = lower, prev_upper = upper) %>%
  mutate(prev_CI = paste0(round(prevalence, 2), " (",
                          round(prev_lower, 2), "–",
                          round(prev_upper, 2), ")"))

# Calculate difference in prevalence between sexes for Eidolon helvum
prop.test(x = flies_prev_sex$bats_w_flies,
          n = flies_prev_sex$n_bats,
          alternative = "two.sided") %>%
  tidy()

# Calculate difference in fly counts between sexes for Eidolon helvum
kruskal.test(Nycteribiid_count ~ sex_recode, data = Eh_flies) %>%
  tidy()
glm(Nycteribiid_count ~ sex_recode, data = Eh_flies, family = poisson(link = "log")) %>%
  tidy()
ggplot(Eh_flies, aes(x = sex_recode, y = Nycteribiid_count)) +
  geom_boxplot(outlier.color = NA, color = "black", fill = NA, width = 0.5) +
  geom_jitter(shape = 1, size = 0.5, width = 0.25) +
  labs(x = "Age",
       y = "Number of bat flies") +
  theme_cowplot(font_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

#############################
### Bartonella prevalence ###
#############################

# Bat fly Bartonella prevalence
prev <- read.csv('Data/GoG_prevalence.csv', header=T)
prev.CI <- binom.wilson(prev$Pos, prev$Tested, conf.level=0.95)
prev.CI <- cbind(prev[,1:2], prev.CI)
write.csv(prev.CI, 'Results/prev_CI.csv')

# Calculate difference in prevalence between bat species
prop.test(x=c(sum(prev.CI[1:5,]$x), sum(sum(prev.CI[6:9,]$x))),
          n=c(sum(prev.CI[1:5,]$n), sum(sum(prev.CI[6:9,]$n))),
          alternative='two.sided')

# Calculate difference in prevalence between locations for C. greefi
prop.test(x=prev.CI[1:5,]$x, n=prev.CI[1:5,]$n,
          alternative='two.sided')

# Calculate difference in prevalence between locations for E. africana
prop.test(x=prev.CI[6:8,]$x, n=prev.CI[6:8,]$n,
          alternative='two.sided')

# Plot prevalence
ggplot(data=prev.CI, aes(x=Location, y=mean)) +
  geom_pointrange(aes(ymin=lower, ymax=upper), size=1) +
  scale_x_discrete(limits=c('Ghana', 'Bioko', 'Principe', 'Sao Tome', 'Annobon'),
                   labels=c('Ghana', 'Bioko', 'Príncipe', 'São Tomé', 'Annobón'),
                   name='Location') +
  ylab('Prevalence') +
  theme_cowplot(font_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  facet_grid(~Species)
ggsave('Results/GoG_prevalence_by_species.pdf', height=4, width=5, units='in')
ggsave('Results/GoG_prevalence_by_species.png', height=4, width=5, units='in', dpi=300)
ggsave('Results/GoG_prevalence_by_species.tiff', height=4, width=5, units='in', dpi=300)

############################
### Bartonella diversity ###
############################

# Calculation of diversity measures
counts <- read.csv(file='Data/GoG_flies_counts.csv', header=T)
relabund <- counts[,3:11]/apply(counts[,3:11], 1, sum)
richness <- apply(relabund, 1, function(i) sum(i > 0))
shannon <- exp(diversity(relabund, index='shannon'))
simpson <- diversity(relabund, index='invsimpson')
diversity <- data.frame(sp=counts$Species, loc=counts$Location, rich=richness, shan=shannon, simp=simpson)
m.div <- melt(diversity[1:5,], id_vars=c('sp', 'loc'))

# Reshape relative abundance data into long format
colnames(counts) <- c('Species', 'Location', 'E1', 'E2', 'E3', 'E4', 'E5', 'Ew', 'Eh6', 'Eh7', 'Ra1')
m.counts <- melt(counts, id_vars=c('Species', 'Location')) %>%
   filter(variable != "Ra1")

# Stacked bar plot of relative abundance data by Bartonella species
(div_A <- ggplot(data=m.counts, aes(x=Location, y=value, fill=variable)) +
    geom_col(position='fill') +
    scale_fill_viridis_d(direction=-1, option='C', name='Species',
                       labels=c('E1', 'E2', 'E3', 'E4', 'E5', 'Ew', 'Eh6', 'Eh7')) +
    scale_x_discrete(limits=c('Ghana', 'Bioko', 'Principe', 'Sao Tome', 'Annobon'),
                     labels=c('Ghana', 'Bioko', 'Príncipe', 'São Tomé', 'Annobón'),
                     name='Location') +
    ylab('Relative abundance') +
    theme_cowplot(font_size=12) +
    theme(axis.text.x=element_text(angle=45, hjust=1)))

# Alpha diversity
(div_B <- ggplot(data=m.div, aes(x=loc, y=value, fill=variable)) +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(values=c('#1b9e77', '#d95f02', '#7570b3'), name='Index',
                      labels=c('Richness', 'Shannon', 'Simpson')) +
    scale_x_discrete(limits=c('Ghana', 'Bioko', 'Principe', 'Sao Tome', 'Annobon'),
                     labels=c('Ghana', 'Bioko', 'Príncipe', 'São Tomé', 'Annobón'),
                     name='Location') +
    ylab('Alpha diversity') +
    theme_cowplot(font_size=12) +
    theme(axis.text.x=element_text(angle=45, hjust=1)))

# Combine diversity plots
plot_grid(div_A, div_B, labels=c('A', 'B'), label_size=16,
          ncol=2, rel_widths=c(0.49, 0.51), align='v')
ggsave('Results/GoG_div.pdf', height=4, width=10, units='in')
ggsave('Results/GoG_div.png', height=4, width=10, units='in', dpi=300, bg = "white")
ggsave('Results/GoG_div.tiff', height=4, width=10, units='in', dpi=300, bg = "white")

######################################################
### Demographic predictors of Bartonella infection ###
######################################################

# Read in predictor data
predictors <- read_excel('Data/GoG_predictors_redo.xlsx', sheet = 1) %>%
  mutate(Location = factor(case_when(Country_code == "AN" ~ "Annobón",
                                     Country_code == "BI" ~ "Bioko",
                                     Country_code == "GH" ~ "Ghana",
                                     Country_code == "PR" ~ "Príncipe",
                                     Country_code == "ST" ~ "São Tomé",
                                     TRUE ~ Country_code),
                           levels = c("Ghana", "Bioko", "Príncipe", "São Tomé", "Annobón")),
         Age2_recode = factor(case_when(Age2 == "A" ~ "adult",
                                        Age2 == "SI" ~ "sexually immature",
                                        Age2 == "J" ~ "juvenile",
                                        Age2 == "N" ~ "neonate", 
                                        TRUE ~ Age2),
                              levels = c("neonate", "juvenile", "sexually immature", "adult")),
         Sex_recode = factor(case_when(Sex == "F" ~ "female",
                                       Sex == "M" ~ "male",
                                       TRUE ~ Sex),
                             levels = c("female", "male")))
# Remove Bioko
predictors_BI <- predictors %>%
  filter(Location != "Bioko")

# Bartonella prevalence across locations for Eidolon helvum
bart_prev_locations <- predictors %>%
  group_by(Location) %>%
  summarize(n_flies = n(),
            n_pos = sum(Pos)) %>%
  mutate(binom.wilson(x = n_pos, n = n_flies)) %>%
  select(!c(x, n)) %>%
  rename(prevalence = mean, prev_lower = lower, prev_upper = upper) %>%
  mutate(prev_CI = paste0(round(prevalence, 2), " (",
                          round(prev_lower, 2), "–",
                          round(prev_upper, 2), ")"))

# Calculate difference in prevalence between locations for Eidolon helvum
prop.test(x = bart_prev_locations$n_pos,
          n = bart_prev_locations$n_flies,
          alternative = "two.sided") %>%
  tidy()

# Repeat test without Bioko
prop.test(x = predictors_BI %>%
            group_by(Location) %>%
            summarize(n_flies = n(),
                      n_pos = sum(Pos)) %>%
            pull(n_pos),
          n = predictors_BI %>%
            group_by(Location) %>%
            summarize(n_flies = n(),
                      n_pos = sum(Pos)) %>%
            pull(n_flies),
          alternative = "two.sided") %>%
  tidy()

# Plot prevalence by location
location_prev_plot <- ggplot(data = bart_prev_locations, aes(x = Location, y = prevalence)) +
  geom_pointrange(aes(ymin = prev_lower, ymax = prev_upper), size=1) +
  scale_y_continuous(name = "Prevalence", limits = c(0, 1)) +
  theme_cowplot(font_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# Bartonella prevalence across age groups for Eidolon helvum
bart_prev_age <- predictors %>%
  group_by(Age2_recode) %>%
  summarize(n_flies = n(),
            n_pos = sum(Pos)) %>%
  mutate(binom.wilson(x = n_pos, n = n_flies)) %>%
  select(!c(x, n)) %>%
  rename(prevalence = mean, prev_lower = lower, prev_upper = upper) %>%
  mutate(prev_CI = paste0(round(prevalence, 2), " (",
                          round(prev_lower, 2), "–",
                          round(prev_upper, 2), ")"))

# Calculate difference in prevalence between ages for Eidolon helvum
prop.test(x = bart_prev_age$n_pos,
          n = bart_prev_age$n_flies,
          alternative = "two.sided") %>%
  tidy()

# Repeat test without Bioko
prop.test(x = predictors_BI %>%
            group_by(Age2_recode) %>%
            summarize(n_flies = n(),
                      n_pos = sum(Pos)) %>%
            pull(n_pos),
          n = predictors_BI %>%
            group_by(Age2_recode) %>%
            summarize(n_flies = n(),
                      n_pos = sum(Pos)) %>%
            pull(n_flies),
          alternative = "two.sided") %>%
  tidy()

# Plot prevalence by age
age_prev_plot <- ggplot(data = bart_prev_age, aes(x = Age2_recode, y = prevalence)) +
  geom_pointrange(aes(ymin = prev_lower, ymax = prev_upper), size=1) +
  xlab("Age class") +
  scale_y_continuous(name = "Prevalence", limits = c(0, 1)) +
  theme_cowplot(font_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# Bartonella prevalence across sexes for Eidolon helvum
bart_prev_sex <- predictors %>%
  group_by(Sex_recode) %>%
  summarize(n_flies = n(),
            n_pos = sum(Pos)) %>%
  mutate(binom.wilson(x = n_pos, n = n_flies)) %>%
  select(!c(x, n)) %>%
  rename(prevalence = mean, prev_lower = lower, prev_upper = upper) %>%
  mutate(prev_CI = paste0(round(prevalence, 2), " (",
                          round(prev_lower, 2), "–",
                          round(prev_upper, 2), ")"))

# Calculate difference in prevalence between ages for Eidolon helvum
prop.test(x = bart_prev_sex$n_pos,
          n = bart_prev_sex$n_flies,
          alternative = "two.sided") %>%
  tidy()

# Repeat test without Bioko
prop.test(x = predictors_BI %>%
            group_by(Sex_recode) %>%
            summarize(n_flies = n(),
                      n_pos = sum(Pos)) %>%
            pull(n_pos),
          n = predictors_BI %>%
            group_by(Sex_recode) %>%
            summarize(n_flies = n(),
                      n_pos = sum(Pos)) %>%
            pull(n_flies),
          alternative = "two.sided") %>%
  tidy()

# Plot prevalence by age
sex_prev_plot <- ggplot(data = bart_prev_sex, aes(x = Sex_recode, y = prevalence)) +
  geom_pointrange(aes(ymin = prev_lower, ymax = prev_upper), size=1) +
  xlab("Sex") +
  scale_y_continuous(name = "Prevalence", limits = c(0, 1)) +
  theme_cowplot(font_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# Import age distribution data (data from Peel et al., 2017)
ages_Peel <- read.csv(file='Data/GoG_age_distrib.csv', header=T) %>%
  select(!Total) %>%
  pivot_longer(N:A, names_to = "Age", values_to = "n") %>%
  mutate(Location = factor(case_when(Location == "Annobon" ~ "Annobón",
                                     Location == "Principe" ~ "Príncipe",
                                     Location == "Sao Tome" ~ "São Tomé",
                                     TRUE ~ Location),
                           levels = c("Ghana", "Bioko", "Príncipe", "São Tomé", "Annobón")),
         Age2_recode = factor(case_when(Age == "A" ~ "adult",
                                        Age == "SI" ~ "sexually immature",
                                        Age == "J" ~ "juvenile",
                                        Age == "N" ~ "neonate", 
                                        TRUE ~ Age),
                              levels = c("neonate", "juvenile", "sexually immature", "adult")),
         Origin = "Bats censused")

# Age distribution by location (sampled bats)
ages_flies <- predictors %>%
  group_by(Location, Age2_recode, Sex_recode, Bat_ID) %>%
  summarize(n_flies = n(),
            n_pos = sum(Pos)) %>%
  group_by(Location, Age2_recode) %>%
  count() %>%
  mutate(Origin = "Bats sampled with flies")

# Combine age distribution data
ages <- full_join(ages_Peel, ages_flies)

# Stacked bar plot of E. helvum age distribution
age_distrib_plot <- ggplot(data = ages, aes(x = Location, y = n, fill = Age2_recode)) +
  geom_col(position = "fill") +
  scale_fill_viridis_d(direction = -1, option = "D", name = "Age class") +
  ylab("Relative abundance") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(color = "black", fill = "white")) +
  facet_wrap(~Origin)

# Combine prevalence plots
top_line <- plot_grid(location_prev_plot, age_prev_plot, sex_prev_plot,
                      labels = c("A", "B", "C"), label_size = 16,
                      ncol = 3, align = "h")
plot_grid(top_line, age_distrib_plot,
          labels = c("", "D"), label_size = 16,
          nrow = 2)
ggsave('Results/GoG_prev.pdf', height=8, width=10, units='in')
ggsave('Results/GoG_prev.png', height=8, width=10, units='in', dpi=300, bg = "white")
ggsave('Results/GoG_prev.tiff', height=8, width=10, units='in', dpi=300, bg = "white")

#############################
### Isolation by distance ###
#############################

# Input physical distances and bat genetic data
bat_dist <- read.csv('Data/Bat_phys_genet_dist.csv', head=T)

# Calculate Bartonella community dissimilarity
bart.spearman <- as.matrix(1-cor(t(counts[1:5, 3:10]), method='spearman'))
write.csv(bart.spearman, 'Data/spearman_matrix.csv', row.names = FALSE)

bart.bray <- as.matrix(vegdist(counts[1:5, 3:10], method = "bray"))
write.csv(bart.bray, 'Data/bray_matrix.csv', row.names = FALSE)

bart.jaccard <- as.matrix(vegdist(counts[1:5, 3:10], method = "jaccard"))
write.csv(bart.jaccard, 'Data/jaccard_matrix.csv', row.names = FALSE)

# Remove lower triangle (redundant because symmetrical)
bart.spearman[upper.tri(bart.spearman)] <- NA
bart.bray[upper.tri(bart.bray)] <- NA
bart.jaccard[upper.tri(bart.jaccard)] <- NA
# Melt data into long format
spear <- melt(bart.spearman[1:5, 1:5])
bray <- melt(bart.bray[1:5, 1:5])
jaccard <- melt(bart.jaccard[1:5, 1:5])
# Filter out values for self dissimilarity and missing values (from lower triangle)
spear <- spear[which(spear$Var1!=spear$Var2 & spear$value!='NA'),]
bray <- bray[which(bray$Var1!=bray$Var2 & bray$value!='NA'),]
jaccard <- jaccard[which(jaccard$Var1!=jaccard$Var2 & jaccard$value!='NA'),]
# Pull all distances together
distance <- data.frame(label=c('BI-MA', 'PR-MA', 'ST-MA', 'AN-MA', 'PR-BI', 'ST-BI', 'AN-BI', 'ST-PR', 'AN-PR', 'AN-ST'),
                       loc1=spear$Var1, loc2=spear$Var2,
                       spear=spear$value, bray=bray$value, jaccard=jaccard$value,
                       main=bat_dist$Dist, bat_mtDNA=bat_dist$mtDNA, bat_microsats=bat_dist$microsats)
write.csv(distance, 'Results/distance.csv')

# Mainland model fit for Bartonella dissimilarity
spear_fit <- lm(spear~main, data=distance)
sum_spear_fit <- summary(spear_fit); sum_spear_fit
cor.test(distance$spear, distance$main)

bray_fit <- lm(bray~main, data=distance)
sum_bray_fit <- summary(bray_fit); sum_bray_fit
cor.test(distance$bray, distance$main)

jaccard_fit <- lm(jaccard~main, data=distance)
sum_jaccard_fit <- summary(jaccard_fit); sum_jaccard_fit
cor.test(distance$jaccard, distance$main)

# Plot mainland model fit for Bartonella dissimilarity
(bart_dissim <- ggplot(data=distance) +
  geom_smooth(aes(x=main, y=spear), method='lm', colour='grey50', fill='grey50', alpha=0.5) +
  geom_text(aes(x=main, y=spear, label=label), size=3, colour='black', angle = 30) +
  annotate('text', x=300, y=0.25, colour='grey50',
           label=paste('R =', round(sqrt(sum_spear_fit$r.squared), 2),
                       '\nP =', round(sum_spear_fit$coefficients[2, 4], 3))) +
  scale_x_continuous(name='Distance', limits=c(0, 605), breaks=seq(0, 600, 200)) +
  scale_y_continuous(name='Community dissimilarity', limits=c(-0.01, 0.3), breaks=seq(0, 0.3, 0.1)) +
  theme_cowplot(font_size=12))
ggsave('Results/GoG_model_fit.pdf', height=4, width=5, units='in')
ggsave('Results/GoG_model_fit.png', height=4, width=5, units='in', bg='white', dpi=300)
ggsave('Results/GoG_model_fit.tiff', height=4, width=5, units='in', bg='white', dpi=300)

# Mainland model fit for bat mtDNA
mtDNA_fit <- lm(bat_mtDNA~main, data=distance)
sum_mtDNA_fit <- summary(mtDNA_fit); sum_mtDNA_fit
cor.test(distance$bat_mtDNA, distance$main)

# Plot mainland model fit for bat mtDNA
(batIBD_A <- ggplot(data=distance) +
  geom_smooth(aes(x=main, y=bat_mtDNA), method='lm', colour='grey50', fill='grey50', alpha=0.5) +
  geom_text(aes(x=main, y=bat_mtDNA, label=label), size=3, colour='black', angle = 30) +
  annotate('text', x=300, y=0.95, colour='grey50',
           label=paste('R =', round(sqrt(sum_mtDNA_fit$r.squared), 2),
                       '\nP =', round(sum_mtDNA_fit$coefficients[2, 4], 3))) +
  scale_x_continuous(name='Distance', limits=c(0, 605), breaks=seq(0, 600, 200)) +
  scale_y_continuous(name='Genetic distance, bat mtDNA', limits=c(-0.1, 1.1), breaks=seq(0, 1, 0.2)) +
  theme_cowplot(font_size=12))

# Mainland model fit for bat microsatellites
microsats_fit <- lm(bat_microsats~main, data=distance)
sum_microsats_fit <- summary(microsats_fit); sum_microsats_fit
cor.test(distance$bat_microsats, distance$main)

# Plot mainland model fit for bat microsatellites
(batIBD_B <- ggplot(data=distance) +
  geom_smooth(aes(x=main, y=bat_microsats), method='lm', colour='grey50', fill='grey50', alpha=0.5) +
  geom_text(aes(x=main, y=bat_microsats, label=label), size=3, colour='black', angle = 30) +
  annotate('text', x=300, y=0.15, colour='grey50',
           label=paste('R =', round(sqrt(sum_microsats_fit$r.squared), 2),
                       '\nP =', round(sum_microsats_fit$coefficients[2, 4], 3))) +
  scale_x_continuous(name='Distance', limits=c(0, 605), breaks=seq(0, 600, 200)) +
  scale_y_continuous(name='Genetic distance, bat microsatellites', limits=c(-0.03, 0.17), breaks=seq(0, 0.16, 0.04)) +
  theme_cowplot(font_size=12))

# Bartonella fit for bat mtDNA
spear_mtDNA <- lm(spear~bat_mtDNA, data=distance)
sum_spear_mtDNA <- summary(spear_mtDNA); spear_mtDNA
cor.test(distance$spear, distance$bat_mtDNA)

# Plot Bartonella fit for bat mtDNA
(batIBD_C <- ggplot(data=distance) +
  geom_smooth(aes(x=bat_mtDNA, y=spear), method='lm', colour='grey50', fill='grey50', alpha=0.5) +
  geom_text(aes(x=bat_mtDNA, y=spear, label=label), size=3, colour='black', angle = 30) +
  annotate('text', x=0.4, y=0.25, colour='grey50',
           label=paste('R =', round(sqrt(sum_spear_mtDNA$r.squared), 2),
                       '\nP =', round(sum_spear_mtDNA$coefficients[2, 4], 2))) +
  scale_x_continuous(name='Genetic distance, bat mtDNA', limits=c(0, 0.8), breaks=seq(0, 0.8, 0.2)) +
  scale_y_continuous(name='Community dissimilarity', limits=c(0, 0.3), breaks=seq(0, 0.3, 0.06)) +
  theme_cowplot(font_size=12))

# Bartonella fit for bat microsatellites
spear_microsats <- lm(spear~bat_microsats, data=distance)
sum_spear_microsats <- summary(spear_microsats); sum_spear_microsats
cor.test(distance$spear, distance$bat_microsats)

# Plot Bartonella fit for bat microsatellites
(batIBD_D <- ggplot(data=distance) +
  geom_smooth(aes(x=bat_microsats, y=spear), method='lm', colour='grey50', fill='grey50', alpha=0.5) +
  geom_text(aes(x=bat_microsats, y=spear, label=label), size=3, colour='black', angle = 30) +
  annotate('text', x=0.0625, y=0.25, colour='grey50',
           label=paste('R =', round(sqrt(sum_spear_microsats$r.squared), 2),
                       '\nP =', round(sum_spear_microsats$coefficients[2, 4], 2))) +
  scale_x_continuous(name='Genetic distance, bat microsatellites', limits=c(0, 0.13), breaks=seq(0, 0.16, 0.04)) +
  scale_y_continuous(name='Community dissimilarity', limits=c(0, 0.3), breaks=seq(0, 0.3 , 0.06)) +
  theme_cowplot(font_size=12))

# Combine bat IBD plots
plot_grid(batIBD_A, batIBD_B, batIBD_C, batIBD_D, labels=c('A', 'B', 'C', 'D'), label_size=16,
          nrow=2, ncol=2, align='v')
ggsave('Results/GoG_bat_IBD.pdf', height=8, width=10, units='in')
ggsave('Results/GoG_bat_IBD.png', height=8, width=10, bg='white', units='in', dpi=300)
ggsave('Results/GoG_bat_IBD.tiff', height=8, width=10, bg='white', units='in', dpi=300)

# Read in data for Mantel tests
spear.mat <- as.matrix(read.csv('Data/spearman_matrix.csv', head=T))
main.mat <- as.matrix(read.csv('Data/main_matrix.csv', head=F))
mtDNA.mat <- as.matrix(read.csv('Data/mtDNA_matrix.csv', head=F))
micro.mat <- as.matrix(read.csv('Data/microsatellites_matrix.csv', head=F))

# Perform Mantel tests
mantel(spear.mat, main.mat)
mantel(mtDNA.mat, main.mat)
mantel(micro.mat, main.mat)
mantel(spear.mat, mtDNA.mat)
mantel(spear.mat, micro.mat)

################################################
### Phylogenetic analysis of genotyping data ###
################################################

# Color palettes
C_greefi_cols <- colorRampPalette(c('#009E73', '#000000'))(3)
E_africana_cols <- colorRampPalette(c('#E69F00', '#000000'))(6)
D_biannulata_cols <- '#CC79A7'

## Bat fly 16S haplotypes ##

# Read in tree data
e16S_tree_data <- read.csv("Data/tree_metadata_16Se.csv", header = TRUE)

# Read in tree
e16S_tree <- read.tree("Sequencing/16Se model selection rerun/AICc/16Se_seqs_mafft_gb_haplotypes.fasta.contree")

# Root tree at midpoint
e16S_root_tree <- midpoint.root(e16S_tree)

# Combine tip labels and tree data
df_e16S_tree_data <- data.frame(label = sort(e16S_root_tree$tip.label),
                                species = e16S_tree_data$Species,
                                accession = e16S_tree_data$Accession,
                                group = factor(e16S_tree_data$Group))

# Plot tree and save to file
(gen_A <- ggtree(e16S_root_tree) %<+% df_e16S_tree_data +
  geom_tiplab(aes(label = paste(species, accession), color = group),
              key_glyph = "point", size = 3) +
  geom_text2(aes(subset = !isTip & as.numeric(label) >= 0, label = label),
             size = 2, hjust = 1.1, vjust = -0.3, color = "grey50") +
  geom_treescale(x = 0, y = -0.5, width = 0.1, fontsize = 3) +
  xlim(0, 0.3) +
  scale_color_manual(values=c('#000000', '#007E5C', '#B87F00', '#A36085')) +
  theme_tree() +
  theme(legend.position='none',
        plot.subtitle=element_text(size = 10, hjust = 0.125)) +
  labs(subtitle = "ectoparasite mitochondrial 16S rRNA"))

# Read in ectoparasite mitochondrial 16S rRNA data
mt16Se <- read.csv('Data/GoG_genotyping_16Se.csv', head=T)
m.mt16Se <- melt(mt16Se, id.vars=c('Species', 'Location'), measure.vars='Count') %>%
  filter(Location != "Nigeria")

# Plot ectoparasite mitochondrial 16S rRNA data
(gen_B <- ggplot(data=m.mt16Se, aes(x=Location, y=value, fill=Species)) +
    geom_col(position='stack') +
    scale_fill_manual(name='Bat fly\nhaplotype', values=c(C_greefi_cols[1], D_biannulata_cols, E_africana_cols[1:2])) +
    scale_x_discrete(limits=c('Ghana', 'Bioko', 'Principe', 'Sao Tome', 'Annobon'),
                     labels=c('Ghana', 'Bioko', 'Príncipe', 'São Tomé', 'Annobón'),
                     name='Location') +
    scale_y_continuous(name='Count', limits=c(0, 100), breaks=seq(0, 100, 25)) +
    theme_cowplot(font_size=12) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          plot.subtitle=element_text(size = 10)) +
    labs(subtitle = "ectoparasite mitochondrial 16S rRNA"))

## Bat fly cytB haplotypes ##

# Read in tree data
cytB_tree_data <- read.csv("Data/tree_metadata_cytB.csv")

# Read in tree
cytB_tree <- read.tree("Sequencing/cytB model selection/cytB AICc/cytB_seqs_mafft_gb_haplotypes.fasta.contree")

# Root tree at midpoint
cytB_root_tree <- midpoint.root(cytB_tree)

# Combine tip labels and tree data
df_cytB_tree_data <- data.frame(label = sort(cytB_root_tree$tip.label),
                                species = cytB_tree_data$Species,
                                accession = cytB_tree_data$Accession,
                                group = factor(cytB_tree_data$Group))

# Plot tree and save to file
(gen_C <- ggtree(cytB_root_tree) %<+% df_cytB_tree_data +
    geom_tiplab(aes(label = paste(species, accession), color = group),
                key_glyph = "point", size = 3) +
    geom_text2(aes(subset = !isTip & as.numeric(label) >= 0, label = label),
               size = 2, hjust = 1.1, vjust = -0.3, color = "grey50") +
    geom_treescale(x = 0, y = -0.5, width = 0.1, fontsize = 3) +
    xlim(0, 0.6) +
    scale_color_manual(values=c('#000000', '#007E5C', '#B87F00')) +
    theme_tree() +
    theme(legend.position='none',
          plot.subtitle=element_text(size = 10, hjust = 0.125)) +
    labs(subtitle = "ectoparasite mitochondrial cytb"))

# Read in ectoparasite mitochondrial cytB data
mtcytB <- read.csv('Data/GoG_genotyping_cytB.csv', head=T)
m.mtcytB <- melt(mtcytB, id.vars=c('Species', 'Location'), measure.vars='Count') %>%
  filter(Location != "Nigeria")

# Plot ectoparasite mitochondrial cytB data
(gen_D <- ggplot(data=m.mtcytB, aes(x=Location, y=value, fill=Species)) +
    geom_col(position='stack') +
    scale_fill_manual(name='Bat fly\nhaplotype', values=c(C_greefi_cols[1:2], E_africana_cols[1:5])) +
    scale_x_discrete(limits=c('Ghana', 'Bioko', 'Principe', 'Sao Tome', 'Annobon'),
                     labels=c('Ghana', 'Bioko', 'Príncipe', 'São Tomé', 'Annobón'),
                     name='Location') +
    scale_y_continuous(name='Count', limits=c(0, 75), breaks=seq(0, 75, 25)) +
    theme_cowplot(font_size=12) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          plot.subtitle=element_text(size = 10)) +
    labs(subtitle = "ectoparasite mitochondrial cytb"))

## Bat fly symbiont 16S haplotypes ##

# Read in tree data
b16S_tree_data <- read.csv("Data/tree_metadata_16Sb.csv")

# Read in tree
b16S_tree <- read.tree("Sequencing/16Sb model selection/16Sb AICc/16Sb_seqs_mafft_gb_haplotypes.fasta.contree")

# Root tree at midpoint
b16S_root_tree <- midpoint.root(b16S_tree)

# Combine tip labels and tree data
df_b16S_tree_data <- data.frame(label = sort(b16S_root_tree$tip.label),
                                species = b16S_tree_data$Species,
                                accession = b16S_tree_data$Accession,
                                group = factor(b16S_tree_data$Group))

# Plot tree and save to file
(gen_E <- ggtree(b16S_root_tree) %<+% df_b16S_tree_data +
    geom_tiplab(aes(label = paste(species, accession), color = group),
                key_glyph = "point", size = 3) +
    geom_text2(aes(subset = !isTip & as.numeric(label) >= 0, label = label),
               size = 2, hjust = 1.1, vjust = -0.3, color = "grey50") +
    geom_treescale(x = 0, y = -0.5, width = 0.1, fontsize = 3) +
    xlim(0, 0.4) +
    scale_color_manual(values=c('#000000', '#007E5C', '#B87F00')) +
    theme_tree() +
    theme(legend.position='none',
          plot.subtitle=element_text(size = 10, hjust = 0.125)) +
    labs(subtitle = "bacterial symbiont 16S rRNA"))

# Read in ectoparasite bacterial symbiont 16S rRNA data
mt16Sb <- read.csv('Data/GoG_genotyping_16Sb.csv', head=T)
m.mt16Sb <- melt(mt16Sb, id.vars=c('Species', 'Location'), measure.vars='Count') %>%
  filter(Location != "Nigeria")

# Plot ectoparasite bacterial symbiont 16S rRNA data
(gen_F <- ggplot(data=m.mt16Sb, aes(x=Location, y=value, fill=Species)) +
    geom_col(position='stack') +
    scale_fill_manual(name='Symbiont\nhaplotype', values=c(C_greefi_cols[1], E_africana_cols[1:2])) +
    scale_x_discrete(limits=c('Ghana', 'Bioko', 'Principe', 'Sao Tome', 'Annobon'),
                     labels=c('Ghana', 'Bioko', 'Príncipe', 'São Tomé', 'Annobón'),
                     name='Location') +
    scale_y_continuous(name='Count', limits=c(0, 30), breaks=seq(0, 30, 10)) +
    ylab('Relative abundance') +
    theme_cowplot(font_size=12) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          plot.subtitle=element_text(size = 10)) +
    labs(subtitle = "bacterial symbiont 16S rRNA"))

# Combine diversity correlate plots
trees <- plot_grid(gen_A, gen_C, gen_E, labels=c('A', 'C', 'E'), label_size=16,
                   ncol=1, nrow=3)
counts <- plot_grid(gen_B, gen_D, gen_F, labels=c('B', 'D', 'F'), label_size=16,
                    ncol=1, nrow=3, align='v')
plot_grid(trees, counts, ncol=2)
ggsave('Results/GoG_haplotyping.pdf', height=12, width=10, units='in')
ggsave('Results/GoG_haplotyping.png', height=12, width=10, units='in', dpi=300, bg = "white")
ggsave('Results/GoG_haplotyping.tiff', height=12, width=10, units='in', dpi=300, bg = "white")

# Bat fly bacterial symbiont prevalence
sym <- read.csv('Data/GoG_symbionts.csv', header=T)
sym.CI <- binom.wilson(sym$Pos, sym$Tested, conf.level=0.95)
sym.CI <- cbind(sym[,1:2], sym.CI)
write.csv(sym.CI, 'Results/sym_CI.csv')

###-----------------------
### Graphical abstract ---
###-----------------------

abstract_row1 <- plot_grid(map_B, gen_C,
                           ncol = 2, axis = "lb", align = "hv")
abstract_row2 <- plot_grid(div_A +
                             scale_fill_viridis_d(direction=-1, option='C', name='Bartonella\nspecies',
                                                  labels=c('E1', 'E2', 'E3', 'E4', 'E5', 'Ew', 'Eh6', 'Eh7')) +
                             labs(y = "Bartonella relative abundance"),
                           bart_dissim +
                             scale_x_continuous(name='Distance between locations', limits=c(0, 605), breaks=seq(0, 600, 200)) +
                             scale_y_continuous(name='Bartonella community dissimilarity', limits=c(-0.01, 0.3), breaks=seq(0, 0.3, 0.1)),
                           ncol = 2, axis = "lb", align = "hv")

abstract <- plot_grid(abstract_row1, abstract_row2, nrow = 2, rel_heights = c(5/9, 4/9), align = "v")
ggsave('Results/GoG_graphical_abstract.pdf', height = 9, width = 10, units = "in")
ggsave('Results/GoG_graphical_abstract.png', height = 9, width = 10, units = "in", dpi = 300, bg = "white")
ggsave('Results/GoG_graphical_abstract.tiff', height = 9, width = 10, units = "in", dpi = 300, bg = "white")

###-------------------------------------------------------------------
### Supplementary Figure: Bartonella concatenated ftsZ + gltA tree ---
###-------------------------------------------------------------------

# Read in tree data
Bart_conc_tree_data <- read.csv("Data/Bart_conc_tree_metadata.csv")

# Read in tree
Bart_conc_tree <- read.tree("Sequencing/Bartonella conc model selection/BIC/Bartonella_conc_sequences.fasta.contree")

# Root tree at midpoint
root_Bart_conc_tree <- midpoint.root(Bart_conc_tree)

# Order tree data by tip label
o_Bart_conc_tree_data <- Bart_conc_tree_data[match(root_Bart_conc_tree$tip.label, Bart_conc_tree_data$Nickname),]

# Combine tip labels and tree data
df_Bart_conc_tree_data <- data.frame(label = root_Bart_conc_tree$tip.label,
                                     species = o_Bart_conc_tree_data$Bartonella.species,
                                     isolate = o_Bart_conc_tree_data$Isolate.Strain.Clone,
                                     group = factor(o_Bart_conc_tree_data$Group))

# Plot tree
ggtree(root_Bart_conc_tree) %<+% df_Bart_conc_tree_data +
  geom_point(aes(size = as.numeric(label))) +
  geom_tiplab(aes(label = paste(species, isolate), color = group), size = 2.5) +
  geom_text2(aes(subset = !isTip & as.numeric(label) >= 0, label = label),
             size = 2, hjust = 1.1, vjust = -0.3, color = "grey50") +
  geom_treescale(x = 0, y = -2, width = 0.1, fontsize = 3) +
  xlim(0, 0.9) +
  scale_size_binned(range = c(0, 1), breaks = c(50, 80)) +
  scale_color_manual(values = c('black', 'gray50', colorblind_pal()(8)[c(6, 3, 7)])) +
  theme_tree() +
  theme(legend.position = "none")
ggsave("Results/Bart_conc_tree.pdf", height = 12, width = 10, units = "in")
ggsave("Results/Bart_conc_tree.png", height = 12, width = 10, units = "in", dpi = 300)
ggsave("Results/Bart_conc_tree.tiff", height = 12, width = 10, units = "in", dpi = 300)
