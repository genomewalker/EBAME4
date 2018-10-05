library(tidyverse)
library(ggpubr)

load(url("http://ebame4.metagenomics.eu/data/ebame4_r-tutorial.Rda"), verbose = TRUE)


# Explore hits taxonomy ---------------------------------------------------

cluster_taxonomy %>%
  group_by(kingdom, category) %>% # we group by the taxonomic rank and the cluster category
  count() %>% # we count how many of the groups
  ungroup() %>%
  mutate(kingdom = fct_reorder(kingdom, n, .fun = "sum"), # we order the taxonomy based on the number of clusters
         category = fct_relevel(category, c("k", "kwp", "gu", "eu"))) %>% # order the categories for visualization
  ggplot(aes(kingdom, n)) +
  geom_col() +
  facet_wrap(~category, nrow = 1) +
  theme_light() +
  rotate() +
  scale_y_continuous(labels = scales::comma) +
  ylab("Number of clusters") +
  xlab("Kingdom")

cluster_taxonomy %>% filter(category == "eu")


# Explore hits temperature range ------------------------------------------

cluster_data %>%
  inner_join(tara_cdata) %>%
  mutate(category = fct_relevel(category, c("k", "kwp", "gu", "eu"))) %>% # order the categories for visualization
  ggplot(aes(temperature)) +
  geom_density(fill = "grey") +
  facet_wrap(~category, nrow = 1) +
  theme_light() +
  ylab("Density") +
  xlab(expression("Temperature " ( degree*C)))
