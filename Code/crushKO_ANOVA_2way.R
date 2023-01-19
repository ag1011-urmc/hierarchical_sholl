library(rstatix)

aggregate_dat = retina_dat %>% 
    select(Radius, Inters., genotype, group, Animal_ID, Cell_ID) %>%
    group_by(Cell_ID) %>%
    mutate(curve_obs = seq_along(Cell_ID)) %>%
    ungroup() %>%
    select(-Cell_ID) %>% 
    group_by(genotype, Animal_ID, group, curve_obs) %>%
    summarise(Inters. = mean(`Inters.`),
              Radius = mean(Radius)) %>%
    ungroup

max_dat = aggregate_dat %>% 
    select(-curve_obs) %>%
    group_by(Animal_ID, group, genotype) %>%
    summarise(N_max = max(Inters.),
              R_max = Radius[which.max(Inters.)])



m = anova_test(data = as.data.frame(max_dat), formula = N_max ~ group * genotype + Error(Animal_ID / group))
m

m = anova_test(data = as.data.frame(max_dat), formula = R_max ~ group * genotype + Error(Animal_ID / group))
m


gg_dat = retina_dat %>% 
    select(Radius, Inters., genotype, group, Animal_ID, Cell_ID) %>%
    group_by(Cell_ID) %>%
    mutate(curve_obs = seq_along(Cell_ID)) %>%
    select(-Cell_ID) %>% 
    group_by(genotype, group, curve_obs) %>%
    summarise(Inters. = mean(`Inters.`),
              Radius = mean(Radius)) %>%
    ungroup

p = ggplot(data = gg_dat, aes(x = Radius, y = Inters., group = paste(group, genotype, sep = '/'), color = paste(group, genotype, sep = '/'))) + geom_line()
ggsave("tmp.png", p)


