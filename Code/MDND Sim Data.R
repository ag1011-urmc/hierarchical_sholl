library(truncnorm)
library(dplyr)
library(ggplot2)

# Function to simulate data under model 2, modified to include cell-level data
# Arguments:
# theta_pop: numeric vector of length 4. true population level parameters (alpha_1, alpha_2, gamma, tau)
# group_effect: numeric vector of length 4. true group effects on (alpha_1, alpha_2, gamma, tau) at group level
# interaction_effect: numeric vector of length 4. true interaction effects on (alpha_1, alpha_2, gamma, tau) at group level
# sigma_pop: numeric vector of length 4. population level standard deviation parameters corresponding to (alpha_1, alpha_2, gamma, tau)
# sigma_cond_side: numeric vector of length 4. group level standard deviation parameters corresponding to (alpha_1, alpha_2, gamma, tau)
# sigma_animal: numeric vector of length 4. animal level standard deviation parameters corresponding to (alpha_1, alpha_2, gamma, tau)
# sigma_cell: numeric vector of length 4. cell level standard deviation parameters corresponding to (alpha_1, alpha_2, gamma, tau)
# N.Cond_Side: Number of groups
# Animal.Per.Cond_Side: Animals per group
# Cell.Per.Animal: cells per animal
# theta_lb: bounds on parameter space, no need to stray from default
simulate_crushknockout = function(theta_pop, 
                                  group_effect = NULL,
                                  interaction_effect = NULL, 
                                  sigma_pop, 
                                  sigma_cond_side,
                                  sigma_animal,
                                  sigma_cell,
                                  N.Cond_Side, 
                                  Animal.Per.Cond_Side,
                                  Cell.Per.Animal,
                                  theta_lb = c(-Inf, -Inf, 0, 0), 
                                  theta_ub = c(0, 0, 100, Inf)){
    #parameter order:
    #a1
    #a2
    #g
    #t
    
    Cond_ID = c(0, 0, 1, 1)
    Side_ID = c(0, 1, 0, 1)
    
    #group effects
    effect_cond_1 = t(apply(replicate(4, Cond_ID), 1, function(rowvec) rowvec * group_effect[['cond_effect_1']]))
    effect_cond_0 = t(apply(replicate(4, (1 - Cond_ID)), 1, function(rowvec) rowvec * group_effect[['cond_effect_0']]))
    effect_cond = effect_cond_1 + effect_cond_0
    
    effect_side_1 = t(apply(replicate(4, Side_ID), 1, function(rowvec) rowvec * group_effect[['side_effect_1']]))
    effect_side_0 = t(apply(replicate(4, (1 - Side_ID)), 1, function(rowvec) rowvec * group_effect[['side_effect_0']]))
    effect_side = effect_side_1 + effect_side_0
    
    #interaction effect
    effect_int_00 = t(apply(replicate(4, (1 - Cond_ID) * (1 - Side_ID)), 1, function(rowvec) rowvec * interaction_effect[['interaction_effect_00']]))
    effect_int_01 = t(apply(replicate(4, (1 - Cond_ID) * (Side_ID)), 1, function(rowvec) rowvec * interaction_effect[['interaction_effect_01']]))
    effect_int_10 = t(apply(replicate(4, (Cond_ID) * (1 - Side_ID)), 1, function(rowvec) rowvec * interaction_effect[['interaction_effect_10']]))
    effect_int_11 = t(apply(replicate(4, Cond_ID * (Side_ID)), 1, function(rowvec) rowvec * interaction_effect[['interaction_effect_11']]))
    effect_int = effect_int_00 + effect_int_01 + effect_int_10 + effect_int_11
    
    #overall effect
    effect = effect_cond + effect_side + effect_int
    
    theta_cond_side = c()
    theta_animal = c()
    theta_cell = c()
    Cond_Side_ID = c()
    for(i in 1:4){
        theta_cond_side = rbind(theta_cond_side, rtruncnorm(n = 4, 
                                                            mean = theta_pop, 
                                                            sd = sigma_cond_side, 
                                                            a = theta_lb,
                                                            b = theta_ub) + effect[i,])
        theta_animal = rbind(theta_animal, t(replicate(Animal.Per.Cond_Side, 
                                                       rtruncnorm(n = 4,
                                                                  mean = theta_cond_side[i,],
                                                                  sd = sigma_animal, 
                                                                  a = theta_lb,
                                                                  b = theta_ub))))
        Cond_Side_ID = c(Cond_Side_ID, rep(i, Animal.Per.Cond_Side))
    }
    N.Animal = nrow(theta_animal)
    
    
    theta_cell = c()
    Animal_ID = c()
    N.Cell = Cell.Per.Animal * N.Animal
    #loop through N.cell here because effect is defined for each cell
    Animal_ID = rep(seq_len(N.Animal), each = Cell.Per.Animal)
    Cond_Side_ID_forCell = rep(seq_len(N.Cond_Side), each = Cell.Per.Animal * Animal.Per.Cond_Side)
    for(i in seq_len(N.Cell)){
        
        theta_cell = rbind(theta_cell, rtruncnorm(n = 4, 
                                                  mean = theta_animal[Animal_ID[i],], 
                                                  sd = sigma_cell, 
                                                  a = theta_lb, 
                                                  b = theta_ub))
    }
    
    N.Per.Cell = theta_ub[3]
    x = seq_len(N.Per.Cell)
    y = c()
    Cell_ID = c()
    for(i in seq_len(N.Cell)){
        a1 = theta_cell[i, 1]
        a2 = theta_cell[i, 2]
        gg = theta_cell[i, 3]
        tt = theta_cell[i, 4]
        eta = ifelse(x <= gg, a1 * (gg - x)^2 + tt, a2 * (x - gg)^2 + tt)
        y = c(y, rpois(N.Per.Cell, exp(eta)))
        Cell_ID = c(Cell_ID, rep(i, N.Per.Cell))
    }
    x = rep(x, N.Cell)
    N = length(y)
    
    d = list(y = y,
             x = x,
             Cell_ID = Cell_ID,
             Animal_ID = Animal_ID,
             Cond_Side_ID = Cond_Side_ID,
             Side_ID = Side_ID,
             Cond_ID = Cond_ID,
             Cond_Side_ID_forCell = Cond_Side_ID_forCell,
             N = N,
             N.Cell = N.Cell,
             N.Animal = N.Animal,
             N.Cond_Side = N.Cond_Side)
    
    true_pars = list(theta_pop = theta_pop, theta_cond_side = theta_cond_side, theta_animal = theta_animal, theta_cell = theta_cell, sigma_pop = sigma_pop, sigma_cond_side = sigma_cond_side, sigma_animal = sigma_animal, sigma_cell = sigma_cell)
    
    return(list(d = d, true_pars = true_pars))
}

################################################
### Simulate Data
################################################

if(!dir.exists('./Simulated_Data/MDND')) dir.create('./Simulated_Data/MDND')

## define baseline SD parameters
sigma_pop_baseline = c(0.0001, 0.0001, 1, 0.1)
sigma_cond_side_baseline = c(0.000025, 0.000025, 0.5, 0.05)
sigma_animal_baseline = c(0.0001, 0.0001, 1, 0.1)
sigma_cell_baseline = c(0.0001, 0.0001, 1, 0.1)



grandparent.seed = 983724952
set.seed(grandparent.seed)
n.scenarios = 20
parent.seeds = sample(seq_len(1e8), n.scenarios)
n_datasets = 50


### Scenario 1 (Scenario 2 in paper)
## Condition effect only
scen = 1
set.seed(parent.seeds[scen])
if(!dir.exists(paste0("./Simulated_Data/MDND/scenario_", scen))) dir.create(paste0("./Simulated_Data/MDND/scenario_", scen))
if(!dir.exists(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots"))) dir.create(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots"))
seeds = sample(seq_len(1e8), n_datasets)
for(i in seq_len(n_datasets)){
    set.seed(seeds[i])
    sim_out = simulate_crushknockout(theta_pop = c(-0.002, -0.002, 30, 2), 
                                     group_effect = list(cond_effect_0 = c(0, 0, 0, 0),
                                                         cond_effect_1 = c(0, 0, 0, 0.5),
                                                         side_effect_0 = c(0, 0, 0, 0),
                                                         side_effect_1 = c(0, 0, 0, 0)),
                                     interaction_effect = list(interaction_effect_00 = c(0, 0, 0, 0),
                                                               interaction_effect_01 = c(0, 0, 0, 0),
                                                               interaction_effect_10 = c(0, 0, 0, 0),
                                                               interaction_effect_11 = c(0, 0, 0, 0)),
                                     sigma_pop = sigma_pop_baseline, #this doesnt do anything
                                     sigma_cond_side = sigma_cond_side_baseline,
                                     sigma_animal = sigma_animal_baseline,
                                     sigma_cell = sigma_cell_baseline,
                                     N.Cond_Side = 4, 
                                     Animal.Per.Cond_Side = 5,
                                     Cell.Per.Animal = 10,
                                     theta_lb = c(-Inf, -Inf, 0, 0), 
                                     theta_ub = c(0, 0, 100, Inf))
    saveRDS(sim_out, paste0("./Simulated_Data/MDND/scenario_", scen, '/', seeds[i], '.RDS'))
    
    #save every 10 plots for paper or general viewing later
    if(i %% 50 == 0){
        d = sim_out$d
        id_dat = data.frame(Cell_ID = seq_len(d$N.Cell), Cond_Side_ID_forCell = d$Cond_Side_ID_forCell, Animal_ID = d$Animal_ID)
        ggdat = data.frame(x = d$x, y = d$y, Cell_ID = d$Cell_ID) %>% left_join(x = ., y = id_dat, by = "Cell_ID")
        
        p = ggplot(data = ggdat, aes(x = x, y = y, group = Cell_ID, color = Cond_Side_ID_forCell)) + 
            geom_point(size = 0.4, alpha = 0.5) + 
            facet_wrap(~Cond_Side_ID_forCell, nrow = 1)+  
            guides(color = guide_legend(title = "Cond / Side", override.aes = list(size = 1))) + 
            theme(legend.position="bottom")
        
        ggsave(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots/", seeds[i], ".png"), p)
    }
}


### Scenario 2 (Scenario 3 in paper)
## Side effect only
scen = 2
set.seed(parent.seeds[scen])
if(!dir.exists(paste0("./Simulated_Data/MDND/scenario_", scen))) dir.create(paste0("./Simulated_Data/MDND/scenario_", scen))
if(!dir.exists(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots"))) dir.create(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots"))
seeds = sample(seq_len(1e8), n_datasets)
for(i in seq_len(n_datasets)){
    set.seed(seeds[i])
    sim_out = simulate_crushknockout(theta_pop = c(-0.002, -0.002, 30, 2), 
                                     group_effect = list(cond_effect_0 = c(0, 0, 0, 0),
                                                         cond_effect_1 = c(0, 0, 0, 0),
                                                         side_effect_0 = c(0, 0, 0, 0),
                                                         side_effect_1 = c(0, 0, 0, -0.25)),
                                     interaction_effect = list(interaction_effect_00 = c(0, 0, 0, 0),
                                                               interaction_effect_01 = c(0, 0, 0, 0),
                                                               interaction_effect_10 = c(0, 0, 0, 0),
                                                               interaction_effect_11 = c(0, 0, 0, 0)),
                                     sigma_pop = sigma_pop_baseline, #this doesnt do anything
                                     sigma_cond_side = sigma_cond_side_baseline,
                                     sigma_animal = sigma_animal_baseline,
                                     sigma_cell = sigma_cell_baseline,
                                     N.Cond_Side = 4, 
                                     Animal.Per.Cond_Side = 5,
                                     Cell.Per.Animal = 10,
                                     theta_lb = c(-Inf, -Inf, 0, 0), 
                                     theta_ub = c(0, 0, 100, Inf))
    saveRDS(sim_out, paste0("./Simulated_Data/MDND/scenario_", scen, '/', seeds[i], '.RDS'))
    
    #save every 10 plots for paper or general viewing later
    if(i %% 50 == 0){
        d = sim_out$d
        id_dat = data.frame(Cell_ID = seq_len(d$N.Cell), Cond_Side_ID_forCell = d$Cond_Side_ID_forCell, Animal_ID = d$Animal_ID)
        ggdat = data.frame(x = d$x, y = d$y, Cell_ID = d$Cell_ID) %>% left_join(x = ., y = id_dat, by = "Cell_ID")
        
        p = ggplot(data = ggdat, aes(x = x, y = y, group = Cell_ID, color = Cond_Side_ID_forCell)) + 
            geom_point(size = 0.4, alpha = 0.5) + 
            facet_wrap(~Cond_Side_ID_forCell, nrow = 1)+  
            guides(color = guide_legend(title = "Cond / Side", override.aes = list(size = 1))) + 
            theme(legend.position="bottom")
        
        ggsave(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots/", seeds[i], ".png"), p)
    }
}



### Scenario 3 (Not included in paper)
## Side and Condition effect only
scen = 3
set.seed(parent.seeds[scen])
if(!dir.exists(paste0("./Simulated_Data/MDND/scenario_", scen))) dir.create(paste0("./Simulated_Data/MDND/scenario_", scen))
if(!dir.exists(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots"))) dir.create(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots"))
seeds = sample(seq_len(1e8), n_datasets)
for(i in seq_len(n_datasets)){
    set.seed(seeds[i])
    sim_out = simulate_crushknockout(theta_pop = c(-0.002, -0.002, 30, 2), 
                                     group_effect = list(cond_effect_0 = c(0, 0, 0, 0),
                                                         cond_effect_1 = c(0, 0, 0, 0.5),
                                                         side_effect_0 = c(0, 0, 0, 0),
                                                         side_effect_1 = c(0, 0, 0, -0.25)),
                                     interaction_effect = list(interaction_effect_00 = c(0, 0, 0, 0),
                                                               interaction_effect_01 = c(0, 0, 0, 0),
                                                               interaction_effect_10 = c(0, 0, 0, 0),
                                                               interaction_effect_11 = c(0, 0, 0, 0)),
                                     sigma_pop = sigma_pop_baseline, #this doesnt do anything
                                     sigma_cond_side = sigma_cond_side_baseline,
                                     sigma_animal = sigma_animal_baseline,
                                     sigma_cell = sigma_cell_baseline,
                                     N.Cond_Side = 4, 
                                     Animal.Per.Cond_Side = 5,
                                     Cell.Per.Animal = 10,
                                     theta_lb = c(-Inf, -Inf, 0, 0), 
                                     theta_ub = c(0, 0, 100, Inf))
    saveRDS(sim_out, paste0("./Simulated_Data/MDND/scenario_", scen, '/', seeds[i], '.RDS'))
    
    #save every 10 plots for paper or general viewing later
    if(i %% 50 == 0){
        d = sim_out$d
        id_dat = data.frame(Cell_ID = seq_len(d$N.Cell), Cond_Side_ID_forCell = d$Cond_Side_ID_forCell, Animal_ID = d$Animal_ID)
        ggdat = data.frame(x = d$x, y = d$y, Cell_ID = d$Cell_ID) %>% left_join(x = ., y = id_dat, by = "Cell_ID")
        
        p = ggplot(data = ggdat, aes(x = x, y = y, group = Cell_ID, color = Cond_Side_ID_forCell)) + 
            geom_point(size = 0.4, alpha = 0.5) + 
            facet_wrap(~Cond_Side_ID_forCell, nrow = 1)+  
            guides(color = guide_legend(title = "Cond / Side", override.aes = list(size = 1))) + 
            theme(legend.position="bottom")
        
        ggsave(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots/", seeds[i], ".png"), p)
    }
}

### Scenario 4 (scenario 4 in paper)
## Side, Condition, and Interaction effects
scen = 4
set.seed(parent.seeds[scen])
if(!dir.exists(paste0("./Simulated_Data/MDND/scenario_", scen))) dir.create(paste0("./Simulated_Data/MDND/scenario_", scen))
if(!dir.exists(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots"))) dir.create(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots"))
seeds = sample(seq_len(1e8), n_datasets)
for(i in seq_len(n_datasets)){
    set.seed(seeds[i])
    sim_out = simulate_crushknockout(theta_pop = c(-0.002, -0.002, 30, 2), 
                                     group_effect = list(cond_effect_0 = c(0, 0, 0, 0),
                                                         cond_effect_1 = c(0, 0, 0, 0.5),
                                                         side_effect_0 = c(0, 0, 0, 0),
                                                         side_effect_1 = c(0, 0, 0, -0.25)),
                                     interaction_effect = list(interaction_effect_00 = c(0, 0, 0, 0),
                                                               interaction_effect_01 = c(0, 0, 0, 0),
                                                               interaction_effect_10 = c(0, 0, 0, 0),
                                                               interaction_effect_11 = c(0, 0, 0, 0.5)),
                                     sigma_pop = sigma_pop_baseline, #this doesnt do anything
                                     sigma_cond_side = sigma_cond_side_baseline,
                                     sigma_animal = sigma_animal_baseline,
                                     sigma_cell = sigma_cell_baseline,
                                     N.Cond_Side = 4, 
                                     Animal.Per.Cond_Side = 5,
                                     Cell.Per.Animal = 10,
                                     theta_lb = c(-Inf, -Inf, 0, 0), 
                                     theta_ub = c(0, 0, 100, Inf))
    saveRDS(sim_out, paste0("./Simulated_Data/MDND/scenario_", scen, '/', seeds[i], '.RDS'))
    
    #save every 10 plots for paper or general viewing later
    if(i %% 50 == 0){
        d = sim_out$d
        id_dat = data.frame(Cell_ID = seq_len(d$N.Cell), Cond_Side_ID_forCell = d$Cond_Side_ID_forCell, Animal_ID = d$Animal_ID)
        ggdat = data.frame(x = d$x, y = d$y, Cell_ID = d$Cell_ID) %>% left_join(x = ., y = id_dat, by = "Cell_ID")
        
        p = ggplot(data = ggdat, aes(x = x, y = y, group = Cell_ID, color = Cond_Side_ID_forCell)) + 
            geom_point(size = 0.4, alpha = 0.5) + 
            facet_wrap(~Cond_Side_ID_forCell, nrow = 1)+  
            guides(color = guide_legend(title = "Cond / Side", override.aes = list(size = 1))) + 
            theme(legend.position="bottom")
        
        ggsave(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots/", seeds[i], ".png"), p)
    }
}

### Scenario 5 (Scenario 5 in paper)
## Side, Condition, and Interaction effects, double the cells-per-animal
scen = 5
set.seed(parent.seeds[scen])
if(!dir.exists(paste0("./Simulated_Data/MDND/scenario_", scen))) dir.create(paste0("./Simulated_Data/MDND/scenario_", scen))
if(!dir.exists(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots"))) dir.create(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots"))
seeds = sample(seq_len(1e8), n_datasets)
for(i in seq_len(n_datasets)){
    set.seed(seeds[i])
    sim_out = simulate_crushknockout(theta_pop = c(-0.002, -0.002, 30, 2), 
                                     group_effect = list(cond_effect_0 = c(0, 0, 0, 0),
                                                         cond_effect_1 = c(0, 0, 0, 0.5),
                                                         side_effect_0 = c(0, 0, 0, 0),
                                                         side_effect_1 = c(0, 0, 0, -0.25)),
                                     interaction_effect = list(interaction_effect_00 = c(0, 0, 0, 0),
                                                               interaction_effect_01 = c(0, 0, 0, 0),
                                                               interaction_effect_10 = c(0, 0, 0, 0),
                                                               interaction_effect_11 = c(0, 0, 0, 0.5)),
                                     sigma_pop = sigma_pop_baseline, #this doesnt do anything
                                     sigma_cond_side = sigma_cond_side_baseline,
                                     sigma_animal = sigma_animal_baseline,
                                     sigma_cell = sigma_cell_baseline,
                                     N.Cond_Side = 4, 
                                     Animal.Per.Cond_Side = 5,
                                     Cell.Per.Animal = 20,
                                     theta_lb = c(-Inf, -Inf, 0, 0), 
                                     theta_ub = c(0, 0, 100, Inf))
    saveRDS(sim_out, paste0("./Simulated_Data/MDND/scenario_", scen, '/', seeds[i], '.RDS'))
    
    #save every 10 plots for paper or general viewing later
    if(i %% 50 == 0){
        d = sim_out$d
        id_dat = data.frame(Cell_ID = seq_len(d$N.Cell), Cond_Side_ID_forCell = d$Cond_Side_ID_forCell, Animal_ID = d$Animal_ID)
        ggdat = data.frame(x = d$x, y = d$y, Cell_ID = d$Cell_ID) %>% left_join(x = ., y = id_dat, by = "Cell_ID")
        
        p = ggplot(data = ggdat, aes(x = x, y = y, group = Cell_ID, color = Cond_Side_ID_forCell)) + 
            geom_point(size = 0.4, alpha = 0.5) + 
            facet_wrap(~Cond_Side_ID_forCell, nrow = 1)+  
            guides(color = guide_legend(title = "Cond / Side", override.aes = list(size = 1))) + 
            theme(legend.position="bottom")
        
        ggsave(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots/", seeds[i], ".png"), p)
    }
}


### Scenario 6 (scenario 6 in paper)
## Side, Condition, and Interaction effects, larger SD at animal level for tau
scen = 6
set.seed(parent.seeds[scen])
if(!dir.exists(paste0("./Simulated_Data/MDND/scenario_", scen))) dir.create(paste0("./Simulated_Data/MDND/scenario_", scen))
if(!dir.exists(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots"))) dir.create(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots"))
seeds = sample(seq_len(1e8), n_datasets)
for(i in seq_len(n_datasets)){
    set.seed(seeds[i])
    sim_out = simulate_crushknockout(theta_pop = c(-0.002, -0.002, 30, 2), 
                                     group_effect = list(cond_effect_0 = c(0, 0, 0, 0),
                                                         cond_effect_1 = c(0, 0, 0, 0.5),
                                                         side_effect_0 = c(0, 0, 0, 0),
                                                         side_effect_1 = c(0, 0, 0, -0.25)),
                                     interaction_effect = list(interaction_effect_00 = c(0, 0, 0, 0),
                                                               interaction_effect_01 = c(0, 0, 0, 0),
                                                               interaction_effect_10 = c(0, 0, 0, 0),
                                                               interaction_effect_11 = c(0, 0, 0, 0.5)),
                                     sigma_pop = sigma_pop_baseline, #this doesnt do anything
                                     sigma_cond_side = sigma_cond_side_baseline,
                                     sigma_animal = c(0.0001, 0.0001, 1, 0.25),
                                     sigma_cell = sigma_cell_baseline,
                                     N.Cond_Side = 4, 
                                     Animal.Per.Cond_Side = 5,
                                     Cell.Per.Animal = 10,
                                     theta_lb = c(-Inf, -Inf, 0, 0), 
                                     theta_ub = c(0, 0, 100, Inf))
    saveRDS(sim_out, paste0("./Simulated_Data/MDND/scenario_", scen, '/', seeds[i], '.RDS'))
    
    #save every 10 plots for paper or general viewing later
    if(i %% 50 == 0){
        d = sim_out$d
        id_dat = data.frame(Cell_ID = seq_len(d$N.Cell), Cond_Side_ID_forCell = d$Cond_Side_ID_forCell, Animal_ID = d$Animal_ID)
        ggdat = data.frame(x = d$x, y = d$y, Cell_ID = d$Cell_ID) %>% left_join(x = ., y = id_dat, by = "Cell_ID")
        
        p = ggplot(data = ggdat, aes(x = x, y = y, group = Cell_ID, color = Cond_Side_ID_forCell)) + 
            geom_point(size = 0.4, alpha = 0.5) + 
            facet_wrap(~Cond_Side_ID_forCell, nrow = 1)+  
            guides(color = guide_legend(title = "Cond / Side", override.aes = list(size = 1))) + 
            theme(legend.position="bottom")
        
        ggsave(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots/", seeds[i], ".png"), p)
    }
}

### Scenario 7 (scenario 1 in paper)
## No effects
scen = 7
set.seed(parent.seeds[scen])
if(!dir.exists(paste0("./Simulated_Data/MDND/scenario_", scen))) dir.create(paste0("./Simulated_Data/MDND/scenario_", scen))
if(!dir.exists(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots"))) dir.create(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots"))
seeds = sample(seq_len(1e8), n_datasets)
for(i in seq_len(n_datasets)){
    set.seed(seeds[i])
    sim_out = simulate_crushknockout(theta_pop = c(-0.002, -0.002, 30, 2), 
                                     group_effect = list(cond_effect_0 = c(0, 0, 0, 0),
                                                         cond_effect_1 = c(0, 0, 0, 0),
                                                         side_effect_0 = c(0, 0, 0, 0),
                                                         side_effect_1 = c(0, 0, 0, 0)),
                                     interaction_effect = list(interaction_effect_00 = c(0, 0, 0, 0),
                                                               interaction_effect_01 = c(0, 0, 0, 0),
                                                               interaction_effect_10 = c(0, 0, 0, 0),
                                                               interaction_effect_11 = c(0, 0, 0, 0)),
                                     sigma_pop = sigma_pop_baseline, #this doesnt do anything
                                     sigma_cond_side = sigma_cond_side_baseline,
                                     sigma_animal = sigma_animal_baseline,
                                     sigma_cell = sigma_cell_baseline,
                                     N.Cond_Side = 4, 
                                     Animal.Per.Cond_Side = 5,
                                     Cell.Per.Animal = 10,
                                     theta_lb = c(-Inf, -Inf, 0, 0), 
                                     theta_ub = c(0, 0, 100, Inf))
    saveRDS(sim_out, paste0("./Simulated_Data/MDND/scenario_", scen, '/', seeds[i], '.RDS'))
    
    #save ever 10 plots for paper or general viewing later
    if(i %% 50 == 0){
        d = sim_out$d
        id_dat = data.frame(Cell_ID = seq_len(d$N.Cell), Cond_Side_ID_forCell = d$Cond_Side_ID_forCell, Animal_ID = d$Animal_ID)
        ggdat = data.frame(x = d$x, y = d$y, Cell_ID = d$Cell_ID) %>% left_join(x = ., y = id_dat, by = "Cell_ID")
        
        p = ggplot(data = ggdat, aes(x = x, y = y, group = Cell_ID, color = Cond_Side_ID_forCell)) + 
            geom_point(size = 0.4, alpha = 0.5) + 
            facet_wrap(~Cond_Side_ID_forCell, nrow = 1)+  
            guides(color = guide_legend(title = "Cond / Side", override.aes = list(size = 1))) + 
            theme(legend.position="bottom")
        
        ggsave(paste0("./Simulated_Data/MDND/scenario_", scen, "/plots/", seeds[i], ".png"), p)
    }
}




