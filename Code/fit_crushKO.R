#Script to fit model 3 to knockout data and reproduce corresponding figures

################################################
### load packages
################################################

library(tidyverse)
library(rjags)
library(truncnorm)
library(invgamma)
library(R2jags)
library(parallel)
library(latex2exp)
library(HDInterval)

#######################################
### Define Model
#######################################


m = "model{
        for(i in 1:N){
            y[i] ~ dpois(lambda[i])
            log(lambda[i]) <- eta[i]
            eta[i] = ifelse(x[i] <= g[Cell_ID[i]], a1[Cell_ID[i]] * (g[Cell_ID[i]] - x[i])^2 + t[Cell_ID[i]], a2[Cell_ID[i]] * (x[i] - g[Cell_ID[i]])^2 + t[Cell_ID[i]])
        }
        
        #Cell Level
        for(j in 1:N.Cell){
            ##CODING:
            #KO coded as 1, WT coded as 0
            #Conrol coded as 0, Crush coded as 1
            
            a1[j] ~ dnorm(a1.animal[Animal_ID[j]] + (Group_ID[j] * a1.crush_effect) + (Group_ID[j] * Genotype_ID_4Int[j] * a1.int_effect), a1.tau.cell) T(, 0)
            a2[j] ~ dnorm(a2.animal[Animal_ID[j]] + (Group_ID[j] * a2.crush_effect) + (Group_ID[j] * Genotype_ID_4Int[j] * a2.int_effect), a2.tau.cell) T(, 0)
            g[j] ~ dnorm(g.animal[Animal_ID[j]] + (Group_ID[j] * g.crush_effect) + (Group_ID[j] * Genotype_ID_4Int[j] * g.int_effect), g.tau.cell) T(0, 100)
            t[j] ~ dnorm(t.animal[Animal_ID[j]] + (Group_ID[j] * t.crush_effect) + (Group_ID[j] * Genotype_ID_4Int[j] * t.int_effect), t.tau.cell) T(0,)
        }
        
        #Animal Level
        for(k in 1:N.Animal){
            a1.animal[k] ~ dnorm(a1.genotype[Genotype_ID[k]], a1.tau.animal) T(, 0)
            a2.animal[k] ~ dnorm(a2.genotype[Genotype_ID[k]], a2.tau.animal) T(, 0)
            g.animal[k] ~ dnorm(g.genotype[Genotype_ID[k]], g.tau.animal) T(0, 100)
            t.animal[k] ~ dnorm(t.genotype[Genotype_ID[k]], t.tau.animal) T(0,)
        }
        
        #Genotype Level
        #Hard coded to speed things up
            #assuming 1 = WT, 2 = KO
        a1.genotype[1] ~ dnorm(a1.pop, a1.tau.genotype) T(, 0)
        a2.genotype[1] ~ dnorm(a2.pop, a2.tau.genotype) T(, 0)
        g.genotype[1] ~ dnorm(g.pop, g.tau.genotype) T(0, 100)
        t.genotype[1] ~ dnorm(t.pop, t.tau.genotype) T(0,)
        
        a1.genotype[2] ~ dnorm(a1.pop + a1.ko_effect, a1.tau.genotype) T(, 0)
        a2.genotype[2] ~ dnorm(a2.pop + a2.ko_effect, a2.tau.genotype) T(, 0)
        g.genotype[2] ~ dnorm(g.pop + g.ko_effect, g.tau.genotype) T(0, 100)
        t.genotype[2] ~ dnorm(t.pop + t.ko_effect, t.tau.genotype) T(0,)
        
        
        ### Effects
            ## a1
            a1.ko_effect ~ dnorm(a1.ko_effect_priorMean, a1.ko_effect_priorTau) T(, -a1.pop)
            a1.ko_effect_priorMean ~ dnorm(0, 70000000) T(, -a1.pop)
            a1.ko_effect_priorTau <- pow(a1.ko_effect_priorSD, -2)
            a1.ko_effect_priorSD ~ dt(0, 100000000, 4) T(0,)
            
            a1.crush_effect ~ dnorm(a1.crush_effect_priorMean, a1.crush_effect_priorTau) T(, -max(a1.animal))
            a1.crush_effect_priorMean ~ dnorm(0, 70000000) T(, -max(a1.animal))
            a1.crush_effect_priorTau <- pow(a1.crush_effect_priorSD, -2)
            a1.crush_effect_priorSD ~ dt(0, 100000000, 4) T(0,)
            
            a1.int_effect ~ dnorm(a1.int_effect_priorMean, a1.int_effect_priorTau) T(, -(max(a1.animal) + a1.crush_effect))
            a1.int_effect_priorMean ~ dnorm(0, 70000000) T(, -(max(a1.animal) + a1.crush_effect))
            a1.int_effect_priorTau <- pow(a1.int_effect_priorSD, -2)
            a1.int_effect_priorSD ~ dt(0, 100000000, 4) T(0,)
            
            
            ## a2
            a2.ko_effect ~ dnorm(a2.ko_effect_priorMean, a2.ko_effect_priorTau) T(, -a2.pop)
            a2.ko_effect_priorMean ~ dnorm(0, 70000000) T(, -a2.pop)
            a2.ko_effect_priorTau <- pow(a2.ko_effect_priorSD, -2)
            a2.ko_effect_priorSD ~ dt(0, 100000000, 4) T(0,)
            
            a2.crush_effect ~ dnorm(a2.crush_effect_priorMean, a2.crush_effect_priorTau) T(, -max(a2.animal))
            a2.crush_effect_priorMean ~ dnorm(0, 70000000) T(, -max(a2.animal))
            a2.crush_effect_priorTau <- pow(a2.crush_effect_priorSD, -2)
            a2.crush_effect_priorSD ~ dt(0, 100000000, 4) T(0,)
            
            a2.int_effect ~ dnorm(a2.int_effect_priorMean, a2.int_effect_priorTau) T(, -(max(a2.animal) + a2.crush_effect))
            a2.int_effect_priorMean ~ dnorm(0, 70000000) T(, -(max(a2.animal) + a2.crush_effect))
            a2.int_effect_priorTau <- pow(a2.int_effect_priorSD, -2)
            a2.int_effect_priorSD ~ dt(0, 100000000, 4) T(0,)
            
            
            ## g
            g.ko_effect ~ dnorm(g.ko_effect_priorMean, g.ko_effect_priorTau) T(-g.pop, )
            g.ko_effect_priorMean ~ dnorm(0, 0.007) T(-g.pop, )
            g.ko_effect_priorTau <- pow(g.ko_effect_priorSD, -2)
            g.ko_effect_priorSD ~ dt(0, 0.07, 4) T(0,)
            
            g.crush_effect ~ dnorm(g.crush_effect_priorMean, g.crush_effect_priorTau) T(-min(g.animal), )
            g.crush_effect_priorMean ~ dnorm(0, 0.007) T(-min(g.animal), )
            g.crush_effect_priorTau <- pow(g.crush_effect_priorSD, -2)
            g.crush_effect_priorSD ~ dt(0, 0.07, 4) T(0,)
            
            g.int_effect ~ dnorm(g.int_effect_priorMean, g.int_effect_priorTau) T(-(min(g.animal) + g.crush_effect), )
            g.int_effect_priorMean ~ dnorm(0, 0.007) T(-(min(g.animal) + g.crush_effect), )
            g.int_effect_priorTau <- pow(g.int_effect_priorSD, -2)
            g.int_effect_priorSD ~ dt(0, 0.07, 4) T(0,)
            
            
            ## t
            t.ko_effect ~ dnorm(t.ko_effect_priorMean, t.ko_effect_priorTau) T(-t.pop, )
            t.ko_effect_priorMean ~ dnorm(0, 2) T(-t.pop, )
            t.ko_effect_priorTau <- pow(t.ko_effect_priorSD, -2)
            t.ko_effect_priorSD ~ dt(0, 50, 4) T(0,)
            
            t.crush_effect ~ dnorm(t.crush_effect_priorMean, t.crush_effect_priorTau) T(-min(t.animal), )
            t.crush_effect_priorMean ~ dnorm(0, 2) T(-min(t.animal), )
            t.crush_effect_priorTau <- pow(t.crush_effect_priorSD, -2)
            t.crush_effect_priorSD ~ dt(0, 50, 4) T(0,)
            
            t.int_effect ~ dnorm(t.int_effect_priorMean, t.int_effect_priorTau) T(-(min(t.animal) + t.crush_effect), )
            t.int_effect_priorMean ~ dnorm(0, 2) T(-(min(t.animal) + t.crush_effect), )
            t.int_effect_priorTau <- pow(t.int_effect_priorSD, -2)
            t.int_effect_priorSD ~ dt(0, 50, 4) T(0,)
        
        
        #cell level stuff
        a1.tau.cell <- pow(a1.sd.cell, -2)
        a2.tau.cell <- pow(a2.sd.cell, -2)
        g.tau.cell <- pow(g.sd.cell, -2)
        t.tau.cell <- pow(t.sd.cell, -2)
        
        #half t(1. 4) priors
        a1.sd.cell ~ dt(0, 70000000, 1) T(0, ) # mean, precision, df
        a2.sd.cell ~ dt(0, 70000000, 1) T(0, )
        g.sd.cell ~ dt(0, 0.7, 1) T(0, )
        t.sd.cell ~ dt(0, 70, 1) T(0, )
        
        #animal level stuff
        a1.tau.animal <- pow(a1.sd.animal, -2)
        a2.tau.animal <- pow(a2.sd.animal, -2)
        g.tau.animal <- pow(g.sd.animal, -2)
        t.tau.animal <- pow(t.sd.animal, -2)
        
        a1.sd.animal ~ dt(0, 70000000, 4) T(0, ) # mean, precision, df
        a2.sd.animal ~ dt(0, 70000000, 4) T(0, )
        g.sd.animal ~ dt(0, 0.7, 4) T(0, )
        t.sd.animal ~ dt(0, 70, 4) T(0, )
        
        #genotype level stuff
        a1.tau.genotype <- pow(a1.sd.genotype, -2)
        a2.tau.genotype <- pow(a2.sd.genotype, -2)
        g.tau.genotype <- pow(g.sd.genotype, -2)
        t.tau.genotype <- pow(t.sd.genotype, -2)
        
        a1.sd.genotype ~ dt(0, 1500000000, 4) T(0, ) # mean, precision, df
        a2.sd.genotype ~ dt(0, 1500000000, 4) T(0, )
        g.sd.genotype ~ dt(0, 3, 4) T(0, )
        t.sd.genotype ~ dt(0, 1500, 4) T(0, )

        #pop level stuff
        a1.pop ~ dnorm(0, 100000) T(, 0)
        a2.pop ~ dnorm(0, 100000) T(, 0)
        g.pop ~ dnorm(0, 0.01) T(0, 100)
        t.pop ~ dnorm(0, 0.25) T(0,)
}"


################################################
### Helper functions
################################################

#function to create initial values for jags sampler
inits_func = function(chain, N.Cell, N.Genotype, N.Animal, N.Groups){
    gen_list = function(chain = chain){
        out = list(
            a1 = truncnorm::rtruncnorm(N.Cell, b = 0, mean = 0.002, sd = 0.000001),
            a2 = truncnorm::rtruncnorm(N.Cell, b = 0, mean = 0.002, sd = 0.000001),
            t = truncnorm::rtruncnorm(N.Cell, a = 0, mean = 3, sd = 0.1),
            g = truncnorm::rtruncnorm(N.Cell, a = 0, b = 100, mean = 15, sd = 1)
        )
        
        out[['a1.ko_effect']] = 0
        out[['a2.ko_effect']] = 0
        out[['t.ko_effect']] = 0
        out[['g.ko_effect']] = 0
        
        out[['a1.crush_effect']] = 0
        out[['a2.crush_effect']] = 0
        out[['t.crush_effect']] = 0
        out[['g.crush_effect']] = 0
        
        out[['a1.int_effect']] = 0
        out[['a2.int_effect']] = 0
        out[['t.int_effect']] = 0
        out[['g.int_effect']] = 0
        
        
        out[['a1.ko_effect_priorMean']] = 0
        out[['a2.ko_effect_priorMean']] = 0
        out[['t.ko_effect_priorMean']] = 0
        out[['g.ko_effect_priorMean']] = 0
        
        out[['a1.crush_effect_priorMean']] = 0
        out[['a2.crush_effect_priorMean']] = 0
        out[['t.crush_effect_priorMean']] = 0
        out[['g.crush_effect_priorMean']] = 0
        
        out[['a1.int_effect_priorMean']] = 0
        out[['a2.int_effect_priorMean']] = 0
        out[['t.int_effect_priorMean']] = 0
        out[['g.int_effect_priorMean']] = 0
        
        
        out[['a1.ko_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['a2.ko_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['t.ko_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.1, nu = 4)
        out[['g.ko_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 1, nu = 4)
        
        out[['a1.crush_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['a2.crush_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['t.crush_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.1, nu = 4)
        out[['g.crush_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 1, nu = 4)
        
        out[['a1.int_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['a2.int_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['t.int_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.1, nu = 4)
        out[['g.int_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 1, nu = 4)
        
        
        out[['a1.pop']] = mean(out$a1)
        out[['a2.pop']] = mean(out$a2)
        out[['g.pop']] = mean(out$g)
        out[['t.pop']] = mean(out$t)
        
        out[['a1.sd.genotype']] = LaplacesDemon::rhalft(1, scale = 0.000025, nu = 4)
        out[['a2.sd.genotype']] = LaplacesDemon::rhalft(1, scale = 0.000025, nu = 4)
        out[['g.sd.genotype']] = LaplacesDemon::rhalft(1, scale = 0.5, nu = 4)
        out[['t.sd.genotype']] = LaplacesDemon::rhalft(1, scale = 0.025, nu = 4)
        
        out[['a1.genotype']] = rep(mean(out$a1), N.Genotype)
        out[['a2.genotype']] = rep(mean(out$a2), N.Genotype)
        out[['g.genotype']] = rep(mean(out$g), N.Genotype)
        out[['t.genotype']] = rep(mean(out$t), N.Genotype)
        
        out[['a1.sd.cell']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['a2.sd.cell']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['g.sd.cell']] = LaplacesDemon::rhalft(1, scale = 1, nu = 4)
        out[['t.sd.cell']] = LaplacesDemon::rhalft(1, scale = 0.1, nu = 4)
        
        out[['a1.sd.animal']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['a2.sd.animal']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['g.sd.animal']] = LaplacesDemon::rhalft(1, scale = 1, nu = 4)
        out[['t.sd.animal']] = LaplacesDemon::rhalft(1, scale = 0.1, nu = 4)
        
        out[['a1.animal']] = rep(mean(out$a1), N.Animal)
        out[['a2.animal']] = rep(mean(out$a2), N.Animal)
        out[['g.animal']] = rep(mean(out$g), N.Animal)
        out[['t.animal']] = rep(mean(out$t), N.Animal)
        
        return(out)
    }
    return(switch(chain,
                  "1" = gen_list(chain),
                  "2" = gen_list(chain),
                  "3" = gen_list(chain),
                  "4" = gen_list(chain),
                  "5" = gen_list(chain),
                  "6" = gen_list(chain),
                  "7" = gen_list(chain),
                  "8" = gen_list(chain)
    )
    )
}

#wrapper around jags sampler for parallelization
rjags_samp_wrapper = function(j){
    set.seed(seeds[j])
    temp_model = rjags::jags.model(file = textConnection(m),
                                   data = d,
                                   inits = function(chain) inits_func(chain = chain, N.Cell = d$N.Cell, N.Genotype = d$N.Genotype, N.Animal = d$N.Animal, N.Groups = d$N.Group),
                                   n.adapt = n.adapt,
                                   n.chains = 1)
    effect_name_vec = c("a1.crush.ko.vs.crush.wt", "a2.crush.ko.vs.crush.wt", "g.crush.ko.vs.crush.wt", "t.crush.ko.vs.crush.wt", "a1.control.ko.vs.control.wt", "a2.control.ko.vs.control.wt", "g.control.ko.vs.control.wt", "t.control.ko.vs.control.wt", "a1.crush.vs.control", "a2.crush.vs.control", "g.crush.vs.control", "t.crush.vs.control", "a1.ko.vs.wt", "a2.ko.vs.wt", "g.ko.vs.wt", "t.ko.vs.wt")
    out = rjags::jags.samples(model = temp_model, 
                              variable.names = c(names(inits_func(1, N.Cell = d$N.Cell, N.Genotype = d$N.Genotype, N.Animal = d$N.Animal, N.Groups = d$N.Group)), effect_name_vec), 
                              n.iter = n.samp, 
                              thin = n.thin)
    # saveRDS(out, file = paste0("./MCMC Samples/crushcontrol_hierarchical_reParameterized_chain", j, ".RDS"))
    return(out)
}

########################################################################
### Loading data
########################################################################


dat_dir <- './Data/CrushKO'

#get file names for cell level data
sholl_files <- list.files(dat_dir)[grepl('.csv', list.files(dat_dir))]
sholl_file_names <- strsplit(sholl_files, '\\.') %>% sapply(., FUN = function(list_elem) list_elem[[1]])

retina_dat <- c()
for(j in seq_along(sholl_files)){
    tmp_dat <- read.csv(paste(dat_dir, sholl_files[j], sep = '/')) %>% 
        mutate(Cell_ID = sholl_file_names[j])
    
    tmp_dat2 <- data.frame(Radius = seq(tmp_dat$Radius[length(tmp_dat$Radius)], 100.49, by = 2), 
                           Inters. = 0, 
                           Radius..norm.Area = 0, 
                           Inters..Area = 0, 
                           log.Radius = 0, 
                           log.Inters..Area = 0, 
                           Cell_ID = tmp_dat$Cell_ID[1])
    
    retina_dat <- rbind(retina_dat, tmp_dat)
}

########################################################################
### Reformatting and recoding data
########################################################################

#animal 60268 accidentally got coded as 59976, so replaceing 60268 with 59976 in the lookuptable
animal_ids = as.character(c(59873, 59875, 59876, 59877, 59879, 59968, 60043, 60044, 59976))
spreadsheet_brain_ids = as.character(c(7, 6, 5, 3, 2, 8, "k", 4, 1))
genotype = c("KO", "KO", "WT", "WT", "KO", "WT", "WT", "WT", "KO")
crush_eye = c("left", "left", "left", "right", "right", "left", "left", "left", "left")
lookup_table = data.frame(Animal_ID_braindat = spreadsheet_brain_ids, Animal_ID = animal_ids, genotype, crush_eye)

#need to swtich eye labels (crush/control) on the following animals
#59877
#59879
#60043
#59968 (also exclude this one because the images looked terrible)
for(id in c(59877, 59879, 60043, 59968)){
    old_eye = lookup_table[lookup_table[,"Animal_ID"] == id, 'crush_eye']
    lookup_table[lookup_table[,"Animal_ID"] == id, 'crush_eye'] = ifelse(old_eye == "right", "left", "right")
}






retina_dat = retina_dat %>% 
    mutate(Animal_ID = substr(unlist(lapply(strsplit(Cell_ID, split = 'MAX_C1-'), function(x) x[2])), 1, 5)) %>% #extract animal id from file names
    mutate(which_eye = ifelse(substr(Cell_ID, 14, 14) == "r", "right", "left")) %>% #extract eye from file names
    left_join(., lookup_table, by = "Animal_ID") %>% #join with lookuptable
    mutate(group = ifelse(crush_eye == which_eye, 'Crush', 'Control')) %>% #label as crush or control
    arrange(genotype, Animal_ID, group, Cell_ID, Radius) #sort

#omit animal k
retina_dat = retina_dat %>% filter(Animal_ID != 59968)

### uncomment to perform analysis excluding subjects 60044 and 59879
# retina_dat = retina_dat %>% filter(!(Animal_ID %in% c("60044", "59879")))

cell_id_conversion_table = data.frame(Cell_ID = unique(retina_dat$Cell_ID), Cell_ID_num = seq_along(unique(retina_dat$Cell_ID)))

retina_dat = retina_dat %>% left_join(., y = cell_id_conversion_table, by = "Cell_ID")

retina_dat = retina_dat %>%
    mutate(Animal_ID_factor = as.numeric(as.factor(Animal_ID))) %>%
    mutate(Cell_ID_factor = as.numeric(as.factor(Cell_ID_num))) %>%
    mutate(Genotype_ID_factor = as.numeric(as.factor(genotype))) %>%
    mutate(Group_ID_factor = as.numeric(as.factor(group))) %>% 
    arrange(Genotype_ID_factor, Animal_ID_factor, Group_ID_factor, Cell_ID_factor, Radius)

#manually relabel animal ID so everything is ordered
unq_anID = unique(retina_dat$Animal_ID_factor)
Animal_ID_factor_new = c()
for(i in seq_along(unq_anID)){
    Animal_ID_factor_new = c(Animal_ID_factor_new, rep(i, sum(retina_dat$Animal_ID_factor == unq_anID[i])))
}
retina_dat = retina_dat %>% mutate(Animal_ID_factor_new = Animal_ID_factor_new)

########################################################################
### Restructuring data for jags
########################################################################

### grouping data according to desired levels of hierarchy

retina_dat_obs_level = retina_dat %>% select(Radius, Inters., Cell_ID_factor)
N = nrow(retina_dat_obs_level)

retina_dat_cell_level =  retina_dat %>% distinct(Cell_ID_factor, .keep_all = TRUE) %>% select(Cell_ID_factor, Animal_ID_factor_new, Group_ID_factor, Genotype_ID_factor)
N.cell = nrow(retina_dat_cell_level)

retina_dat_group_level = retina_dat %>% distinct(Group_ID_factor)
N.group = nrow(retina_dat_group_level)

retina_dat_animal_level = retina_dat %>% distinct(Animal_ID, .keep_all = TRUE) %>% select(Animal_ID, Genotype_ID_factor)
N.animal = nrow(retina_dat_animal_level)

retina_dat_genotype_level = retina_dat %>% distinct(Genotype_ID_factor)
N.genotype = nrow(retina_dat_genotype_level)

Cell_ID = as.numeric(as.factor(retina_dat_obs_level$Cell_ID_factor))
Animal_ID = as.numeric(as.factor(retina_dat_cell_level$Animal_ID_factor_new))
Group_ID = as.numeric(as.factor(retina_dat_cell_level$Group_ID_factor)) - 1 # minus 1 to make is 0-1 coding
Genotype_ID_4Int = 1 - (as.numeric(as.factor(retina_dat_cell_level$Genotype_ID_factor)) - 1)  # minus 1 to make is 0-1 coding
Genotype_ID = 3 - as.numeric(as.factor(retina_dat_animal_level$Genotype_ID_factor))

d = list(y = retina_dat_obs_level$Inters.,
         x = retina_dat_obs_level$Radius,
         Cell_ID = Cell_ID,
         Animal_ID = Animal_ID,
         Genotype_ID = Genotype_ID,
         Genotype_ID_4Int = Genotype_ID_4Int,
         
         Group_ID = Group_ID,
         N = N,
         N.Cell = N.cell,
         # N.Group = N.group,
         N.Animal = N.animal,
         N.Genotype = N.genotype)



################################################
### Run Jags Sampler
################################################

#initialize directories and parameters for sampling
n.chains = 4
seed = 92370243
set.seed(seed)
seeds = sample(1e8, n.chains)
n.adapt = 10000
n.burn = 250000
n.samp = 500000
n.thin = 50


# Run chains in parallel
snow.start.time = proc.time()
cl <- makeCluster(n.chains)

##Make sure the rjags library is loaded in each worker
clusterEvalQ(cl, library(rjags))

##Send data to workers, then fit models. One disadvantage of this
##parallelization is that you lose the ability to watch the progress bar.
clusterExport(cl, list("d", "inits_func", "n.adapt", "n.burn", "n.samp", "m", "seeds", "n.thin"))

snow.start.time = proc.time()
samp = clusterApply(cl, 1:n.chains, rjags_samp_wrapper)
snow.end.time = proc.time()
snow.dtime = snow.end.time - snow.start.time

stopCluster(cl)

#resample from the posterior to enforce constraints

cat(getwd())

saveRDS(samp, "./Out/Fits/crushKO_fit.RDS")

# saveRDS(samp, "Output Data/MCMC Samples/crushcontrol_homoskedastic.RDS")

########################################################################
### Formal Analysis
########################################################################


effect_names = names(samp[[1]])[grepl("effect", names(samp[[1]]))]
effect_names = effect_names[!grepl("priorMean", effect_names)]
effect_names = effect_names[!grepl("priorSD", effect_names)]

mean(c(samp[[1]][[effect_names[1]]]) < 0)
mean(c(samp[[1]][[effect_names[1]]]) > 0)


effect_inf_mat = matrix(nrow = length(effect_names), ncol = 2)
rownames(effect_inf_mat) = effect_names
colnames(effect_inf_mat) = c("Prop < 0", "Prop > 0")

for(i in seq_along(effect_names)){
    post_samples = c()
    for(j in seq_along(samp)){
        post_samples = c(post_samples, c(samp[[j]][[effect_names[i]]]))
    }
    effect_inf_mat[i, 1] = mean(post_samples < 0)
    effect_inf_mat[i, 2] = mean(post_samples > 0)
}
effect_inf_mat


########################################################################
### Create Plots
########################################################################

#function to reformat mcmc samples run in parallel to rjags output when runin sequence
fits2df = function(samp){
    #initialize dataframe for fitted model
    fitted_dat = data.frame()
    
    #loop over parallel samples
    for(i in seq_along(samp)){
        
        #check if chain is NA (error in jags sampler, likely infinite density)
        if(any(is.na(samp[[i]]))) next
        
        #initialize sub-sample and sub-label vectors
        samp_vec = c()
        label_vec = c()
        num_vec = c()
        
        #loop over parameters
        for(p in seq_along(samp[[i]])){
            
            #check if only 1 chains (there should always only be 1 chain)
            if(seq_len(dim(samp[[i]][[p]])[3]) != 1) stop(paste0("Parallel sample ", i, ", parameter ", p, ", does not have exactly 1 chain."))
            
            #initialize sub-sample and sub-label vectors
            sub_samp_vec = c()
            sub_label_vec = c()
            sub_num_vec = c()
            
            #loop over parameter vector
            for(v in seq_len(dim(samp[[i]][[p]])[1])){
                
                #concatenate samples for this element of the parameter vector
                sub_samp_vec = c(sub_samp_vec, samp[[i]][[p]][v,,1])
                
                #concatenate label
                sub_label_vec = c(sub_label_vec, rep(names(samp[[i]])[p], dim(samp[[i]][[p]])[2]))
                sub_num_vec = c(sub_num_vec, rep(v, dim(samp[[i]][[p]])[2]))
            }
            
            samp_vec = c(samp_vec, sub_samp_vec)
            label_vec = c(label_vec, sub_label_vec)
            num_vec = c(num_vec, sub_num_vec)
        }
        fitted_dat = rbind(fitted_dat, data.frame(samp = samp_vec, par = label_vec, par_num = num_vec, chain = i))
    }
    return(fitted_dat)
}

#reformat data
samp = fits2df(samp)



#Fitted Curve Plot

par_dat = samp %>% group_by(par, par_num) %>% summarise(samp = mean(samp))

nn = 200
xx = seq(0, 100, length.out = 300)
yy = c()
n.cell = par_dat %>% filter(par == "a1") %>% pull(par_num) %>% unique
for(i in seq_along(n.cell)){
    a1_tmp = par_dat %>% filter(par == "a1") %>% pull(samp) %>% .[i]
    a2_tmp = par_dat %>% filter(par == "a2") %>% pull(samp) %>% .[i]
    g_tmp = par_dat %>% filter(par == "g") %>% pull(samp) %>% .[i]
    t_tmp = par_dat %>% filter(par == "t") %>% pull(samp) %>% .[i]
    
    
    yy = c(yy, exp(ifelse(xx < g_tmp, 
                          a1_tmp * (g_tmp - xx)^2 + t_tmp, 
                          a2_tmp * (xx - g_tmp)^2 + t_tmp)))
}
ggdat = data.frame(x = rep(xx, length(n.cell)), y = yy, id = rep(seq_along(n.cell), each = 300), Genotype_ID = rep(d$Genotype_ID_4Int, each = 300), Group_ID = rep(d$Group_ID, each = 300), Animal_ID_factor_new = rep(d$Animal_ID, each = 300))
ggdat = ggdat %>% left_join(., retina_dat %>% select(Animal_ID, Animal_ID_factor_new) %>% distinct, by = "Animal_ID_factor_new")

ggdat = ggdat %>% mutate(Genotype = ifelse(Genotype_ID == 1, "KO", "WT"), Group = ifelse(Group_ID == 1, "Crush", "Control"), `Genotype/Group` = paste(Genotype, Group, sep = "/"))



cell_level_facet_plot = ggplot(ggdat, aes(x = x, y = y, color = `Genotype/Group`, group = as.factor(id)))+
    geom_line(size = 0.25, alpha = 0.8) + 
    # facet_wrap(~Animal_ID, nrow = 2) + 
    facet_wrap(~factor(Animal_ID, levels = c("59873", "59875", "59879", "59976", "59876", "59877", "60043", "60044")), nrow = 2) +
    theme(legend.position = "bottom") + 
    ylab("Number of Intersections") + 
    xlab("Radius") + 
    guides(color = guide_legend(title = "Genotype / Group", override.aes = list(size = 1))) +
    scale_color_manual(values = c("#B8DE29FF", "#2D708EFF", "#29AF7FFF", "#482677FF")) + 
    ggtitle("Cell Level Fitted Curves by Animal") + 
    theme(text = element_text(size = 14))
ggsave("./Out/Plots/CrushKO/cell_level_facet_NewData.png", cell_level_facet_plot, width = 10, height = 7)


### Effect Plot

ggdat = samp %>% filter(grepl("\\_effect", par))  %>% filter(!grepl("priorMean", par)) %>% filter(!grepl("priorSD", par))

ggdat = ggdat %>% 
    mutate(par_eff = par) %>%
    group_by(par_eff) %>% 
    mutate(par = strsplit(par_eff, split = "\\.")[[1]][1], 
           effect = strsplit(strsplit(par_eff, split = "\\.")[[1]][2], split = "_")[[1]][1])

ggdat$par = factor(ggdat$par, labels = c("A" = parse(text = TeX('$b_{\\alpha_1}')),
                                         "B" = parse(text = TeX('$b_{\\alpha_2}$')),
                                         "C" = parse(text = TeX('$b_{\\gamma}$')),
                                         "D" = parse(text = TeX('$b_{\\tau}$'))))

ggdat = ggdat %>%
    mutate(effect = ifelse(effect == "crush", "Condition", effect)) %>%
    mutate(effect = ifelse(effect == "ko", "Genotype", effect)) %>%
    mutate(effect = ifelse(effect == "int", "Interaction", effect))

ggdat_summ = ggdat %>% 
    group_by(par, effect) %>% 
    summarise(mean = mean(samp), 
              LB = hdi(samp)[1], 
              UB = hdi(samp)[2]) 



posterior_interval_plot = ggplot(ggdat_summ, 
                                 aes(x = mean, y = as.factor(effect))) + 
    geom_violin(data = ggdat, 
                aes(x = samp, y = as.factor(effect), 
                    fill = as.factor(effect), 
                    color = as.factor(effect)), 
                alpha = 0.5, 
                show.legend = F) + 
    geom_vline(xintercept = 0, 
               linetype = 'dotted', 
               color = 'red') + 
    geom_point() +
    geom_text(data = ggdat_summ, aes(label = ifelse(abs(mean) < 0.01, round(mean, 6), round(mean, 2))), nudge_y = 0.16, size = 2.5) +
    geom_errorbar(aes(xmin = LB, xmax = UB)) + 
    facet_wrap(~par, scales = 'free', labeller = label_parsed) +
    ylab(NULL) +
    xlab(NULL) + 
    ggtitle("95% Credible Intervals for Effects by Parameter")

ggsave(filename = "./Out/Plots/CrushKO/effects_NewData.png", 
       plot = posterior_interval_plot, 
       width = 8.5, 
       height = 5.2)






