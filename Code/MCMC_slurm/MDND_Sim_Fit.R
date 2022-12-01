library(truncnorm)
library(parallel)
library(dplyr)
library(ggplot2)
library(purrr)
library(rjags)
library(LaplacesDemon)


setwd("/gpfs/fs2/scratch/evonkaen/Microglia/Sholl Analysis")

scen = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(scen)

# manually initialize true effects
blank_eff_mat = matrix(0, nrow = 3, ncol = 4, dimnames = list(c("true_cond_effect", "true_side_effect", "true_int_effect"), c("a1.cond_side", "a2.cond_side", "g.cond_side", "t.cond_side")))

scen1_effect_mat = blank_eff_mat
scen1_effect_mat[1,4] = 0.5

scen2_effect_mat = blank_eff_mat
scen2_effect_mat[2,4] = -0.25

scen3_effect_mat = blank_eff_mat
scen3_effect_mat[1,4] = 0.5
scen3_effect_mat[2,4] = -0.25

scen4_effect_mat = blank_eff_mat
scen4_effect_mat[1,4] = 0.5
scen4_effect_mat[2,4] = -0.25
scen4_effect_mat[3,4] = 0.5

scen5_effect_mat = scen4_effect_mat

scen6_effect_mat = scen4_effect_mat

scen7_effect_mat = blank_eff_mat

true_effect_list = list(scen1_effect_mat,
                        scen2_effect_mat,
                        scen3_effect_mat,
                        scen4_effect_mat,
                        scen5_effect_mat,
                        scen6_effect_mat,
                        scen7_effect_mat)

true_effect_mat = true_effect_list[[scen]]

m = "model{
        for(i in 1:N){
            y[i] ~ dpois(lambda[i])
            eta[i] = ifelse(x[i] <= g[Cell_ID[i]], a1[Cell_ID[i]] * (g[Cell_ID[i]] - x[i])^2 + t[Cell_ID[i]], a2[Cell_ID[i]] * (x[i] - g[Cell_ID[i]])^2 + t[Cell_ID[i]])
            log(lambda[i]) <- eta[i]
        }
        for(j in 1:N.Cell){
            a1[j] ~ dnorm(a1.animal[Animal_ID[j]], a1.tau.cell) T(, 0)
            a2[j] ~ dnorm(a2.animal[Animal_ID[j]], a2.tau.cell) T(, 0)
            g[j] ~ dnorm(g.animal[Animal_ID[j]], g.tau.cell) T(0, 98)
            t[j] ~ dnorm(t.animal[Animal_ID[j]], t.tau.cell) T(0,)
        }
        for(k in 1:N.Animal){
            a1.animal[k] ~ dnorm(a1.cond_side[Cond_Side_ID[k]], a1.tau.animal) T(, 0)
            a2.animal[k] ~ dnorm(a2.cond_side[Cond_Side_ID[k]], a2.tau.animal) T(, 0)
            g.animal[k] ~ dnorm(g.cond_side[Cond_Side_ID[k]], g.tau.animal) T(0, 98)
            t.animal[k] ~ dnorm(t.cond_side[Cond_Side_ID[k]], t.tau.animal) T(0,)
        }
        ## Group Level Stuff
            ### Easier/better visualization to hard-code
            a1.cond_side[1] ~ dnorm(a1.pop, a1.tau.cond_side) T(, 0)
            a1.cond_side[2] ~ dnorm(a1.pop + a1.side_effect, a1.tau.cond_side) T(, 0)
            a1.cond_side[3] ~ dnorm(a1.pop + a1.cond_effect, a1.tau.cond_side) T(, 0)
            a1.cond_side[4] ~ dnorm(a1.pop + a1.cond_effect + a1.side_effect + a1.int_effect, a1.tau.cond_side) T(, 0)
            
            a2.cond_side[1] ~ dnorm(a2.pop, a2.tau.cond_side) T(, 0)
            a2.cond_side[2] ~ dnorm(a2.pop + a2.side_effect, a2.tau.cond_side) T(, 0)
            a2.cond_side[3] ~ dnorm(a2.pop + a2.cond_effect, a2.tau.cond_side) T(, 0)
            a2.cond_side[4] ~ dnorm(a2.pop + a2.cond_effect + a2.side_effect + a2.int_effect, a2.tau.cond_side) T(, 0)
            
            g.cond_side[1] ~ dnorm(g.pop, g.tau.cond_side) T(0, 98)
            g.cond_side[2] ~ dnorm(g.pop + g.side_effect, g.tau.cond_side) T(0, 98)
            g.cond_side[3] ~ dnorm(g.pop + g.cond_effect, g.tau.cond_side) T(0, 98)
            g.cond_side[4] ~ dnorm(g.pop + g.cond_effect + g.side_effect + g.int_effect, g.tau.cond_side) T(0, 98)
            
            t.cond_side[1] ~ dnorm(t.pop, t.tau.cond_side) T(0, )
            t.cond_side[2] ~ dnorm(t.pop + t.side_effect, t.tau.cond_side) T(0, )
            t.cond_side[3] ~ dnorm(t.pop + t.cond_effect, t.tau.cond_side) T(0, )
            t.cond_side[4] ~ dnorm(t.pop + t.cond_effect + t.side_effect + t.int_effect, t.tau.cond_side) T(0, )
        
       ### Effects 
       
       ## a1
       a1.cond_effect ~ dnorm(a1.cond_effect_priorMean, a1.cond_effect_priorTau) T(, -a1.pop)
       a1.cond_effect_priorMean ~ dnorm(0, 70000000) T(, -a1.pop)
       a1.cond_effect_priorTau <- pow(a1.cond_effect_priorSD, -2)
       a1.cond_effect_priorSD ~ dt(0, 100000000, 4) T(0,)
       
       a1.side_effect ~ dnorm(a1.side_effect_priorMean, a1.side_effect_priorTau) T(, -a1.pop)
       a1.side_effect_priorMean ~ dnorm(0, 70000000) T(, -a1.pop)
       a1.side_effect_priorTau <- pow(a1.side_effect_priorSD, -2)
       a1.side_effect_priorSD ~ dt(0, 100000000, 4) T(0,)
       
       a1.int_effect ~ dnorm(a1.int_effect_priorMean, a1.int_effect_priorTau) T(, -(a1.pop + a1.cond_effect + a1.side_effect))
       a1.int_effect_priorMean ~ dnorm(0, 70000000) T(, -(a1.pop + a1.cond_effect + a1.side_effect))
       a1.int_effect_priorTau <- pow(a1.int_effect_priorSD, -2)
       a1.int_effect_priorSD ~ dt(0, 100000000, 4) T(0,)
       
       ## a2
       a2.cond_effect ~ dnorm(a2.cond_effect_priorMean, a2.cond_effect_priorTau) T(, -a2.pop)
       a2.cond_effect_priorMean ~ dnorm(0, 70000000) T(, -a2.pop)
       a2.cond_effect_priorTau <- pow(a2.cond_effect_priorSD, -2)
       a2.cond_effect_priorSD ~ dt(0, 100000000, 4) T(0,)
       
       a2.side_effect ~ dnorm(a2.side_effect_priorMean, a2.side_effect_priorTau) T(, -a2.pop)
       a2.side_effect_priorMean ~ dnorm(0, 70000000) T(, -a2.pop)
       a2.side_effect_priorTau <- pow(a2.side_effect_priorSD, -2)
       a2.side_effect_priorSD ~ dt(0, 100000000, 4) T(0,)
       
       a2.int_effect ~ dnorm(a2.int_effect_priorMean, a2.int_effect_priorTau) T(, -(a2.pop + a2.cond_effect + a2.side_effect))
       a2.int_effect_priorMean ~ dnorm(0, 70000000) T(, -(a2.pop + a2.cond_effect + a2.side_effect))
       a2.int_effect_priorTau <- pow(a2.int_effect_priorSD, -2)
       a2.int_effect_priorSD ~ dt(0, 100000000, 4) T(0,)
       
       ## g
       g.cond_effect ~ dnorm(g.cond_effect_priorMean, g.cond_effect_priorTau) T( -g.pop, )
       g.cond_effect_priorMean ~ dnorm(0, 0.01) T( -g.pop, )
       g.cond_effect_priorTau <- pow(g.cond_effect_priorSD, -2)
       g.cond_effect_priorSD ~ dt(0, 0.1, 4) T(0,)
       
       g.side_effect ~ dnorm(g.side_effect_priorMean, g.side_effect_priorTau) T( -g.pop, )
       g.side_effect_priorMean ~ dnorm(0, 0.01) T( -g.pop, )
       g.side_effect_priorTau <- pow(g.side_effect_priorSD, -2)
       g.side_effect_priorSD ~ dt(0, 0.1, 4) T(0,)
       
       g.int_effect ~ dnorm(g.int_effect_priorMean, g.int_effect_priorTau)  T(-(g.pop + g.cond_effect + g.side_effect), )
       g.int_effect_priorMean ~ dnorm(0, 0.01) T(-(g.pop + g.cond_effect + g.side_effect), )
       g.int_effect_priorTau <- pow(g.int_effect_priorSD, -2)
       g.int_effect_priorSD ~ dt(0, 0.1, 4) T(0,)
       
       ## t
       t.cond_effect ~ dnorm(t.cond_effect_priorMean, t.cond_effect_priorTau) T( -t.pop,)
       t.cond_effect_priorMean ~ dnorm(0, 3) T( -t.pop,)
       t.cond_effect_priorTau <- pow(t.cond_effect_priorSD, -2)
       t.cond_effect_priorSD ~ dt(0, 70, 4) T(0,)
       
       t.side_effect ~ dnorm(t.side_effect_priorMean, t.side_effect_priorTau) T( -t.pop,)
       t.side_effect_priorMean ~ dnorm(0, 3) T( -t.pop,)
       t.side_effect_priorTau <- pow(t.side_effect_priorSD, -2)
       t.side_effect_priorSD ~ dt(0, 70, 4) T(0,)
       
       t.int_effect ~ dnorm(t.int_effect_priorMean, t.int_effect_priorTau)  T(-(t.pop + t.cond_effect + t.side_effect), )
       t.int_effect_priorMean ~ dnorm(0, 3) T(-(t.pop + t.cond_effect + t.side_effect), )
       t.int_effect_priorTau <- pow(t.int_effect_priorSD, -2)
       t.int_effect_priorSD ~ dt(0, 70, 4) T(0,)
            
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
        
        #cond_side level stuff
        a1.tau.cond_side <- pow(a1.sd.cond_side, -2)
        a2.tau.cond_side <- pow(a2.sd.cond_side, -2)
        g.tau.cond_side <- pow(g.sd.cond_side, -2)
        t.tau.cond_side <- pow(t.sd.cond_side, -2)
        
        a1.sd.cond_side ~ dt(0, 1500000000, 4) T(0, ) # mean, precision, df
        a2.sd.cond_side ~ dt(0, 1500000000, 4) T(0, )
        g.sd.cond_side ~ dt(0, 3, 4) T(0, )
        t.sd.cond_side ~ dt(0, 1500, 4) T(0, )

        #pop level stuff
        a1.pop ~ dnorm(0, 100000) T(, 0)
        a2.pop ~ dnorm(0, 100000) T(, 0)
        g.pop ~ dnorm(0, 0.01) T(0, 100)
        t.pop ~ dnorm(0, 0.25) T(0,)
    }"

inits_func = function(chain, N.Animal, N.Cond_Side, N.Cell, true_effect_mat){
    gen_list = function(chain = chain){
        out = list(
            a1 = truncnorm::rtruncnorm(N.Cell, b = 0, mean = 0.002, sd = 0.00001),
            a2 = truncnorm::rtruncnorm(N.Cell, b = 0, mean = 0.002, sd = 0.00001),
            t = truncnorm::rtruncnorm(N.Cell, a = 0, mean = 2, sd = 0.01),
            g = truncnorm::rtruncnorm(N.Cell, a = 0, b = 100, mean = 30, sd = 0.1)
        )
        
        out[['a1.cond_effect']] = true_effect_mat[1,1]
        out[['a2.cond_effect']] = true_effect_mat[1,2]
        out[['t.cond_effect']] = true_effect_mat[1,4]
        out[['g.cond_effect']] = true_effect_mat[1,3]
        
        out[['a1.side_effect']] = true_effect_mat[2,1]
        out[['a2.side_effect']] = true_effect_mat[2,2]
        out[['t.side_effect']] = true_effect_mat[2,4]
        out[['g.side_effect']] = true_effect_mat[2,3]
        
        out[['a1.int_effect']] = true_effect_mat[3,1]
        out[['a2.int_effect']] = true_effect_mat[3,2]
        out[['t.int_effect']] = true_effect_mat[3,4]
        out[['g.int_effect']] = true_effect_mat[3,3]
        
        
        out[['a1.cond_effect_priorMean']] = 0
        out[['a2.cond_effect_priorMean']] = 0
        out[['t.cond_effect_priorMean']] = 0
        out[['g.cond_effect_priorMean']] = 0
        
        out[['a1.side_effect_priorMean']] = 0
        out[['a2.side_effect_priorMean']] = 0
        out[['t.side_effect_priorMean']] = 0
        out[['g.side_effect_priorMean']] = 0
        
        out[['a1.int_effect_priorMean']] = 0
        out[['a2.int_effect_priorMean']] = 0
        out[['t.int_effect_priorMean']] = 0
        out[['g.int_effect_priorMean']] = 0
        
        
        out[['a1.cond_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['a2.cond_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['t.cond_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.1, nu = 4)
        out[['g.cond_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 1, nu = 4)
        
        out[['a1.side_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['a2.side_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['t.side_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.1, nu = 4)
        out[['g.side_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 1, nu = 4)
        
        out[['a1.int_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['a2.int_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['t.int_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 0.1, nu = 4)
        out[['g.int_effect_priorSD']] = LaplacesDemon::rhalft(1, scale = 1, nu = 4)


        out[['a1.pop']] = mean(out$a1)
        out[['a2.pop']] = mean(out$a2)
        out[['g.pop']] = mean(out$g)
        out[['t.pop']] = mean(out$t)
        
        out[['a1.sd.cond_side']] = LaplacesDemon::rhalft(1, scale = 0.000025, nu = 4)
        out[['a2.sd.cond_side']] = LaplacesDemon::rhalft(1, scale = 0.000025, nu = 4)
        out[['g.sd.cond_side']] = LaplacesDemon::rhalft(1, scale = 0.5, nu = 4)
        out[['t.sd.cond_side']] = LaplacesDemon::rhalft(1, scale = 0.025, nu = 4)
        
        out[['a1.cond_side']] = rep(mean(out$a1), N.Cond_Side)
        out[['a2.cond_side']] = rep(mean(out$a2), N.Cond_Side)
        out[['g.cond_side']] = rep(mean(out$g), N.Cond_Side)
        out[['t.cond_side']] = rep(mean(out$t), N.Cond_Side)
        
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


#function to collapse simulated data at animal level
collapse_at_animal = function(d){
  avg_curve = c()
  Animal_ID_aggregate = c()
  x_avg = c()
  for(i in seq_len(d$N.Animal)){
    bool_ind = d$Cell_ID %in% which(d$Animal_ID == i)
    y_tmp = d$y[bool_ind]
    id_tmp = d$Cell_ID[bool_ind]
    sum_curve = rep(0, 100)
    for(j in unique(id_tmp)){
      sum_curve = sum_curve + y_tmp[which(id_tmp == j)]
    }
    avg_curve = c(avg_curve, sum_curve / length(unique(id_tmp)))
    x_avg = c(x_avg, 1:100)
    Animal_ID_aggregate = c(Animal_ID_aggregate, rep(i, 100))
  }
  d$y = round(avg_curve, 0)
  d$x = x_avg
  d$Animal_ID = Animal_ID_aggregate
  d$N = length(d$y)
  return(d)
}


rjags_samp_wrapper = function(jj){
  set.seed(chain_seeds[jj])
  temp_model = rjags::jags.model(file = textConnection(m),
                                 data = d,
                                 inits = function(chain) inits_func(chain = chain, N.Animal = d$N.Animal, N.Cond_Side = d$N.Cond_Side, N.Cell = d$N.Cell, true_effect_mat = true_effect_mat),
                                 n.adapt = n.adapt,
                                 n.chains = 1)
  update(temp_model, n.burn) #burning in n.burn samples
  out = rjags::jags.samples(model = temp_model, 
                            variable.names = names(inits_func(chain = 1, N.Animal = d$N.Animal, N.Cond_Side = d$N.Cond_Side, N.Cell = d$N.Cell, true_effect_mat = true_effect_mat)), 
                            n.iter = n.samp, 
                            thin = n.thin)
  # saveRDS(out, file = paste0("./MCMC Samples/crushcontrol_hierarchical_reParameterized_chain", j, ".RDS"))
  return(out)
}

#error catch just in case
#should only trigger if caught at infinite density, then returns list of NAs. Has not happened under current formulation.
rjags_samp_wrapper_safe = possibly(rjags_samp_wrapper, NA)
# rjags_samp_wrapper_safe = rjags_samp_wrapper

# Run chains in parallel
dat_dir = "Simulated_Data/MDND_new/"
fit_dir = "Simulated_Fits/MDND/"
n.chains = 4
n.adapt = 5000
n.burn = 15000
n.samp = 15000
n.thin = 10
scenarios = sapply(strsplit(list.files(dat_dir), split = '_'), function(x) x[2])

#create seeds to mcmc samplers
n.datasets = 0
for(i in scen){
  n.datasets = (n.datasets + length(list.files(paste0(dat_dir, "scenario_", scenarios[i]))[grepl("RDS", list.files(paste0(dat_dir, "scenario_", scenarios[i])))]))
}
grandparent_seed = scen
set.seed(grandparent_seed)
mcmc_seeds = sample(seq_len(1e8), n.datasets)

seed_num = 1
for(i in scen){
  # for(i in 1:2){
  cat("Scenario: ", i, "\n")
  datasets = list.files(paste0(dat_dir, "scenario_", scenarios[i]))[grepl("RDS", list.files(paste0(dat_dir, "scenario_", scenarios[i])))]
  for(j in seq_along(datasets)){
    data_seed = strsplit(datasets[j], split = "\\.")[[1]][1]
    # for(j in 1:2){
    cat("\tFitting dataset", j, " out of ", length(datasets), "\n")
    #generate 1 seed for each chain
    set.seed(mcmc_seeds[seed_num])
    chain_seeds = sample(seq_len(1e8), n.chains)
    
    
    #skip dataset if theres already a fit for it in the output directory
    if(!file.exists(paste0(fit_dir, "scenario_", scenarios[i], "/", data_seed, "_fitted.RDS"))){
        #read in data
        out_list = readRDS(paste0(dat_dir, "scenario_", scenarios[i], "/", datasets[j]))
        d = out_list$d
        # d = collapse_at_animal(d)
        
        #make cluster
        snow.start.time = proc.time()
        cl <- makeCluster(n.chains)
        
        ##Make sure the rjags library is loaded in each worker
        clusterEvalQ(cl, library(rjags))
        
        ##Send data to workers, then fit models. One disadvantage of this
        ##parallelization is that you lose the ability to watch the progress bar.
        clusterExport(cl, list("d", "inits_func", "n.adapt", "n.burn", "n.samp", "m", "chain_seeds", "n.thin", "true_effect_mat"))
        
        snow.start.time = proc.time()
        par.samples = clusterApply(cl, 1:n.chains, rjags_samp_wrapper_safe)
        snow.end.time = proc.time()
        snow.dtime = snow.end.time - snow.start.time
        
        #end cluster
        stopCluster(cl)
        
        #save relevant information
        # data_seed = strsplit(datasets[j], split = "\\.")[[1]][1]
        fit_out = list(par.samples = par.samples, time = snow.dtime, chain_seeds = chain_seeds, mcmc_seed = mcmc_seeds[seed_num], data_seed = as.numeric(data_seed))
        if(!dir.exists(paste0(fit_dir, "scenario_", scenarios[i]))) dir.create(paste0(fit_dir, "scenario_", scenarios[i]))
        saveRDS(fit_out, paste0(fit_dir, "scenario_", scenarios[i], "/", data_seed, "_fitted.RDS"))
    }
    
    seed_num = seed_num + 1
  }
}

