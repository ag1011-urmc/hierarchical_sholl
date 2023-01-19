#Script to fit model 2 to mdnd data and reproduce corresponding figures

################################################
### load packages
################################################

library(LaplacesDemon)
library(tidyverse)
library(rjags)
library(truncnorm)
library(ggplot2)
library(cowplot)
library(latex2exp)
library(HDInterval)
library(cowplot)
library(readxl)

#######################################
### Define Model
#######################################

m = "model{
        for(i in 1:N){
            y[i] ~ dpois(lambda[i])
            eta[i] = ifelse(x[i] <= g[Animal_ID[i]], a1[Animal_ID[i]] * (g[Animal_ID[i]] - x[i])^2 + t[Animal_ID[i]], a2[Animal_ID[i]] * (x[i] - g[Animal_ID[i]])^2 + t[Animal_ID[i]])
            log(lambda[i]) <- eta[i]
        }
        for(k in 1:N.Animal){
            a1[k] ~ dnorm(a1.cond_side[Cond_Side_ID[k]], a1.tau.animal) T(, 0)
            a2[k] ~ dnorm(a2.cond_side[Cond_Side_ID[k]], a2.tau.animal) T(, 0)
            g[k] ~ dnorm(g.cond_side[Cond_Side_ID[k]], g.tau.animal) T(0, 98)
            t[k] ~ dnorm(t.cond_side[Cond_Side_ID[k]], t.tau.animal) T(0,)
        }
        ## Group Level Stuff
            ### Easier/better visualization to hard-code
            a1.cond_side[3] ~ dnorm(a1.pop, a1.tau.cond_side) T(, 0)
            a1.cond_side[1] ~ dnorm(a1.pop + a1.side_effect, a1.tau.cond_side) T(, 0)
            a1.cond_side[4] ~ dnorm(a1.pop + a1.cond_effect, a1.tau.cond_side) T(, 0)
            a1.cond_side[2] ~ dnorm(a1.pop + a1.cond_effect + a1.side_effect + a1.int_effect, a1.tau.cond_side) T(, 0)
            
            a2.cond_side[3] ~ dnorm(a2.pop, a2.tau.cond_side) T(, 0)
            a2.cond_side[1] ~ dnorm(a2.pop + a2.side_effect, a2.tau.cond_side) T(, 0)
            a2.cond_side[4] ~ dnorm(a2.pop + a2.cond_effect, a2.tau.cond_side) T(, 0)
            a2.cond_side[2] ~ dnorm(a2.pop + a2.cond_effect + a2.side_effect + a2.int_effect, a2.tau.cond_side) T(, 0)
            
            g.cond_side[3] ~ dnorm(g.pop, g.tau.cond_side) T(0, 98)
            g.cond_side[1] ~ dnorm(g.pop + g.side_effect, g.tau.cond_side) T(0, 98)
            g.cond_side[4] ~ dnorm(g.pop + g.cond_effect, g.tau.cond_side) T(0, 98)
            g.cond_side[2] ~ dnorm(g.pop + g.cond_effect + g.side_effect + g.int_effect, g.tau.cond_side) T(0, 98)
            
            t.cond_side[3] ~ dnorm(t.pop, t.tau.cond_side) T(0, )
            t.cond_side[1] ~ dnorm(t.pop + t.side_effect, t.tau.cond_side) T(0, )
            t.cond_side[4] ~ dnorm(t.pop + t.cond_effect, t.tau.cond_side) T(0, )
            t.cond_side[2] ~ dnorm(t.pop + t.cond_effect + t.side_effect + t.int_effect, t.tau.cond_side) T(0, )
        
       ### Effects 
       
       ## a1
       a1.cond_effect ~ dnorm(a1.cond_effect_priorMean, a1.cond_effect_priorTau) T(, -a1.pop)
       a1.cond_effect_priorMean ~ dnorm(0, 50000000) T(, -a1.pop)
       a1.cond_effect_priorTau <- pow(a1.cond_effect_priorSD, -2)
       a1.cond_effect_priorSD ~ dt(0, 70000000, 4) T(0,)
       
       a1.side_effect ~ dnorm(a1.side_effect_priorMean, a1.side_effect_priorTau) T(, -a1.pop)
       a1.side_effect_priorMean ~ dnorm(0, 50000000) T(, -a1.pop)
       a1.side_effect_priorTau <- pow(a1.side_effect_priorSD, -2)
       a1.side_effect_priorSD ~ dt(0, 70000000, 4) T(0,)
       
       a1.int_effect ~ dnorm(a1.int_effect_priorMean, a1.int_effect_priorTau) T(, -(a1.pop + a1.cond_effect + a1.side_effect))
       a1.int_effect_priorMean ~ dnorm(0, 50000000) T(, -(a1.pop + a1.cond_effect + a1.side_effect))
       a1.int_effect_priorTau <- pow(a1.int_effect_priorSD, -2)
       a1.int_effect_priorSD ~ dt(0, 70000000, 4) T(0,)
       
       ## a2
       a2.cond_effect ~ dnorm(a2.cond_effect_priorMean, a2.cond_effect_priorTau) T(, -a2.pop)
       a2.cond_effect_priorMean ~ dnorm(0, 50000000) T(, -a2.pop)
       a2.cond_effect_priorTau <- pow(a2.cond_effect_priorSD, -2)
       a2.cond_effect_priorSD ~ dt(0, 70000000, 4) T(0,)
       
       a2.side_effect ~ dnorm(a2.side_effect_priorMean, a2.side_effect_priorTau) T(, -a2.pop)
       a2.side_effect_priorMean ~ dnorm(0, 50000000) T(, -a2.pop)
       a2.side_effect_priorTau <- pow(a2.side_effect_priorSD, -2)
       a2.side_effect_priorSD ~ dt(0, 70000000, 4) T(0,)
       
       a2.int_effect ~ dnorm(a2.int_effect_priorMean, a2.int_effect_priorTau) T(, -(a2.pop + a2.cond_effect + a2.side_effect))
       a2.int_effect_priorMean ~ dnorm(0, 50000000) T(, -(a2.pop + a2.cond_effect + a2.side_effect))
       a2.int_effect_priorTau <- pow(a2.int_effect_priorSD, -2)
       a2.int_effect_priorSD ~ dt(0, 70000000, 4) T(0,)
       
       ## g
       g.cond_effect ~ dnorm(g.cond_effect_priorMean, g.cond_effect_priorTau) T( -g.pop, )
       g.cond_effect_priorMean ~ dnorm(0, 0.007) T( -g.pop, )
       g.cond_effect_priorTau <- pow(g.cond_effect_priorSD, -2)
       g.cond_effect_priorSD ~ dt(0, 0.07, 4) T(0,)
       
       g.side_effect ~ dnorm(g.side_effect_priorMean, g.side_effect_priorTau) T( -g.pop, )
       g.side_effect_priorMean ~ dnorm(0, 0.007) T( -g.pop, )
       g.side_effect_priorTau <- pow(g.side_effect_priorSD, -2)
       g.side_effect_priorSD ~ dt(0, 0.07, 4) T(0,)
       
       g.int_effect ~ dnorm(g.int_effect_priorMean, g.int_effect_priorTau)  T(-(g.pop + g.cond_effect + g.side_effect), )
       g.int_effect_priorMean ~ dnorm(0, 0.007) T(-(g.pop + g.cond_effect + g.side_effect), )
       g.int_effect_priorTau <- pow(g.int_effect_priorSD, -2)
       g.int_effect_priorSD ~ dt(0, 0.07, 4) T(0,)
       
       ## t
       t.cond_effect ~ dnorm(t.cond_effect_priorMean, t.cond_effect_priorTau) T( -t.pop,)
       t.cond_effect_priorMean ~ dnorm(0, 2) T( -t.pop,)
       t.cond_effect_priorTau <- pow(t.cond_effect_priorSD, -2)
       t.cond_effect_priorSD ~ dt(0, 50, 4) T(0,)
       
       t.side_effect ~ dnorm(t.side_effect_priorMean, t.side_effect_priorTau) T( -t.pop,)
       t.side_effect_priorMean ~ dnorm(0, 2) T( -t.pop,)
       t.side_effect_priorTau <- pow(t.side_effect_priorSD, -2)
       t.side_effect_priorSD ~ dt(0, 50, 4) T(0,)
       
       t.int_effect ~ dnorm(t.int_effect_priorMean, t.int_effect_priorTau)  T(-(t.pop + t.cond_effect + t.side_effect), )
       t.int_effect_priorMean ~ dnorm(0, 2) T(-(t.pop + t.cond_effect + t.side_effect), )
       t.int_effect_priorTau <- pow(t.int_effect_priorSD, -2)
       t.int_effect_priorSD ~ dt(0, 50, 4) T(0,)
       
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

################################################
### Helper functions
################################################

#function to create initial values for jags sampler
inits_func = function(chain, N.Animal, N.Cond_Side){
    gen_list = function(chain = chain){
        out = list(
            a1 = truncnorm::rtruncnorm(N.Animal, b = 0, mean = 0.002, sd = 0.00001),
            a2 = truncnorm::rtruncnorm(N.Animal, b = 0, mean = 0.002, sd = 0.00001),
            t = truncnorm::rtruncnorm(N.Animal, a = 0, mean = 2, sd = 0.01),
            g = truncnorm::rtruncnorm(N.Animal, a = 0, b = 100, mean = 30, sd = 0.1)
        )
        
        out[['a1.cond_effect']] = 0
        out[['a2.cond_effect']] = 0
        out[['t.cond_effect']] = 0
        out[['g.cond_effect']] = 0
        
        out[['a1.side_effect']] = 0
        out[['a2.side_effect']] = 0
        out[['t.side_effect']] = 0
        out[['g.side_effect']] = 0
        
        out[['a1.int_effect']] = 0
        out[['a2.int_effect']] = 0
        out[['t.int_effect']] = 0
        out[['g.int_effect']] = 0
        
        
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
        
        out[['a1.sd.cond_side']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['a2.sd.cond_side']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['g.sd.cond_side']] = LaplacesDemon::rhalft(1, scale = 0.5, nu = 4)
        out[['t.sd.cond_side']] = LaplacesDemon::rhalft(1, scale = 0.025, nu = 4)
        
        out[['a1.cond_side']] = rep(mean(out$a1), N.Cond_Side)
        out[['a2.cond_side']] = rep(mean(out$a2), N.Cond_Side)
        out[['g.cond_side']] = rep(mean(out$g), N.Cond_Side)
        out[['t.cond_side']] = rep(mean(out$t), N.Cond_Side)
        
        out[['a1.sd.animal']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['a2.sd.animal']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['g.sd.animal']] = LaplacesDemon::rhalft(1, scale = 1, nu = 4)
        out[['t.sd.animal']] = LaplacesDemon::rhalft(1, scale = 0.1, nu = 4)
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


#######################################
### Load Data
#######################################

dat = read_excel("./Data/MDND/MajewskaRO1-formatt.xlsx", 1, skip = 2)

intersections = as.vector(as.matrix(dat[ ,-1]))
grps = c("Contra ND", "Contra 12HR MD", "Ipsi ND", "Ipsi 12HR MD")
grps = rep(rep(grps, c(4,4,4,5)), each = 25)
ci = sapply(grps, function(x) strsplit(x, " ")[[1]][1])
cond = sapply(grps, function(x) paste(sapply(strsplit(x, NULL), rev)[2:1], collapse = ""))
reps = c(rep(rep(1:4, each = 25), 4), rep(5, 25))
d = data.frame("cond" = factor(cond, levels = c("ND", "MD")), 
                "lateral" = factor(ci, levels = c("Contra", "Ipsi")),
                "replicate" = reps,
                "R" = rep(dat$`Distance from soma (um)`, 17),
                "N_R"=as.vector(as.matrix(dat[ ,-1])))

rm(ci, cond, grps, intersections, reps, dat)


d = d %>% mutate(Animal_ID = paste0(cond, lateral, replicate))

Animal_ID = c()
ind = 1
for(i in unique(d %>% pull(Animal_ID))){
    Animal_ID = c(Animal_ID, rep(ind, d %>% filter(Animal_ID == i) %>% nrow))
    ind = ind + 1
}
N = nrow(d)
N.Animal = d %>% distinct(replicate, cond, lateral) %>% nrow 
Condition.Dummy = ifelse(d %>% distinct(replicate, cond, lateral) %>% pull(cond) == 'ND', 0, 1)
Side.Dummy = ifelse(d %>% distinct(replicate, cond, lateral) %>% pull(lateral) == 'Contra', 0, 1)
Cond_Side_ID = c(rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 5))

x = d$R
y = d$N_R

dat = list(x = x,
           y = round(y, 0),
           N.Animal = N.Animal,
           Animal_ID = Animal_ID,
           Cond_Side_ID = Cond_Side_ID,
           N = N,
           N.Cond_Side = 4)


################################################
### Run Jags Sampler
################################################

#initialize parameters for sampling
n.adapt = 5000
n.burn = 50000
n.samp = 100000
n.thin = 50
n.chains = 4

set.seed(738645923)
jgs = rjags::jags.model(file = textConnection(m),
                        data = dat,
                        inits = function(chain) inits_func(chain = chain, N.Animal = dat$N.Animal, N.Cond_Side = dat$N.Cond_Side),
                        n.adapt = n.adapt,
                        n.chains = n.chains)

update(jgs, n.burn) #burning in 2000 samples

samp = rjags::jags.samples(jgs, names(inits_func(chain = 1, N.Animal = dat$N.Animal, N.Cond_Side = dat$N.Cond_Side)), n.samp, n.thin)

saveRDS(samp, "Out/Fits/mdnd_fit.RDS")

################################################
### Fitted Curve Plot
################################################

a1_hat = apply(samp$a1.pop, 1, mean)
a2_hat = apply(samp$a2.pop, 1, mean)
g_hat = apply(samp$g.pop, 1, mean)
t_hat = apply(samp$t.pop, 1, mean)

a1.cond_effect = apply(samp$a1.cond_effect, 1, mean)
a2.cond_effect = apply(samp$a2.cond_effect, 1, mean)
g.cond_effect = apply(samp$g.cond_effect, 1, mean)
t.cond_effect = apply(samp$t.cond_effect, 1, mean)

a1.side_effect = apply(samp$a1.side_effect, 1, mean)
a2.side_effect = apply(samp$a2.side_effect, 1, mean)
g.side_effect = apply(samp$g.side_effect, 1, mean)
t.side_effect = apply(samp$t.side_effect, 1, mean)

a1.int_effect = apply(samp$a1.int_effect, 1, mean)
a2.int_effect = apply(samp$a2.int_effect, 1, mean)
g.int_effect = apply(samp$g.int_effect, 1, mean)
t.int_effect = apply(samp$t.int_effect, 1, mean)
#ND==0, Contra==0

comb_mat = matrix(c(0, 0, 0, 1, 1, 0, 1, 1), 4, 2, byrow = T)
rownames(comb_mat) = c("NDIpsi", "NDContra", "MDIpsi", "MDContra")
yy = c()
idid = c()
xx = seq(0, 100, length.out = 300)
for(i in seq_len(nrow(comb_mat))){
    a1_tmp = a1_hat + (a1.cond_effect * comb_mat[i, 1]) + (a1.side_effect * comb_mat[i, 2]) + (a1.int_effect * comb_mat[i, 1] * comb_mat[i, 2])
    a2_tmp = a2_hat + (a2.cond_effect * comb_mat[i, 1]) + (a2.side_effect * comb_mat[i, 2]) + (a2.int_effect * comb_mat[i, 1] * comb_mat[i, 2])
    g_tmp = g_hat + (g.cond_effect * comb_mat[i, 1]) + (g.side_effect * comb_mat[i, 2]) + (g.int_effect * comb_mat[i, 1] * comb_mat[i, 2])
    t_tmp = t_hat + (t.cond_effect * comb_mat[i, 1]) + (t.side_effect * comb_mat[i, 2]) + (t.int_effect * comb_mat[i, 1] * comb_mat[i, 2])
    
    yy = c(yy, exp(ifelse(xx < g_tmp, 
                          a1_hat * (g_tmp - xx)^2 + t_tmp, 
                          a2_hat * (xx - g_tmp)^2 + t_tmp)))
    idid = c(idid, rep(rownames(comb_mat)[i], length(xx)))
}
xx = rep(xx, nrow(comb_mat))

d = d %>% mutate(curve_ID = paste0(cond, lateral))
dfdf = data.frame(y = yy, x = xx, curve_ID = idid) %>% left_join(., (d %>% select(cond, lateral, Animal_ID, curve_ID)), by = "curve_ID")

p2 =ggplot(dfdf, aes(x = x, y = y, group = curve_ID, color = curve_ID)) + 
    geom_line() + 
    geom_point(data = d, aes(x = R, y = N_R, color = curve_ID, group = Animal_ID)) +
    xlim(0, 50) +
    ylab("Number of Intersections") + 
    xlab("Radius") + 
    ggtitle("Animal Level\nFitted Sholl Curves\n")  + 
    guides(color = guide_legend(title = "Condition / Side")) + 
    scale_color_discrete(labels = c("MD / Contralateral", "MD / Ipsilateral", "ND / Contralateral", "ND / Ipsilateral"))+ 
    theme(legend.position = "none") + 
    theme(text = element_text(size = 16))


p1 =ggplot(dfdf, aes(x = x, y = y, group = curve_ID, color = curve_ID)) + 
    geom_line() + 
    geom_point(data = d, aes(x = R, y = N_R, color = curve_ID, group = Animal_ID)) +
    xlim(0, 50) +
    ylab("Number of Intersections") + 
    xlab("Radius") + 
    ggtitle("Animal Level\nFitted Sholl Curves\n")  + 
    guides(color = guide_legend(title = "Condition / Side")) + 
    scale_color_discrete(labels = c("MD / Contralateral", "MD / Ipsilateral", "ND / Contralateral", "ND / Ipsilateral")) +
    facet_wrap(~curve_ID) + 
    theme(legend.position = "bottom") + 
    theme(text = element_text(size = 16))

save_dir = "./Out/Plots/mdnd"

legend <- cowplot::get_legend(p1)

top_row = plot_grid(p1 + theme(legend.position="none"),
                    p2 + theme(legend.position="none"),
                    nrow = 1,
                    labels = c("A", "B")
)


p = plot_grid(top_row,
              legend,
              ncol = 1,
              rel_heights = c(5, 0.2)
)


cowplot::save_plot(paste0(save_dir, '/', 'fitted_plots.png'), p, nrow = 2, ncol = 2, base_height = 3, base_width = 7)

################################################
### Effect Plots
################################################

samp_vec = c()
par_vec = c()
effect_names = names(samp)[grepl("\\_effect", names(samp))]
effect_names = effect_names[!grepl("prior", effect_names)]
for(i in effect_names){
    n_samps = length(c(samp[[i]][1,,]))
    samp_vec = c(samp_vec, c(samp[[i]][1,,]))
    par_vec = c(par_vec, rep(i, n_samps))
}
ggdat = data.frame(samp = samp_vec, par_eff = par_vec)
# ggdat = ggdat %>% mutate(samp = ifelse(par == "t.cond_side", exp(samp), samp))


ggdat = ggdat %>% 
    group_by(par_eff) %>% 
    mutate(par = strsplit(par_eff, split = "\\.")[[1]][1], 
           effect = strsplit(strsplit(par_eff, split = "\\.")[[1]][2], split = "_")[[1]][1])

ggdat$par = factor(ggdat$par, labels = c("A" = parse(text = TeX('$b_{\\alpha_1}')),
                                         "B" = parse(text = TeX('$b_{\\alpha_2}$')),
                                         "C" = parse(text = TeX('$b_{\\gamma}$')),
                                         "D" = parse(text = TeX('$b_{\\tau}$'))))

ggdat = ggdat %>%
    mutate(effect = ifelse(effect == "cond", "Condition", effect)) %>%
    mutate(effect = ifelse(effect == "int", "Interaction", effect)) %>%
    mutate(effect = ifelse(effect == "side", "Side", effect))

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
    geom_text(data = ggdat_summ, aes(label = ifelse(abs(mean) < 0.01, round(mean, 8), round(mean, 2))), nudge_y = 0.16, size = 2.5) +
    geom_errorbar(aes(xmin = LB, xmax = UB)) + 
    facet_wrap(~par, scales = 'free', labeller = label_parsed) +
    ylab(NULL) +
    xlab(NULL) + 
    ggtitle("95% Credible Intervals for Effects by Parameter")


save_dir = "./Out/Plots/mdnd"
ggsave(filename = paste0(save_dir, "/MDND_posterior_interval_plot.png"), 
       plot = posterior_interval_plot, 
       width = 8.5, 
       height = 5.2)


########################################################################
### Posterior probability effect is in the positive and negative directions
########################################################################

effect_names = names(samp)[grepl("effect", names(samp))]
effect_names = effect_names[!grepl("priorMean", effect_names)]
effect_names = effect_names[!grepl("priorSD", effect_names)]

effect_inf_mat = matrix(nrow = length(effect_names), ncol = 2)
rownames(effect_inf_mat) = effect_names
colnames(effect_inf_mat) = c("Prop < 0", "Prop > 0")

for(i in seq_along(effect_names)){
    post_samples = c(samp[[effect_names[i]]])
    effect_inf_mat[i, 1] = mean(post_samples < 0)
    effect_inf_mat[i, 2] = mean(post_samples > 0)
}

effect_inf_mat


