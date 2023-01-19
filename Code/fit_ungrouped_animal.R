#Script to fit model 1 to ungrouped animal data and reproduce corresponding figures

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


#######################################
### Define Model
#######################################

h = "model{
        for(i in 1:N){
            y[i] ~ dpois(lambda[i])
            eta[i] = ifelse(x[i] <= g[Cell_ID[i]], a1[Cell_ID[i]] * (g[Cell_ID[i]] - x[i])^2 + t[Cell_ID[i]], a2[Cell_ID[i]] * (x[i] - g[Cell_ID[i]])^2 + t[Cell_ID[i]])
            log(lambda[i]) <- eta[i]
        }
        for(j in 1:N.Cell){
            a1[j] ~ dnorm(a1.image[Image_ID[j]], a1.tau.cell) T(, 0)
            a2[j] ~ dnorm(a2.image[Image_ID[j]], a2.tau.cell) T(, 0)
            g[j] ~ dnorm(g.image[Image_ID[j]], g.tau.cell) T(0, 98)
            t[j] ~ dnorm(t.image[Image_ID[j]], t.tau.cell) T(0,)
        }
        for(k in 1:N.Image){
            a1.image[k] ~ dnorm(a1.animal[Animal_ID[k]], a1.tau.image) T(, 0)
            a2.image[k] ~ dnorm(a2.animal[Animal_ID[k]], a2.tau.image) T(, 0)
            g.image[k] ~ dnorm(g.animal[Animal_ID[k]], g.tau.image) T(0, 98)
            t.image[k] ~ dnorm(t.animal[Animal_ID[k]], t.tau.image) T(0,)
        }
        for(l in 1:N.Animal){
            a1.animal[l] ~ dnorm(a1.pop, a1.tau.animal) T(, 0)
            a2.animal[l] ~ dnorm(a2.pop, a2.tau.animal) T(, 0)
            g.animal[l] ~ dnorm(g.pop, g.tau.animal) T(0, 98)
            t.animal[l] ~ dnorm(t.pop, t.tau.animal) T(0,)
        }
        
        #cell level
        a1.tau.cell <- pow(a1.sd.cell, -2)
        a2.tau.cell <- pow(a2.sd.cell, -2)
        g.tau.cell <- pow(g.sd.cell, -2)
        t.tau.cell <- pow(t.sd.cell, -2)
        
        a1.sd.cell ~ dt(0, 5000000, 4) T(0, ) # mean, precision, df
        a2.sd.cell ~ dt(0, 5000000, 4) T(0, )
        g.sd.cell ~ dt(0, 0.5, 4) T(0, )
        t.sd.cell ~ dt(0, 50, 4) T(0, )
        
        #image level
        a1.tau.image <- pow(a1.sd.image, -2)
        a2.tau.image <- pow(a2.sd.image, -2)
        g.tau.image <- pow(g.sd.image, -2)
        t.tau.image <- pow(t.sd.image, -2)
        
        a1.sd.image ~ dt(0, 5000000, 4) T(0, ) # mean, precision, df
        a2.sd.image ~ dt(0, 5000000, 4) T(0, )
        g.sd.image ~ dt(0, 0.5, 4) T(0, )
        t.sd.image ~ dt(0, 50, 4) T(0, )
        
        #animal level
        a1.tau.animal <- pow(a1.sd.animal, -2)
        a2.tau.animal <- pow(a2.sd.animal, -2)
        g.tau.animal <- pow(g.sd.animal, -2)
        t.tau.animal <- pow(t.sd.animal, -2)
        
        a1.sd.animal ~ dt(0, 5000000, 4) T(0, ) # mean, precision, df
        a2.sd.animal ~ dt(0, 5000000, 4) T(0, )
        g.sd.animal ~ dt(0, 0.5, 4) T(0, )
        t.sd.animal ~ dt(0, 50, 4) T(0, )
        
        #pop level
        a1.pop ~ dnorm(0, 100000) T(, 0)
        a2.pop ~ dnorm(0, 100000) T(, 0)
        g.pop ~ dnorm(0, 0.01) T(0, 98)
        t.pop ~ dnorm(0, 0.25) T(0,)
    }"


################################################
### Helper functions
################################################

#function to create initial values for jags sampler
inits_func = function(chain, N.Animal, N.Image, N.Cell){
    gen_list = function(chain = chain){
        out = list(
            a1 = truncnorm::rtruncnorm(N.Cell, b = 0, mean = 0.002, sd = 0.00001),
            a2 = truncnorm::rtruncnorm(N.Cell, b = 0, mean = 0.002, sd = 0.00001),
            t = truncnorm::rtruncnorm(N.Cell, a = 0, mean = 3, sd = 0.01),
            g = truncnorm::rtruncnorm(N.Cell, a = 0, b = 100, mean = 20, sd = 0.1)
        )
        out[['a1.pop']] = mean(out$a1)
        out[['a2.pop']] = mean(out$a2)
        out[['g.pop']] = mean(out$g)
        out[['t.pop']] = mean(out$t)
        
        out[['a1.sd.cell']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['a2.sd.cell']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['g.sd.cell']] = LaplacesDemon::rhalft(1, scale = 1, nu = 4)
        out[['t.sd.cell']] = LaplacesDemon::rhalft(1, scale = 0.1, nu = 4)
        
        out[['a1.image']] = rep(mean(out$a1), N.Image)
        out[['a2.image']] = rep(mean(out$a2), N.Image)
        out[['g.image']] = rep(mean(out$g), N.Image)
        out[['t.image']] = rep(mean(out$t), N.Image)
        
        out[['a1.sd.image']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['a2.sd.image']] = LaplacesDemon::rhalft(1, scale = 0.0001, nu = 4)
        out[['g.sd.image']] = LaplacesDemon::rhalft(1, scale = 1, nu = 4)
        out[['t.sd.image']] = LaplacesDemon::rhalft(1, scale = 0.1, nu = 4)
        
        out[['a1.animal']] = rep(mean(out$a1), N.Animal)
        out[['a2.animal']] = rep(mean(out$a2), N.Animal)
        out[['g.animal']] = rep(mean(out$g), N.Animal)
        out[['t.animal']] = rep(mean(out$t), N.Animal)
        
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

root_dir = "./Data/ungrouped_animal"
sholl_data = read.csv(file = paste(root_dir, 
                                   '09 27 2019 New Sholl Analysis.csv', 
                                   sep = '/'), 
                      header = TRUE, 
                      sep = ',')

#clean up some columns
sholl_data = sholl_data %>% 
    mutate(Name = trimws(Name)) %>%
    mutate(Location = trimws(Location))

#add image ID
sholl_data = sholl_data %>% 
    mutate(Cell_ID = paste(Name, 
                           Location, 
                           Number, 
                           sep = '_'),
           Image_ID = paste(Name, 
                            Location,
                            sep = '_'),
           Animal_ID = Name)

#sort data
sholl_data = sholl_data %>% 
    arrange(Name, Location, Number, R)


#remove all but 3 animals
sholl_data = sholl_data %>% filter(Name %in% unique(sholl_data$Name))


#create subsets of dataset corresponding to levels of hierarchy
sholl_data_obs_level = sholl_data %>% select(R, N_R, Cell_ID)
N = nrow(sholl_data_obs_level)

sholl_data_cell_level = sholl_data %>% distinct(Cell_ID, .keep_all = TRUE) %>% select(Cell_ID, Image_ID)
N.Cell = nrow(sholl_data_cell_level)

sholl_data_image_level = sholl_data %>% distinct(Image_ID, .keep_all = TRUE) %>% select(Image_ID, Animal_ID)
N.Image = nrow(sholl_data_image_level)

sholl_data_animal_level = sholl_data %>% distinct(Animal_ID)
N.Animal = nrow(sholl_data_animal_level)

d = list(y = sholl_data_obs_level$N_R,
         x = sholl_data_obs_level$R,
         Cell_ID = as.numeric(as.factor(sholl_data_obs_level$Cell_ID)),
         Image_ID = as.numeric(as.factor(sholl_data_cell_level$Image_ID)),
         Animal_ID = as.numeric(as.factor(sholl_data_image_level$Animal_ID)),
         N = N,
         N.Cell = N.Cell,
         N.Image = N.Image,
         N.Animal = N.Animal
)

#######################################
### Fit Model
#######################################

#initialize directories and parameters for sampling
set.seed(2349812)
n_adapt = 5000
n_burn = 50000
n_samp = 150000
n_thin = 50

jgs = rjags::jags.model(file = textConnection(h),
                        data = d,
                        inits = function(chain) inits_func(chain = chain, N.Cell = d$N.Cell, N.Image = d$N.Image, N.Animal = d$N.Animal),
                        n.adapt = n_adapt,
                        n.chains = 4)

update(jgs, n_burn) #burn in

samp = rjags::jags.samples(jgs, c(names(inits_func(1, N.Cell, N.Image, N.Animal))), n.iter =  n_samp, thin =  n_thin)

saveRDS(samp, file = "./Out/Fits/ungroupe_animal_fit.RDS")


#######################################
### Plots
#######################################

### Fitted curve plots

# Cell Level

a1_hat = apply(samp$a1, 1, mean)
a2_hat = apply(samp$a2, 1, mean)
g_hat = apply(samp$g, 1, mean)
t_hat = apply(samp$t, 1, mean)

nn = 1000
maxX = 100
xx = seq(0, maxX, length.out = nn)
yy = c()
for(i in seq_along(a1_hat)){
    yy = c(yy, exp(ifelse(xx < g_hat[i], 
                          a1_hat[i]*(g_hat[i] - xx)^2+t_hat[i], 
                          a2_hat[i]*(xx - g_hat[i])^2+t_hat[i])))
}
idid = rep(seq_along(a1_hat), each = nn)
dfdf = data.frame(y = yy, x = rep(xx, length(a1_hat)), Cell_ID_num = idid) %>% left_join(., (sholl_data %>% mutate(Cell_ID_num = as.numeric(as.factor(Cell_ID))) %>% select(Cell_ID_num, Animal_ID)), by = "Cell_ID_num")


# Stacked Cell Level Curves
cell_level_overlap_plot = ggplot(dfdf, aes(x = x, 
                                           y = y, 
                                           group = Cell_ID_num,
                                           color = as.factor(Animal_ID))) + 
    geom_line() + 
    ylab("Number of Intersections") + 
    xlab("Radius") + 
    ggtitle("Cell-Level Fitted Sholl Curves") + 
    guides(color = guide_legend(title = "Animal")) + 
    ylim(0, 50) +
    theme(legend.position = "bottom")
# cell_level_overlap_plot

# Cell Level Curves by Animal
dfdf2 = dfdf %>% group_by(Animal_ID) %>% mutate(Cell_ID_num =  Cell_ID_num - min(Cell_ID_num) + 1)
dfdf3 = data.frame(x = d$x,
                   y = d$y,
                   Cell_ID_num = as.numeric(as.factor(sholl_data$Cell_ID)),
                   Animal_ID = sholl_data$Animal_ID) %>% group_by(Animal_ID) %>% mutate(Cell_ID_num =  Cell_ID_num - min(Cell_ID_num) + 1)


cell_level_facet_plot = ggplot(dfdf2, aes(x = x, 
                                          y = y, 
                                          group = Cell_ID_num,
                                          color = as.factor(Cell_ID_num))) + 
    geom_line() + 
    geom_point(data = dfdf3, 
               aes(x = x, 
                   y = y, 
                   group = Cell_ID_num,
                   color = as.factor(Cell_ID_num)),
               alpha = 0.5,
               size = 0.8) + 
    facet_wrap(~Animal_ID, ncol = 2) + 
    xlim(0, 70) + 
    guides(color = "none") + 
    ylab("Number of Intersections") + 
    xlab("Radius") + 
    ggtitle("Cell-Level Fitted Sholl Curves\nby Animal") + 
    ylim(0, 50)
# cell_level_facet_plot

# Image Level
a1_hat = apply(samp$a1.image, 1, mean)
a2_hat = apply(samp$a2.image, 1, mean)
g_hat = apply(samp$g.image, 1, mean)
t_hat = apply(samp$t.image, 1, mean)

nn = 1000
maxX = 100
xx = seq(0, maxX, length.out = nn)
yy = c()
for(i in seq_along(a1_hat)){
    yy = c(yy, exp(ifelse(xx < g_hat[i], 
                          a1_hat[i]*(g_hat[i] - xx)^2+t_hat[i], 
                          a2_hat[i]*(xx - g_hat[i])^2+t_hat[i])))
}
idid = rep(seq_along(a1_hat), each = nn)
dfdf = data.frame(y = yy, x = rep(xx, length(a1_hat)), Image_ID_num = idid) %>% left_join(., (sholl_data %>% mutate(Image_ID_num = as.numeric(as.factor(Image_ID))) %>% select(Image_ID_num, Animal_ID)), by = "Image_ID_num")

# Image Level Curves by Animal (Dots are cell level)
dfdf2 = dfdf %>% group_by(Animal_ID) %>% mutate(Image_ID_num =  Image_ID_num - min(Image_ID_num) + 1)
dfdf3 = data.frame(x = d$x,
                   y = d$y,
                   Image_ID_num = as.numeric(as.factor(sholl_data$Image_ID)),
                   Animal_ID = sholl_data$Animal_ID) %>% group_by(Animal_ID) %>% mutate(Image_ID_num =  Image_ID_num - min(Image_ID_num) + 1)


image_level_overlap_plot = ggplot(dfdf, aes(x = x, 
                                            y = y, 
                                            group = Image_ID_num,
                                            color = as.factor(Animal_ID))) + 
    geom_line()  + 
    ylab("Number of Intersections") + 
    xlab("Radius") + 
    ggtitle("Image-Level Fitted Sholl Curves") + 
    guides(color = guide_legend(title = "Animal", nrow = 2))  + 
    ylim(0, 50) +
    theme(legend.position = "bottom")

image_level_facet_plot = ggplot(dfdf2, aes(x = x, 
                                           y = y, 
                                           group = Image_ID_num,
                                           color = as.factor(Image_ID_num))) + 
    geom_line() + 
    geom_point(data = dfdf3, 
               aes(x = x, 
                   y = y, 
                   group = Image_ID_num,
                   color = as.factor(Image_ID_num)),
               alpha = 0.5,
               size = 0.8) + 
    facet_wrap(~Animal_ID, ncol = 2) + 
    xlim(0, 70) + 
    guides(color = "none") + 
    ylab("Number of Intersections") + 
    xlab("Radius") + 
    ggtitle("Image-Level Fitted Sholl Curves\nby Animal")  + 
    ylim(0, 50)
# image_level_facet_plot

# Image Level Curves by Animal (Dots are within image means)
dfdf2 = dfdf %>% group_by(Animal_ID) %>% mutate(Image_ID_num =  Image_ID_num - min(Image_ID_num) + 1)
dfdf3 = data.frame(x = d$x,
                   y = d$y,
                   Image_ID_num = as.numeric(as.factor(sholl_data$Image_ID)),
                   Animal_ID = sholl_data$Animal_ID) %>% 
    group_by(Animal_ID) %>% 
    mutate(Image_ID_num =  Image_ID_num - min(Image_ID_num) + 1) %>%
    ungroup %>%
    group_by(Image_ID_num, Animal_ID, x) %>%
    summarise(y = mean(y))

image_level_facet_plot_mean = ggplot(dfdf2, aes(x = x, 
                                                y = y, 
                                                group = Image_ID_num,
                                                color = as.factor(Image_ID_num))) + 
    geom_line() + 
    geom_point(data = dfdf3, 
               aes(x = x, 
                   y = y, 
                   group = Image_ID_num,
                   color = as.factor(Image_ID_num)),
               alpha = 0.5,
               size = 0.8) + 
    facet_wrap(~Animal_ID, ncol = 2) + 
    xlim(0, 70) + 
    guides(color = "none") + 
    ylab("Number of Intersections") + 
    xlab("Radius") + 
    ggtitle("Image-Level Fitted Sholl Curves\nby Animal")  + 
    ylim(0, 50)
# image_level_facet_plot_mean

# Animal Level
a1_hat = apply(samp$a1.animal, 1, mean)
a2_hat = apply(samp$a2.animal, 1, mean)
g_hat = apply(samp$g.animal, 1, mean)
t_hat = apply(samp$t.animal, 1, mean)

nn = 1000
maxX = 100
xx = seq(0, maxX, length.out = nn)
yy = c()
for(i in seq_along(a1_hat)){
    yy = c(yy, exp(ifelse(xx < g_hat[i], 
                          a1_hat[i]*(g_hat[i] - xx)^2+t_hat[i], 
                          a2_hat[i]*(xx - g_hat[i])^2+t_hat[i])))
}
idid = rep(seq_along(a1_hat), each = nn)
# dfdf = data.frame(y = yy, x = rep(xx, length(a1_hat)), id = idid)

dfdf = data.frame(y = yy, x = rep(xx, length(a1_hat)), Animal_ID_num = idid) %>% left_join(., (sholl_data %>% mutate(Animal_ID_num = as.numeric(as.factor(Animal_ID))) %>% select(Animal_ID_num, Animal_ID)), by = "Animal_ID_num")

# Animal Image Level Curves
animal_level_overlap_plot = ggplot(dfdf, aes(x = x, 
                                             y = y, 
                                             color = Animal_ID)) + 
    geom_line() + 
    ylab("Number of Intersections") + 
    xlab("Radius") + 
    ggtitle("Animal-Level Fitted Sholl Curves") + 
    ylim(0, 50) +
    theme(legend.position = "bottom")
# animal_level_overlap_plot

# Animal Image Level Curves
animal_level_facet_plot = ggplot(dfdf, aes(x = x, 
                                           y = y, 
                                           color = Animal_ID)) + 
    geom_line() + 
    geom_point(data = data.frame(x = d$x, 
                                 y = d$y, 
                                 Animal_ID = sholl_data$Animal_ID), aes(x = x, 
                                                                        y = y, 
                                                                        color = Animal_ID), 
               alpha = 0.5,
               size = 0.8) + 
    facet_wrap(~Animal_ID, ncol = 2) + 
    xlim(0, 70) + 
    guides(color = "none")+ 
    ylab("Number of Intersections") + 
    xlab("Radius") + 
    ggtitle("Animal Level Fitted Sholl Curves\nby Animal")  + 
    ylim(0, 50)
# animal_level_facet_plot

#readjust font size
font_size = 14


### Cowplot for paper
cell_level_overlap_plot = cell_level_overlap_plot + 
    theme(text = element_text(size = font_size))

image_level_overlap_plot = image_level_overlap_plot + 
    theme(text = element_text(size = font_size))

animal_level_overlap_plot = animal_level_overlap_plot + 
    theme(text = element_text(size = font_size))

cell_level_facet_plot = cell_level_facet_plot + 
    theme(text = element_text(size = font_size))

image_level_facet_plot = image_level_facet_plot + 
    theme(text = element_text(size = font_size))

animal_level_facet_plot = animal_level_facet_plot + 
    theme(text = element_text(size = font_size))

# Using the cowplot package
legend <- cowplot::get_legend(image_level_overlap_plot)

top_row = plot_grid(
    cell_level_facet_plot,
    image_level_facet_plot + ylab(NULL),
    animal_level_facet_plot + ylab(NULL),
    nrow = 1,
    labels = c("A", "B", "C")
)

bottom_row = plot_grid(
    cell_level_overlap_plot + theme(legend.position="none"),
    image_level_overlap_plot + theme(legend.position="none") + ylab(NULL),
    animal_level_overlap_plot + theme(legend.position="none") + ylab(NULL),
    nrow = 1,
    labels = c("D", "E", "F")
)


p = plot_grid(top_row,
              bottom_row,
              legend,
              ncol = 1,
              rel_heights = c(5, 1, 0.2)
)

save_dir = "./Out/Plots/ungrouped_animal"
cowplot::save_plot(paste0(save_dir, '/', 'grouped_plot.png'), p, nrow = 3, ncol = 3, base_height = 7, base_width = 4.5)






