library(dplyr)
library(tibble)
options(dplyr.summarise.inform = FALSE)


# dat_dir = "MDND/"
dat_dir = "MDND/"

#initialize overall summary storage
pow_fpr_tbl_t = c()
par_est_tbl_t = c()
pow_fpr_tbl_g = c()
par_est_tbl_g = c()
for(j in seq_along(list.files(paste0("./Simulated_Data/", dat_dir)))){
  scen = j
  
  dat_files = list.files(paste0(paste0("./Simulated_Data/", dat_dir), "/scenario_", scen, "/"))[!grepl("plots", list.files(paste0(paste0("./Simulated_Data/", dat_dir), "/scenario_", scen, "/")))]
  dat_seeds = sapply(strsplit(dat_files, split = "\\."), function(x) x[1])
  
  #initialize summary storage for scenario
  dat_inds = seq_along(dat_files)
  
  p_val_mat_t = matrix(nrow = length(dat_inds), ncol = 3)
  coef_mat_t = matrix(nrow = length(dat_inds), ncol = 4)
  p_val_mat_g = matrix(nrow = length(dat_inds), ncol = 3)
  coef_mat_g = matrix(nrow = length(dat_inds), ncol = 4)
  
  for(i in dat_inds){
    if((i %% 10) == 0) cat("Scenario ", scen, ": Fitting dataset ", i, " out of ", length(dat_files), "\r")
    
    dat = readRDS(paste0(paste0("./Simulated_Data/", dat_dir), "/scenario_", scen, "/", dat_files[i]))
    
    #coerce dat$d into data.frame and then aggregate cell level data on animal
    aggregate_dat = data.frame(y = dat$d$y, x = dat$d$x, Cell_ID = dat$d$Cell_ID) %>% 
      left_join(x = ., y = data.frame(Cell_ID = seq_along(dat$d$Animal_ID), Animal_ID = dat$d$Animal_ID, Cond_Side_ID = dat$d$Cond_Side_ID_forCell), by = "Cell_ID") %>%
      left_join(x = ., y = data.frame(Cond_Side_ID = 1:4, Side_ID = dat$d$Side_ID, Cond_ID = dat$d$Cond_ID), by = "Cond_Side_ID")
    
    aggregate_dat = aggregate_dat %>% group_by(x, Animal_ID, Cond_Side_ID, Side_ID, Cond_ID) %>% summarise(y = mean(y))
    
    #get max intercept for each animal level aggregate curve to feed into ANOVA
    max_int_dat = aggregate_dat %>% ungroup %>% group_by(Animal_ID, Cond_Side_ID, Side_ID, Cond_ID) %>% summarise(max_y = max(y), x = x[which.max(y)])
    
    #t parameter
    m = aov(max_y ~ Side_ID * Cond_ID, data = max_int_dat)
    p_val_mat_t[i,] = summary(m)[[1]]$`Pr(>F)`[-4]
    coef_mat_t[i,] = m$coefficients
    
    #g parameter
    m = aov(x ~ Side_ID * Cond_ID, data = max_int_dat)
    p_val_mat_g[i,] = summary(m)[[1]]$`Pr(>F)`[-4]
    coef_mat_g[i,] = m$coefficients
  }
  
  pow_fpr_tbl_t = rbind(pow_fpr_tbl_t, colMeans(p_val_mat_t < 0.05))
  par_est_tbl_t = rbind(par_est_tbl_t, colMeans(coef_mat_t))
  pow_fpr_tbl_g = rbind(pow_fpr_tbl_g, colMeans(p_val_mat_g < 0.05))
  par_est_tbl_g = rbind(par_est_tbl_g, colMeans(coef_mat_g))
}

#make labels appropriate
rownames(pow_fpr_tbl_t) = paste0("Scenario ", seq_along(list.files(paste0("./Simulated_Data/", dat_dir))))
rownames(par_est_tbl_t) = paste0("Scenario ", seq_along(list.files(paste0("./Simulated_Data/", dat_dir))))
colnames(pow_fpr_tbl_t) = c("Side", "Cond", "Interaction")

rownames(pow_fpr_tbl_g) = paste0("Scenario ", seq_along(list.files(paste0("./Simulated_Data/", dat_dir))))
rownames(par_est_tbl_g) = paste0("Scenario ", seq_along(list.files(paste0("./Simulated_Data/", dat_dir))))
colnames(pow_fpr_tbl_g) = c("Side", "Cond", "Interaction")

#t power/fpr table
pow_fpr_tbl_t

#g fpr table
pow_fpr_tbl_g


