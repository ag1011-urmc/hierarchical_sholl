library(tidyverse)
library(latex2exp)
library(HDInterval)

########################################################################
### Loading Model Fits
########################################################################

samp <- readRDS("/gpfs/fs2/scratch/evonkaen/Microglia/Sholl Analysis/Simulated_Fits/CrushKO/crushcontrol_homoskedastic_updatedData.RDS")

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


samp = fits2df(samp)


########################################################################
### Loading data
########################################################################


dat_dir <- '../../Data/Alexis_Data/Spreadsheet_retinaRGC'

dat_dir <- './Data/Alexis_Data/Spreadsheet_retinaRGC'

#get directories with cell level data
sample_dirs <- list.files(dat_dir)[!grepl('\\.', list.files(dat_dir))]

retina_dat <- c()
for(i in seq_along(sample_dirs)){
    tmp_dir <- paste(dat_dir, sample_dirs[i], sep = '/')
    sholl_files <- list.files(tmp_dir)[grepl('.csv', list.files(tmp_dir))]
    sholl_file_names <- strsplit(sholl_files, '\\.') %>% sapply(., FUN = function(list_elem) list_elem[[1]])
    for(j in seq_along(sholl_files)){
        tmp_dat <- read.csv(paste(tmp_dir, sholl_files[j], sep = '/')) %>% mutate(Subject_ID = sample_dirs[i], Cell_ID = sholl_file_names[j])
        tmp_dat2 <- data.frame(Radius = seq(tmp_dat$Radius[length(tmp_dat$Radius)], 100.49, by = 2), Inters. = 0, Radius..norm.Area = 0, Inters..Area = 0, log.Radius = 0, log.Inters..Area = 0, Subject_ID = tmp_dat$Subject_ID[1], Cell_ID = tmp_dat$Cell_ID[1])
        retina_dat <- rbind(retina_dat, tmp_dat)
    }
}

retina_dat %>% select(Subject_ID, Animal_ID) %>% distinct

#remove cells which had poor Sholl data
# retina_dat = retina_dat %>% filter(!(Cell_ID %in% c("MAX_C1-59875_left_mgROI_9of11_Sholl-Profiles",
#                                                     "MAX_C1-59876_right_mgROI_1of13_Sholl-Profiles", 
#                                                     "MAX_C1-60043_left_mgROI_2of5_Sholl-Profiles")))


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
    mutate(which_eye = ifelse(grepl("left", Subject_ID), "left" , "right")) %>% #create left\right identifier
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


########################################################################
### Create Plots
########################################################################

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
ggsave("Plots/CrushKO/cell_level_facet_NewData.png", cell_level_facet_plot, width = 10, height = 7)




### Sanity checks on fitted curves
tmp_ggdat = ggdat %>% left_join(., y = retina_dat %>% 
                                    select(Animal_ID, Cell_ID_factor, Cell_ID) %>% 
                                    distinct %>% 
                                    mutate(id = Cell_ID_factor) %>% 
                                    select(-Cell_ID_factor), 
                                by = c("Animal_ID", "id"))

for(an_id in unique(tmp_ggdat %>% pull(Animal_ID))){
    animal_level_facet_plot_by_group = ggplot(tmp_ggdat %>% filter(Animal_ID == an_id), aes(x = x, y = y, color = as.factor(id), group = as.factor(id)))+
        geom_line(size = 0.25, alpha = 0.8) + 
        geom_point(data = retina_dat %>% filter(Animal_ID == an_id) %>% mutate(`Genotype/Group` = paste(genotype, group, sep = "/"), id = Cell_ID_factor), aes(x = Radius, y = `Inters.`, group = Cell_ID, color = as.factor(id)))+ 
        geom_vline(xintercept = 50, color = "red", linetype = "dotted", size = 0.5) + 
        facet_grid(Cell_ID ~ `Genotype/Group`) +
        theme(legend.position = "bottom") + 
        ylab("Number of Intersections") + 
        xlab("Radius") + 
        guides(color = guide_legend(title = "Genotype / Group", override.aes = list(size = 1))) +
        # scale_color_manual(values = c("#B8DE29FF", "#2D708EFF", "#29AF7FFF", "#482677FF")) + 
        # ggtitle("Cell Level Fitted Curves by Animal") + 
        theme(text = element_text(size = 14)) + 
        theme(strip.text.y = element_text(size = 2.8)) + 
        theme(legend.position = "none") + 
        theme(plot.title = element_text(size = 14)) + 
        ggtitle(paste0("Cell-level fitted curves for Animal ", an_id, "\nSuperimposed over cell-level Sholl curves")) + 
        labs(caption = "Red dotted line at x = 50\nAppears cell-level Sholl data gets cut off at x = 50")
    ggsave(paste0("Plots/CrushKO/cell_level_curves_faceted/cell_level_",an_id, ".png"), animal_level_facet_plot_by_group, width = 5, height = 40)
}



retina_dat %>% filter(Animal_ID == an_id) %>% filter(Cell_ID_factor == 147)



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

ggsave(filename = "Plots/CrushKO/effects_NewData.png", 
       plot = posterior_interval_plot, 
       width = 8.5, 
       height = 5.2)



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
























