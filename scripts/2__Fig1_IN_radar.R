### Fig: Radar of species ####



library (tidyverse)
library (ggradar)
library (ggpubr)
library (scales)
library (rstatix)
library (flextable)



## radar plots per experiment ####

plot_radar_E1<- 
  All_metrics_E1_8Sp %>%  
  select (-PlantFamily, -Exp) %>% 
  ungroup() %>% 
  mutate(across(!c(PlantSpeciesfull,`β(core)`, CU), rescale)) %>% ### achtung,`β(core)` aus de Berechnung genommen da es schon auf o bis 1 scala ist
  mutate (`β(core)` = `β(core)`/100) %>% 
  ggradar(legend.title = "Plant species", 
          plot.title = "Experiment 1")

plot_radar_E2 <- 
All_metrics_E2 %>%  ungroup () %>% 
  select (-PlantFamily, -Exp) %>% 
  mutate(across(!c(PlantSpeciesfull,`β(core)`, CU), rescale)) %>% ### achtung,`β(core)` aus de Berechnung genommen da es schon auf o bis 1 scala ist
  mutate (`β(core)` = `β(core)`/100) %>% 
  ggradar(legend.title = "Plant species", plot.title = "Experiment 2") 


ggarrange (plot_radar_E1, plot_radar_E2, common.legend = T, legend = "bottom", label = "AUTO")

# radar plot for each plant species ######
radar_plots <- 
  All_metrics_E1_8Sp %>%  
  ungroup() %>% 
  mutate(across(!c(PlantSpeciesfull,PlantFamily,Exp, `β(core)`, CU), rescale)) %>% 
        bind_rows (All_metrics_E2 %>% 
                   ungroup() %>% 
                   mutate(across(!c(PlantSpeciesfull,PlantFamily,Exp,`β(core)`, CU ), rescale))) %>% 
                   select (-PlantFamily) %>% 
  mutate (`β(core)` = `β(core)`/100) %>% 
   # select (!uniqueASV) %>% 
  group_split(PlantSpeciesfull) %>% 
  map (~select (., -PlantSpeciesfull) %>% relocate (Exp)) %>% 
  map (~ggradar(.,)) 


ggarrange (radar_plots[[1]], radar_plots[[2]], radar_plots[[3]], radar_plots[[4]], 
           radar_plots[[5]], radar_plots[[6]], radar_plots[[7]], radar_plots[[8]], 
           labels = c("A. capillaris", "A. millefolium", "B. willdenowii", "C. intybus", 
                      "H. lanatus", "P. cita", "P. lanceolata", "S. arundinaceus"),
           common.legend = T, nrow = 2, ncol=4)

# sorted for generalism according to table below
ggarrange (radar_plots[[5]], radar_plots[[8]], radar_plots[[7]], radar_plots[[1]], 
           radar_plots[[3]], radar_plots[[6]], radar_plots[[2]], radar_plots[[4]], 
           labels = c("H. lanatus", "S. arundinaceus", "P. lanceolata","A. capillaris",
                      "B. willdenowii","P. cita", "A. millefolium",  "C. intybus" 
                        ),
           common.legend = T, nrow = 2, ncol=4)

# Metric of generalism , relative #####
RelGen_E1_E2 <- 
All_metrics_E1_8Sp %>%  
  ungroup() %>% 
  mutate(across(!c(PlantSpeciesfull,PlantFamily,Exp, `β(core)`, CU), rescale)) %>% 
  bind_rows (All_metrics_E2 %>% 
               ungroup() %>% 
               mutate(across(!c(PlantSpeciesfull,PlantFamily,Exp,`β(core)`, CU), rescale))) %>% 
  select (-PlantFamily)  %>% 
  mutate (RelGenSpec = (Shannon + Richness + uniqueASV + CU + `β(core)`/100 + PD + MPD) /7 )

RelGen_E1_E2 %>%  reorder () %>% 
  mutate (round (across (), 2))
                          
(flextable ()

RelGen_E1_E2_sample <-
All_Metrics_E2_sample %>% 
  ungroup() %>% 
  mutate(across(!c(PlantSpeciesfull,PlantFamily,PlaSpe, Exp,`β(core)`, CU, sampleID), rescale)) %>% 
  bind_rows (All_metrics_E1_samples %>%    
               ungroup() %>% 
               mutate(across(!c(PlantSpeciesfull,PlantFamily,PlaSpe,`β(core)`, CU, Exp, sampleID), rescale))) %>% 
  mutate (RelGenSpec = (Shannon + Richness + uniqueASV + CU + `β(core)`/100 + PD + MPD) /7 ) %>% 
  relocate (where (is.numeric)) 

# experiment 2 sorting for highest generalist
RelGen_E2_mean <- 
  RelGen_E1_E2%>% 
  group_by  (Exp, PlantSpeciesfull) %>%  summarize (meanGen = mean (RelGenSpec)) %>% 
  arrange (-meanGen) %>%  
  filter (Exp == "E2")


RelGen_E1_E2_mean <- 
  RelGen_E1_E2%>% 
  group_by  (Exp, PlantSpeciesfull) %>%  summarize (meanGen = mean (RelGenSpec)) %>% 
  arrange (-meanGen) 

# Plot 
# plot using the values from radarplot ##
RelGen_E1_E2 %>% 
  left_join (RelGen_E2_mean %>% ungroup () %>%   select (!Exp)) %>%
  ggplot (aes (x = PlantSpeciesfull, y = RelGenSpec, fill = Exp)) + 
  geom_col(position = position_dodge ())


RelGen_E1_E2 %>% 
  left_join (RelGen_E2_mean %>% ungroup () %>%   select (!Exp)) %>%
  ggplot (aes (x = reorder (PlantSpeciesfull, -meanGen), y = RelGenSpec, fill = Exp)) + 
  geom_col() +
  facet_grid(~Exp)



## plot using the sample values ###
 RelGen_E1_E2_sample %>% ggplot (aes (x = PlantSpeciesfull, y = RelGenSpec, fill = Exp)) + 
  geom_boxplot (position = position_dodge ())



 ### I think sorting alphabetically looks better!!!
  
 RelGen_E1_E2_sample %>%  group_by ( PlantSpeciesfull ) %>%  t_test (RelGenSpec ~ Exp)


 RelGen_E1_E2_sample %>%  
   group_by ( PlantSpeciesfull ) %>%  
   t_test (RelGenSpec ~ Exp) %>%  
   flextable() %>%  save_as_image("figures/ttest_boxplot_InNi.png")



