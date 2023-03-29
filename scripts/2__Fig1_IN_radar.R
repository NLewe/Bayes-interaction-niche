### Fig: Radar of species ####




library (ggradar)

library (scales)



## radar plots per experiment ####

plot_radar_E1<- 
  All_metrics_E1_8Sp %>%  
  select (-PlantFamily, -Exp) %>% 
  ungroup() %>% 
  mutate(across(!PlantSpeciesfull, rescale)) %>% 
  ggradar(legend.title = "Plant species", 
          plot.title = "Experiment 1")

plot_radar_E2 <- 
All_metrics_E2 %>%  ungroup () %>% 
  select (-PlantFamily, -Exp) %>% 
  mutate (across (!PlantSpeciesfull, rescale)) %>% 
  ggradar(legend.title = "Plant species", plot.title = "Experiment 2") 


ggarrange (plot_radar_E1, plot_radar_E2, common.legend = T, legend = "bottom", label = "AUTO")

# radar plot for each plant species ######
radar_plots <- 
  All_metrics_E1_8Sp %>%  
  ungroup() %>% 
  mutate(across(!c(PlantSpeciesfull,PlantFamily,Exp), rescale)) %>% 
      bind_rows (All_metrics_E2 %>% 
                   ungroup() %>% 
                   mutate(across(!c(PlantSpeciesfull,PlantFamily,Exp), rescale))) %>% 
                   select (-PlantFamily) %>% 
  group_split(PlantSpeciesfull) %>% 
  map (~select (., -PlantSpeciesfull) %>% relocate (Exp)) %>% 
  map (~ggradar(.,)) 


ggarrange (radar_plots[[1]], radar_plots[[2]], radar_plots[[3]], radar_plots[[4]], 
           radar_plots[[5]], radar_plots[[6]], radar_plots[[7]], radar_plots[[8]], 
           labels = c("A. capillaris", "A. millefolium", "B. willdenowii", "C. intybus", 
                      "H. lanatus", "P. cita", "P. lanceolata", "S. arundinaceus"),
           common.legend = T, nrow = 2, ncol=4)





  
  






