### Fig: Radar of species ####



library (tidyverse)
library (ggradar)
library (ggpubr)
library (scales)
library (rstatix)
library (flextable)
library (RVAideMemoire)


## radar plots per experiment ####

plot_radar_E1<- 
  All_metrics_E1_8Sp %>%  
  dplyr::select (-PlantFamily, -Exp) %>% 
  ungroup() %>% 
  mutate(across(!c(PlantSpeciesfull), rescale)) %>% ### achtung,`β(core)` kann man aus de Berechnung genommen da es schon auf o bis 1 scala ist
  #mutate (`β(core)` = `β(core)`/100) %>% 
  ggradar(legend.title = "Plant species", 
          plot.title = "Experiment 1")

plot_radar_E2 <- 
  All_metrics_E2 %>%  ungroup () %>% 
  dplyr::select (-PlantFamily, -Exp) %>% 
  mutate(across(!c(PlantSpeciesfull), rescale)) %>% ### achtung,`β(core)` aus de Berechnung genommen da es schon auf o bis 1 scala ist
  #mutate (`β(core)` = `β(core)`/100) %>% 
  ggradar(legend.title = "Plant species", plot.title = "Experiment 2", 
          background.circle.colour = "white",
          background.circle.transparency = 0, 
          axis.line.colour = "grey30", 
          gridline.min.colour = "grey30",
          gridline.mid.colour = "grey30",
          gridline.max.colour = "grey30")





ggarrange (plot_radar_E1, plot_radar_E2, common.legend = T, legend = "bottom", label = "AUTO")

extra <- tribble (~PlantSpeciesfull,~Richness, ~Shannon, ~MPD, ~y_MPD, ~uniqueASV, ~CU, ~`β(core)`,  "leer", 0, 0, 0, 0, 0, 0, 0 )
radar_empty  <- ggradar(extra,          background.circle.colour = "white",
                        background.circle.transparency = 0, 
                        axis.line.colour = "grey30", 
                        gridline.min.colour = "grey30",
                        gridline.mid.colour = "grey30",
                        gridline.max.colour = "grey30", 
                        group.line.width=0, 
                        grid.label.size = 5,
                        group.point.size	=0)

# radar plot for each plant species ######
radar_plots <- 
  All_metrics_E1_8Sp %>%  
  ungroup() %>% 
  mutate(across(!c(PlantSpeciesfull,PlantFamily,Exp), rescale)) %>% 
        bind_rows (All_metrics_E2 %>% 
                   ungroup() %>% 
                   mutate(across(!c(PlantSpeciesfull,PlantFamily,Exp ), rescale))) %>% 
                   dplyr:: select (!PlantFamily) %>% 
  select (PlantSpeciesfull, `Richness S`, `Shannon's H'`, MPD, `Phyl. γ-diversity`, ,`γ-diversity`, `β(CU)` ,`β(core)`, Exp ) %>% 
  #mutate (`β(core)` = `β(core)`/100) %>% 
   # select (!uniqueASV) %>% 
  group_split(PlantSpeciesfull ) %>% 
  map (~dplyr::select (., -PlantSpeciesfull) %>% relocate (Exp)) %>% 
  map (~ggradar(.,          background.circle.colour = "white",
                background.circle.transparency = 0, 
                axis.line.colour = "grey30", 
                gridline.min.colour = "grey30",
                gridline.mid.colour = "grey30",
                gridline.max.colour = "grey30",
                values.radar = c(" ", " ", " "),
                axis.labels = c("", "","", " ", "", "", "")))


ggarrange (NULL, radar_plots[[1]], radar_plots[[2]], radar_plots[[3]], radar_plots[[4]], 
           radar_plots[[5]], radar_plots[[6]], radar_plots[[7]], radar_plots[[8]], 
           labels = c("", "A. capillaris", "A. millefolium", "B. willdenowii", "C. intybus", 
                      "H. lanatus", "P. cita", "P. lanceolata", "S. arundinaceus"),
           font.label = list (face = "italic"),
           common.legend = T, nrow = 3, ncol=3)



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
  mutate(across(!c(PlantSpeciesfull,PlantFamily,Exp), rescale)) %>% 
  bind_rows (All_metrics_E2 %>% 
               ungroup() %>% 
               mutate(across(!c(PlantSpeciesfull,PlantFamily,Exp), rescale))) %>% 
  dplyr::select (-PlantFamily)  %>% 
 # mutate (RelGenSpec = (Shannon + Richness + uniqueASV + CU + `β(core)` + PD + MPD) /7 )
  mutate (RelGenSpec = (Shannon + Richness + uniqueASV + CU + `β(core)` + y_MPD + MPD) /7 )

## table for plant species ###
RelGen_E1_E2 %>%  
  mutate (RelGenSpec = round (RelGenSpec, 3)) %>% 
  dplyr::select (PlantSpeciesfull, RelGenSpec, Exp) %>% 
  pivot_wider (values_from = RelGenSpec, names_from = Exp) %>% 
  flextable ()
                          

## Rel Gen for each sample ##
RelGen_E1_E2_sample <-
All_Metrics_E2_sample %>% 
  ungroup() %>% 
  mutate(across(!c(PlantSpeciesfull,PlantFamily,PlaSpe, Exp, sampleID), rescale)) %>% 
  bind_rows (All_metrics_E1_samples %>%    
               ungroup() %>% 
               mutate(across(!c(PlantSpeciesfull,PlantFamily,PlaSpe, Exp, sampleID), rescale))) %>% 
  #mutate (RelGenSpec = (Shannon + Richness + uniqueASV + CU + `β(core)` + PD + MPD) /7 ) %>% 
  mutate (RelGenSpec = (Shannon + Richness + uniqueASV + CU + `β(core)` + y_MPD + MPD) /7 )
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
  left_join (RelGen_E2_mean %>% ungroup () %>%   dplyr::select (!Exp)) %>%
  ggplot (aes (x = PlantSpeciesfull, y = RelGenSpec, fill = Exp)) + 
  geom_col(position = position_dodge ())


RelGen_E1_E2 %>% 
  left_join (RelGen_E2_mean %>% ungroup () %>%   dplyr::select (!Exp)) %>%
  ggplot (aes (x = reorder (PlantSpeciesfull, -meanGen), y = RelGenSpec, fill = Exp)) + 
  geom_col() +
  facet_grid(~Exp)



## plot using the sample values ###
 RelGen_E1_E2_sample %>% ggplot (aes (x = PlantSpeciesfull, y = RelGenSpec, fill = Exp)) + 
  geom_boxplot (position = position_dodge ())

 

# t test #####
 ### I think sorting alphabetically looks better!!!
  
 RelGen_E1_E2_sample %>%  
   group_by ( PlantSpeciesfull ) %>%  
 summarise(model = list(perm.t.test (RelGenSpec ~ Exp),
                                      data = cur_data(), nperm=999)) %>%
   mutate(res = map(model, broom::tidy)) %>%
   unnest(res)

 RelGen_E1_E2_sample %>%  
   group_by ( PlantSpeciesfull ) %>%  
   t_test (RelGenSpec ~ Exp) %>%  
   flextable() %>%  save_as_image("figures/ttest_boxplot_InNi.png")

 ## Permutational t test #####
Permut_Ttest_RelGenSpec <-
 RelGen_E1_E2_sample %>%  
    filter (!is.na (PD)) %>% 
group_split (PlantSpeciesfull) %>%  # split into each plant species
   map (~perm.t.test (.$RelGenSpec ~ .$Exp, nperm= 999)) 
   


