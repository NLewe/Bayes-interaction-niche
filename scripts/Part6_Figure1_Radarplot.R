### Part 6: Figure 1: Radar of species' interaction generalism ####


## packages ####
library (tidyverse)
library (ggradar)
library (ggpubr)
library (scales)
library (rstatix)
library (flextable)
library (RVAideMemoire)


## radar plots per experiment (not shown) ####

plot_radar_E1<- 
  All_metrics_E1 %>%  
  dplyr::select (-PlantFamily, -Exp) %>% 
  ungroup() %>% 
  mutate(across(!c(PlantSpeciesfull), rescale)) %>% 
  ggradar(legend.title = "Plant species", 
          plot.title = "Experiment 1")

plot_radar_E2 <- 
  All_metrics_E2 %>%  ungroup () %>% 
  dplyr::select (-PlantFamily, -Exp) %>% 
  mutate(across(!c(PlantSpeciesfull), rescale)) %>% # all metrics get scaled to obtain their relative values to each other, 
  #i.e. from most specialist to most generalist
  ggradar(legend.title = "Plant species", plot.title = "Experiment 2", 
          background.circle.colour = "white",
          background.circle.transparency = 0, 
          axis.line.colour = "grey30", 
          gridline.min.colour = "grey30",
          gridline.mid.colour = "grey30",
          gridline.max.colour = "grey30") 





ggarrange (plot_radar_E1, plot_radar_E2, common.legend = T, legend = "bottom", label = "AUTO")

# Prepare an empty radarplot for the figure ##
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

# Radar plot for each plant species ######
radar_plots <- 
  All_metrics_E1 %>%  
  ungroup() %>% 
  mutate(across(!c(PlantSpeciesfull,PlantFamily,Exp), rescale)) %>% 
  bind_rows (All_metrics_E2 %>% 
               ungroup() %>% 
               mutate(across(!c(PlantSpeciesfull,PlantFamily,Exp ), rescale))) %>% 
  dplyr:: select (!PlantFamily) %>% 
  select (PlantSpeciesfull, `Richness S`, `Shannon's H'`, MPD, `Phyl. γ-diversity`, ,`γ-diversity`, `β(CU)` ,`β(core)`, Exp ) %>% 
    group_split(PlantSpeciesfull) %>% 
  map (~dplyr::select (., -PlantSpeciesfull) %>% relocate (Exp)) %>% 
  map (~ggradar(.,          background.circle.colour = "white",
                background.circle.transparency = 0, 
                axis.line.colour = "grey30", 
                gridline.min.colour = "grey30",
                gridline.mid.colour = "grey30",
                gridline.max.colour = "grey30",
                values.radar = c(" ", " ", " "),
                axis.labels = c("", "","", " ", "", "", "")))

## Arrange the plots ##
ggarrange (NULL, radar_plots[[1]], radar_plots[[2]], radar_plots[[3]], radar_plots[[4]], 
           radar_plots[[5]], radar_plots[[6]], radar_plots[[7]], radar_plots[[8]], 
           labels = c("", "A. capillaris", "A. millefolium", "B. willdenowii", "C. intybus", 
                      "H. lanatus", "P. cita", "P. lanceolata", "S. arundinaceus"),
           font.label = list (face = "italic"),
           common.legend = T, nrow = 3, ncol=3)

# Further labelling for the plot was added manually using inkscape


# sorted for generalism according to table below (not shown) 
ggarrange (radar_plots[[5]], radar_plots[[8]], radar_plots[[7]], radar_plots[[1]], 
           radar_plots[[3]], radar_plots[[6]], radar_plots[[2]], radar_plots[[4]], 
           labels = c("H. lanatus", "S. arundinaceus", "P. lanceolata","A. capillaris",
                      "B. willdenowii","P. cita", "A. millefolium",  "C. intybus" 
           ),
           common.legend = T, nrow = 2, ncol=4)




library (FactoMineR)
library(factoextra)
library(ggrepel)
library (ggpubr)
library(vegan)



## Procrustes analysis of the absolute values of the diversity metrics  #####
procr_test<- 
  All_metrics_E1 %>%  select (-Exp, -PlantFamily) %>% 
  pivot_longer(!PlantSpeciesfull, names_to = "metric", values_to = "Mvalue") %>%  
  left_join(All_metrics_E2 %>%  ungroup () %>% select (-Exp, -PlantFamily) %>% 
              #mutate_at (vars (-PlantSpeciesfull), rescale) %>% 
              pivot_longer(!PlantSpeciesfull, names_to = "metric", values_to = "valuesE2")) %>%
  split(~PlantSpeciesfull) %>% 
  map (~protest(X= .$Mvalue, Y = .$valuesE2, scale = F, scores = "sites", permutations = 9999, symmetric = F))


# randomisation (permutational) test estimate the significance of an observed statistic
# relative to a large number of values for the same statistic that are generataed 
# by permuting the original data #
# get results in table #
map_dfr(procr_test, ~unlist (c(.$ss, .$svd$d, .$signif)) ) %>% 
  add_column ("procr" = c("M2", "correl", "sig")) %>%  
  mutate (across(where (is.numeric), ~round(.x, digits = 5))) %>% 
  write.csv("results/Procrustes_test.csv")


