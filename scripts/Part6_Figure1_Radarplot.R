### Part 6: Figure 1: Radarplots: species' interaction generalism ####


## Packages ####
library (tidyverse)
library (ggradar)
library (ggpubr)
library (scales)
library (rstatix)
library (flextable)
library (RVAideMemoire)


## Radar plots per experiment (not shown in publication) ####

# plot_radar_E1<- 
#   All_metrics_E1 %>%  
#   dplyr::select (-PlantFamily, -Exp) %>% 
#   ungroup() %>% 
#   mutate(across(!c(PlantSpeciesfull), rescale)) %>% 
#   ggradar(legend.title = "Plant species", 
#           plot.title = "Experiment 1")
# 
# plot_radar_E2 <- 
#   All_metrics_E2 %>%  ungroup () %>% 
#   dplyr::select (-PlantFamily, -Exp) %>% 
#   mutate(across(!c(PlantSpeciesfull), rescale)) %>% # all metrics get scaled to obtain their relative values to each other, 
#   #i.e. from most specialist to most generalist
#   ggradar(legend.title = "Plant species", plot.title = "Experiment 2", 
#           background.circle.colour = "white",
#           background.circle.transparency = 0, 
#           axis.line.colour = "grey30", 
#           gridline.min.colour = "grey30",
#           gridline.mid.colour = "grey30",
#           gridline.max.colour = "grey30") 
# 
# 
# ggarrange (plot_radar_E1, plot_radar_E2, common.legend = T, legend = "bottom", label = "AUTO")


# Empty radarplot to include in Figure 1 ####
extra <- tribble (~PlantSpeciesfull,~Richness, ~Shannon, ~MPD, ~y_MPD, ~uniqueASV, ~CU, ~`β(core)`, ~Unifrac, "leer", 0, 0, 0, 0, 0, 0, 0,0)
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
  select (PlantSpeciesfull, `Richness S`, `Shannon's H'`, MPD, `Phyl. γ-diversity`, ,`γ-diversity`, `β(CU)` ,`β(core)`, Unifrac, Exp ) %>% 
    group_split(PlantSpeciesfull) %>% 
  map (~dplyr::select (., -PlantSpeciesfull) %>% relocate (Exp)) %>% 
  map (~ggradar(.,          background.circle.colour = "white",
                background.circle.transparency = 0, 
                axis.line.colour = "grey30", 
                gridline.min.colour = "grey30",
                gridline.mid.colour = "grey30",
                gridline.max.colour = "grey30",
                values.radar = c(" ", " ", " "),
                axis.labels = c("", "","", " ", "", "", "", ""))) # remove labels

## Arrange the plots ##
ggarrange (NULL, radar_plots[[1]], radar_plots[[2]], radar_plots[[3]], radar_plots[[4]], 
           radar_plots[[5]], radar_plots[[6]], radar_plots[[7]], radar_plots[[8]], 
           labels = c("", "A. capillaris", "A. millefolium", "B. willdenowii", "C. intybus", 
                      "H. lanatus", "P. cita", "P. lanceolata", "S. arundinaceus"),
           font.label = list (face = "italic"),
           common.legend = T, nrow = 3, ncol=3)

ggsave ("figures/testFig1.svg", width = 18, height = 18)

# Note: Further formatting of Figure 1 was applied manually using inkscape.


# # sorted for generalism according to table below (not shown) 
# ggarrange (radar_plots[[5]], radar_plots[[8]], radar_plots[[7]], radar_plots[[1]], 
#            radar_plots[[3]], radar_plots[[6]], radar_plots[[2]], radar_plots[[4]], 
#            labels = c("H. lanatus", "S. arundinaceus", "P. lanceolata","A. capillaris",
#                       "B. willdenowii","P. cita", "A. millefolium",  "C. intybus" 
#            ),
#            common.legend = T, nrow = 2, ncol=4)



## Procrustes analysis of the absolute values of the diversity metrics  #####

## Packages ####
library(vegan)


# Permutational Procrustes test #####
procr_test<- 
  All_metrics_E1 %>%  select (-Exp, -PlantFamily) %>% 
  pivot_longer(!PlantSpeciesfull, names_to = "metric", values_to = "Mvalue") %>%  
  left_join(All_metrics_E2 %>%  ungroup () %>% select (-Exp, -PlantFamily) %>% 
              pivot_longer(!PlantSpeciesfull, names_to = "metric", values_to = "valuesE2")) %>%
  split(~PlantSpeciesfull) %>% 
  map (~protest(X= .$Mvalue, Y = .$valuesE2, scale = F, scores = "sites", 
                permutations = 9999, symmetric = T))


# Randomisation, i.e. (permutational) test estimate the significance of an observed statistic
# relative to a large number of values for the same statistic that are generataed 
# by permuting the original data #
# 
# Save a table with the results #
map_dfr(procr_test, ~unlist (c(.$ss, .$svd$d, .$signif)) ) %>% 
  add_column ("procr" = c("M2", "correl", "sig")) %>%  
  mutate (across(where (is.numeric), ~round(.x, digits = 5))) %>% 
  write.csv("results/2024SymTRUE_Procrustes_test.csv")



Test: set.seed(123)  # for reproducibility


# Function to calculate mean Euclidean distance
calc_mean_euclidean <- function(x, y) {
  mean(sqrt((x - y)^2))
}

# Observed mean Euclidean distance
observed_distance <-   All_metrics_E1 %>%  select (-Exp, -PlantFamily) %>% 
  pivot_longer(!PlantSpeciesfull, names_to = "metric", values_to = "Mvalue") %>%  
  left_join(All_metrics_E2 %>%  ungroup () %>% select (-Exp, -PlantFamily) %>% 
              pivot_longer(!PlantSpeciesfull, names_to = "metric", values_to = "valuesE2")) %>%
  split(~PlantSpeciesfull) %>% 
  map (~calc_mean_euclidean( .$Mvalue,  .$valuesE2))

calc_mean_euclidean <- function(df) {
  mean(sqrt((df$Mvalue - df$valuesE2)^2))
}

# Function to perform randomization test on one data frame
randomization_test <- function(df, n_permutations = 1000) {
  observed_distance <- calc_mean_euclidean(df)
  
  # Perform randomization (shuffle column2 and recalculate distances)
  permuted_distances <- replicate(n_permutations, {
    shuffled_column2 <- sample(df$valuesE2)
    calc_mean_euclidean(data.frame(Mvalue = df$Mvalue , valuesE2 = shuffled_column2))
  })
  
  # Calculate p-value (one-sided test)
  p_value <- mean(permuted_distances <= observed_distance)
  
  list(observed_distance = observed_distance, p_value = p_value)
}

data_list <- All_metrics_E1 %>%  select (-Exp, -PlantFamily) %>% 
  pivot_longer(!PlantSpeciesfull, names_to = "metric", values_to = "Mvalue") %>%  
  left_join(All_metrics_E2 %>%  ungroup () %>% select (-Exp, -PlantFamily) %>% 
              pivot_longer(!PlantSpeciesfull, names_to = "metric", values_to = "valuesE2")) %>%
  split(~PlantSpeciesfull)

# Apply randomization test to each object in the list
results <- lapply(data_list, randomization_test)

# View results
results


#The p_value will tell you the probability that the observed mean Euclidean distance is smaller than or equal to the distances from the random shuffles. A small p-value (e.g., < 0.05) indicates that the distances between the columns are significantly 
#smaller than expected by chance, meaning they are more similar than random.
