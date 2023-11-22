# part 7 - PCAs #####

## packages ####
library (tidyverse)
library (ggradar)
library (ggrepel)
library (ggpubr)
library (scales)
library (factoextra)
library (FactoMineR)

# PCA E2 per sample  (supplementary material) #####

# All_Metrics_E2_sample <- readRDS ("data/All_metrics_E2_sample.rds")
All_Metrics_E2_sample_df <- 
  All_Metrics_E2_sample %>% 
  relocate(where(is.numeric), .after = where(is.character)) %>%  data.frame(row.names = "sampleID")

## calculation of PCA ##
PCA_all_metrics_E2_sample <- PCA(All_Metrics_E2_sample_df, quali.sup = c(1:4),   scale.unit = T, graph = F)

#PCA plot Exp 2 -  for each sample #
PCA_arrows_metrics_sample<-
  PCA_all_metrics_E2_sample$var$coord %>%  
  as_tibble (rownames = "metric") %>% 
  dplyr::rename("D1end" = "Dim.1", "D2end"= "Dim.2", "D3end" = "Dim.3")

PCA_metric_data_sample   <- PCA_all_metrics_E2_sample$ind$coord  %>%  as_tibble(rownames = "sampleID") %>% 
  left_join(meta_M1 %>%  select (sampleID, PlantSpeciesfull, PlantFamily) ) 

PCA_metric_data_sample %>%  write_csv ("results/PCA_metric_PCs.csv")

saveRDS(PCA_metric_data_sample, "results/PCA_metric_data_sample.rds")
## for Appendix tbale 

PCA_eig_m_Dim1 <- round (PCA_all_metrics_E2_sample$eig[1,2],1)
PCA_eig_m_Dim1
PCA_eig_m_Dim2 <- round (PCA_all_metrics_E2_sample$eig[2,2],1)
PCA_eig_m_Dim2
PCA_eig_m_Dim3 <- round (PCA_all_metrics_E2_sample$eig[3,2],1)
PCA_eig_m_Dim3


#Plot PCA E2 per sample - dimensions 1 and 2 #
plotPca_dim1_2  <- 
  PCA_metric_data_sample %>%  
  mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
  ggplot (aes (x = Dim.1, y = Dim.2, color = PlantSpeciesfull)) + 
  geom_point (size =3) +
  geom_hline(yintercept = 0, lty = 2, color = "black", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "black", alpha = 0.9) + 
  #stat_ellipse(aes (x= Dim.1, y = Dim.2, color = group))  +   # thisis 95%confidence
  geom_segment(data = PCA_arrows_metrics_sample, aes (x=0, xend= D1end*4.5, y = 0, yend = D2end*4.5), 
               arrow = arrow(length = unit(0.3, "picas")), color = "blue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics_sample, aes (x  = D1end*4.5, y = D2end*4.5, label = metric), 
                    color = "blue", inherit.aes = F , force = 0.6) + 
  theme_classic() + 
  xlab(label = "PC1 (52.6 %)") +
  ylab ("PC3 (22.0 %)") + 
  guides (color= guide_legend( "Plant species")) +
  theme (legend.position = "bottom") 

# Plot contributions to PCs ##
plot1 <- fviz_contrib(PCA_all_metrics_E2_sample, choice = "var", axes = 1, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC1") +
  theme_classic()
plot2 <- fviz_contrib(PCA_all_metrics_E2_sample, choice = "var", axes = 2, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC2")
plot3 <- fviz_contrib(PCA_all_metrics_E2_sample, choice = "var", axes = 3, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC3")



# Plot PCA E2 per sample - dimensions 2 and 3 ##
plotPca_dim2_3  <- 
  PCA_metric_data_sample %>%  
  mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
  ggplot (aes (x = Dim.2, y = Dim.3, color = PlantSpeciesfull)) + 
  geom_point (size =3) +
  geom_hline(yintercept = 0, lty = 2, color = "black", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "black", alpha = 0.9) + 
  geom_segment(data = PCA_arrows_metrics_sample, aes (x=0, xend= D2end*4.5, y = 0, yend = D3end*4.5), 
               arrow = arrow(length = unit(0.3, "picas")), color = "blue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics_sample, aes (x  = D2end*4.5, y = D3end*4.5, label = metric), 
                    color = "blue", inherit.aes = F , force = 0.6) + 
  theme_classic() + 
  xlab(label = "PC2 (22.0 %)") +
  ylab ("PC3 (12.2 %)") + 
  guides (color= guide_legend( "Plant species")) +
  theme (legend.position = "bottom") 

# Full plot 
ggarrange (ggarrange (plotPca_dim1_2,plotPca_dim2_3, nrow = 1, ncol = 2, labels = c ("A", "B"), common.legend = T),
           ggarrange (plot1, plot2, plot3, nrow = 1, ncol= 3, labels= c( "C", "D", "E")), 
           nrow=2, ncol=1, heights  =  c(2, 1)) 






## additional plots PCA for supplementary ##

# PCA E1 per plant species ####
# 

PCA_all_metrics_E1 <- PCA(All_metrics_E1_df, quali.sup = c(2,9),   scale.unit = T, graph = F)

# get data for plot ###
PCA_arrows_metrics_E1<-
  PCA_all_metrics_E1$var$coord %>%  
  as_tibble (rownames = "metric") %>% 
  dplyr::rename("D1end" = "Dim.1", "D2end"= "Dim.2")

PCA_metric_data_E1   <- PCA_all_metrics_E1$ind$coord  %>%  as_tibble(rownames = "PlantSpeciesfull") %>% 
  left_join(metaM0 %>%  select (PlantSpeciesfull, PlantFamily) %>%  unique ())

PCA_eig_m_Dim1 <- round (PCA_all_metrics_E1$eig[1,2],1)
PCA_eig_m_Dim1
PCA_eig_m_Dim2 <- round (PCA_all_metrics_E1$eig[2,2],1)
PCA_eig_m_Dim2


#Plot  PCA for E1 - plant species (in supplementary) ###

plotPCA_E1  <- 
  PCA_metric_data_E1 %>%  
  mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
  ggplot (aes (x = Dim.1, y = Dim.2, color = PlantFamily)) + 
  geom_point (size =3) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  #stat_ellipse(aes (x= Dim.1, y = Dim.2, color = group))  +   # thisis 95%confidence
  geom_segment(data = PCA_arrows_metrics_E1, aes (x=0, xend= D1end*4.5, y = 0, yend = D2end*4.5), 
               arrow = arrow(length = unit(0.3, "picas")), color = "blue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics_E1, aes (x  = D1end*4.5, y = D2end*4.5, label = metric), 
                    color = "blue", inherit.aes = F , force = 0.6) + 
  geom_text_repel(data = PCA_metric_data_E1, aes (x = Dim.1, y = Dim.2, label = PlantSpeciesfull), fontface = "italic", inherit.aes = F) +
  theme_minimal() + 
  xlab(label = "PC1 (86.3 %)") +
  ylab ("PC2 (6.1 %)") + 
  guides (color= guide_legend( "Plant family")) +
  theme (legend.position = "bottom") +
  scale_color_manual(values=c("Asteraceae"= "#edc948" ,"Cyperaceae" ="#5A6351" ,
                              "Fabaceae" = "#f28e2b" , 
                              "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
                              "Poaceae" = "#7F9A65" )) 

# Plots for th most contributing variables ##
plot1 <- fviz_contrib(PCA_all_metrics_E1, choice = "var", axes = 1, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC1")
plot2 <- fviz_contrib(PCA_all_metrics_E1, choice = "var", axes = 2, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC2")
plot3 <- fviz_contrib(PCA_all_metrics_E1, choice = "var", axes = 3, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC3")

ggarrange (plotPCA_E1,
           ggarrange (plot1, plot2, plot3, nrow = 3, ncol= 1, labels= c("B", "C", "D")), 
           nrow=1, ncol=2, widths =  c(2, 1), labels= "A") 



# PCA E2 per plant species ####

PCA_all_metrics_E2 <- PCA(All_metrics_E2_df, quali.sup = c(8,9,10),   scale.unit = T, graph = F)

# data for the PCA plot ###
PCA_arrows_metrics_E2<-
  PCA_all_metrics_E2$var$coord %>%  
  as_tibble (rownames = "metric") %>% 
  dplyr::rename("D1end" = "Dim.1", "D2end"= "Dim.2")

PCA_metric_data_E2   <- PCA_all_metrics_E2$ind$coord  %>%  as_tibble(rownames = "PlantSpeciesfull") %>% 
  left_join(metaM0 %>%  select (PlantSpeciesfull, PlantFamily) %>%  unique ())


PCA_eig_m_Dim1 <- round (PCA_all_metrics_E2$eig[1,2],1)
PCA_eig_m_Dim1
PCA_eig_m_Dim2 <- round (PCA_all_metrics_E2$eig[2,2],1)
PCA_eig_m_Dim2

# PCA Plot E2 # 

plotPCA_E2  <- 
  PCA_metric_data_E2 %>%  
  mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
  ggplot (aes (x = Dim.1, y = Dim.2, color = PlantFamily)) + 
  geom_point (size =3) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_segment(data = PCA_arrows_metrics_E2, aes (x=0, xend= D1end*4.5, y = 0, yend = D2end*4.5), 
               arrow = arrow(length = unit(0.3, "picas")), color = "blue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics_E2, aes (x  = D1end*4.5, y = D2end*4.5, label = metric), 
                    color = "blue", inherit.aes = F , force = 0.6) + 
  geom_text_repel(data = PCA_metric_data_E2, aes (x = Dim.1, y = Dim.2, label = PlantSpeciesfull), fontface = "italic", inherit.aes = F) +
  theme_minimal() + 
  xlab(label = "PC1 (63.7 %)") +
  ylab ("PC2 (22.0 %)") + 
  guides (color= guide_legend( "Plant family")) +
  theme (legend.position = "bottom") +
  scale_color_manual(values=c("Asteraceae"= "#edc948" ,"Cyperaceae" ="#5A6351" ,
                              "Fabaceae" = "#f28e2b" , 
                              "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
                              "Poaceae" = "#7F9A65" )) 

plot1 <- fviz_contrib(PCA_all_metrics_E2, choice = "var", axes = 1, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC1")
plot2 <- fviz_contrib(PCA_all_metrics_E2, choice = "var", axes = 2, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC2")
plot3 <- fviz_contrib(PCA_all_metrics_E2, choice = "var", axes = 3, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC3")

ggarrange (plotPCA_E2,
           ggarrange (plot1, plot2, plot3, nrow = 3, ncol= 1, labels= c("B", "C", "D")), 
           nrow=1, ncol=2, widths =  c(2, 1), labels= "A") 


# PCA E1 per sample  (supplementary material) #####

# All_Metrics_E1_sample <- readRDS ("data/All_metrics_E1_sample.rds")
All_Metrics_E1_sample_df <- 
  All_metrics_E1_samples %>% 
  relocate(where(is.numeric), .after = where(is.character)) %>%  data.frame(row.names = "sampleID")

## calculation of PCA ##
PCA_all_metrics_E1_sample <- PCA(All_Metrics_E1_sample_df, quali.sup = c(1:4),   scale.unit = T, graph = F)

#PCA plot Exp 2 -  for each sample #
PCA_arrows_metrics_sample<-
  PCA_all_metrics_E1_sample$var$coord %>%  
  as_tibble (rownames = "metric") %>% 
  dplyr::rename("D1end" = "Dim.1", "D2end"= "Dim.2", "D3end" = "Dim.3")

PCA_metric_data_sample   <- PCA_all_metrics_E1_sample$ind$coord  %>%  as_tibble(rownames = "sampleID") %>% 
  left_join(metaM0 %>%  select (sampleID, PlantSpeciesfull, PlantFamily) )

saveRDS(PCA_metric_data_sample, "results/PCA_metric_data_sample.rds")

PCA_eig_m_Dim1 <- round (PCA_all_metrics_E1_sample$eig[1,2],1)
PCA_eig_m_Dim1
PCA_eig_m_Dim2 <- round (PCA_all_metrics_E1_sample$eig[2,2],1)
PCA_eig_m_Dim2
PCA_eig_m_Dim3 <- round (PCA_all_metrics_E1_sample$eig[3,2],1)
PCA_eig_m_Dim3


#Plot PCA E1 per sample - dimensions 1 and 2 #
plotPca_dim1_2  <- 
  PCA_metric_data_sample %>%  
  mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
  ggplot (aes (x = Dim.1, y = Dim.2, color = PlantSpeciesfull)) + 
  geom_point (size =3) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_segment(data = PCA_arrows_metrics_sample, aes (x=0, xend= D1end*4.5, y = 0, yend = D2end*4.5), 
               arrow = arrow(length = unit(0.3, "picas")), color = "blue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics_sample, aes (x  = D1end*4.5, y = D2end*4.5, label = metric), 
                    color = "blue", inherit.aes = F , force = 0.6) + 
  theme_minimal() + 
  xlab(label = "PC1 (70.4 %)") +
  ylab ("PC3 (12.0 %)") + 
  guides (color= guide_legend( "Plant species")) +
  theme (legend.position = "bottom") 

# Plot contributions to PCs ##
plot1 <- fviz_contrib(PCA_all_metrics_E1_sample, choice = "var", axes = 1, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC1")
plot2 <- fviz_contrib(PCA_all_metrics_E1_sample, choice = "var", axes = 2, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC2")
plot3 <- fviz_contrib(PCA_all_metrics_E1_sample, choice = "var", axes = 3, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC3")



# Plot PCA E1 per sample - dimensions 2 and 3 ##
plotPca_dim2_3  <- 
  PCA_metric_data_sample %>%  
  mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
  ggplot (aes (x = Dim.2, y = Dim.3, color = PlantSpeciesfull)) + 
  geom_point (size =3) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_segment(data = PCA_arrows_metrics_sample, aes (x=0, xend= D2end*4.5, y = 0, yend = D3end*4.5), 
               arrow = arrow(length = unit(0.3, "picas")), color = "blue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics_sample, aes (x  = D2end*4.5, y = D3end*4.5, label = metric), 
                    color = "blue", inherit.aes = F , force = 0.6) + 
  theme_minimal() + 
  xlab(label = "PC2 (12.0 %)") +
  ylab ("PC3 (7.9 %)") + 
  guides (color= guide_legend( "Plant species")) +
  theme (legend.position = "bottom") 

# Full plot 
ggarrange (ggarrange (plotPca_dim1_2,plotPca_dim2_3, nrow = 1, ncol = 2, labels = c ("A", "B"), common.legend = T),
           ggarrange (plot1, plot2, plot3, nrow = 1, ncol= 3, labels= c( "C", "D", "E")), 
           nrow=2, ncol=1, heights  =  c(2, 1)) 




