
library (FactoMineR)
library(factoextra)
library(ggrepel)
library (ggpubr)



#PCA Exp 1 ####

PCA_all_metrics_E1 <- PCA(All_metrics_E1_8Sp_df, quali.sup = c(2,9),   scale.unit = T, graph = F)

PCA_arrows_metrics<-
  PCA_all_metrics_E1$var$coord %>%  
  as_tibble (rownames = "metric") %>% 
  dplyr::rename("D1end" = "Dim.1", "D2end"= "Dim.2")

PCA_metric_data   <- PCA_all_metrics_E1$ind$coord  %>%  as_tibble(rownames = "PlantSpeciesfull") %>% 
  left_join(metaM0 %>%  select (PlantSpeciesfull, PlantFamily) %>%  unique ())

PCA_eig_m_Dim1 <- round (PCA_all_metrics_E1$eig[1,2],1)
PCA_eig_m_Dim1
PCA_eig_m_Dim2 <- round (PCA_all_metrics_E1$eig[2,2],1)
PCA_eig_m_Dim2

#Plot  PCA E1 ####

plotPca  <- 
  PCA_metric_data %>%  
  mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
  ggplot (aes (x = Dim.1, y = Dim.2, color = PlantFamily)) + 
  geom_point (size =3) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  #stat_ellipse(aes (x= Dim.1, y = Dim.2, color = group))  +   # thisis 95%confidence
  geom_segment(data = PCA_arrows_metrics, aes (x=0, xend= D1end*4.5, y = 0, yend = D2end*4.5), 
               arrow = arrow(length = unit(0.3, "picas")), color = "blue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics, aes (x  = D1end*4.5, y = D2end*4.5, label = metric), 
                    color = "blue", inherit.aes = F , force = 0.6) + 
  geom_text_repel(data = PCA_metric_data, aes (x = Dim.1, y = Dim.2, label = PlantSpeciesfull), fontface = "italic", inherit.aes = F) +
  theme_minimal() + 
  xlab(label = "PC1 (88 %)") +
  ylab ("PC2 (6 %)") + 
  guides (color= guide_legend( "Plant family")) +
  theme (legend.position = "bottom") +
  scale_color_manual(values=c("Asteraceae"= "#edc948" ,"Cyperaceae" ="#5A6351" ,
                              "Fabaceae" = "#f28e2b" , 
                              "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
                              "Poaceae" = "#7F9A65" )) 



# plotPca  <- fviz_pca_biplot(PCA_all_metrics,axes = c(1,2), 
#                             col.var = "contrib", 
#                             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
#                             repel =T)

plot1 <- fviz_contrib(PCA_all_metrics_E1, choice = "var", axes = 1, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC1")
plot2 <- fviz_contrib(PCA_all_metrics_E1, choice = "var", axes = 2, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC2")
plot3 <- fviz_contrib(PCA_all_metrics_E1, choice = "var", axes = 3, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC3")

ggarrange (plotPca,
           ggarrange (plot1, plot2, plot3, nrow = 3, ncol= 1, labels= c("B", "C", "D")), 
           nrow=1, ncol=2, widths =  c(2, 1), labels= "A") 



# Exp 2 - all metrics ####


##calc of eigenvalues for PC1 for text
PCA_all_metrics_E2 <- PCA(All_metrics_E2_df, quali.sup = c(8,9,10),   scale.unit = T, graph = F)

#PCA Exp 2 ####


PCA_arrows_metrics<-
  PCA_all_metrics_E2$var$coord %>%  
  as_tibble (rownames = "metric") %>% 
  dplyr::rename("D1end" = "Dim.1", "D2end"= "Dim.2")

PCA_metric_data   <- PCA_all_metrics_E2$ind$coord  %>%  as_tibble(rownames = "PlantSpeciesfull") %>% 
  left_join(metaM0 %>%  select (PlantSpeciesfull, PlantFamily) %>%  unique ())


PCA_eig_m_Dim1 <- round (PCA_all_metrics_E2$eig[1,2],1)
PCA_eig_m_Dim1
PCA_eig_m_Dim2 <- round (PCA_all_metrics_E2$eig[2,2],1)
PCA_eig_m_Dim2
#Plot 

plotPca  <- 
  PCA_metric_data %>%  
  mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
  ggplot (aes (x = Dim.1, y = Dim.2, color = PlantFamily)) + 
  geom_point (size =3) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  #stat_ellipse(aes (x= Dim.1, y = Dim.2, color = group))  +   # thisis 95%confidence
  geom_segment(data = PCA_arrows_metrics, aes (x=0, xend= D1end*4.5, y = 0, yend = D2end*4.5), 
               arrow = arrow(length = unit(0.3, "picas")), color = "blue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics, aes (x  = D1end*4.5, y = D2end*4.5, label = metric), 
                    color = "blue", inherit.aes = F , force = 0.6) + 
  geom_text_repel(data = PCA_metric_data, aes (x = Dim.1, y = Dim.2, label = PlantSpeciesfull), fontface = "italic", inherit.aes = F) +
  theme_minimal() + 
  xlab(label = "PC1 (74 %)") +
  ylab ("PC2 (17 %)") + 
  guides (color= guide_legend( "Plant family")) +
  theme (legend.position = "bottom") +
  scale_color_manual(values=c("Asteraceae"= "#edc948" ,"Cyperaceae" ="#5A6351" ,
                              "Fabaceae" = "#f28e2b" , 
                              "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
                              "Poaceae" = "#7F9A65" )) 



# plotPca  <- fviz_pca_biplot(PCA_all_metrics,axes = c(1,2), 
#                             col.var = "contrib", 
#                             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
#                             repel =T)

plot1 <- fviz_contrib(PCA_all_metrics_E2, choice = "var", axes = 1, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC1")
plot2 <- fviz_contrib(PCA_all_metrics_E2, choice = "var", axes = 2, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC2")
plot3 <- fviz_contrib(PCA_all_metrics_E2, choice = "var", axes = 3, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC3")

ggarrange (plotPca,
           ggarrange (plot1, plot2, plot3, nrow = 3, ncol= 1, labels= c("B", "C", "D")), 
           nrow=1, ncol=2, widths =  c(2, 1), labels= "A") 


##Procrustes analysis between PCAs#####


## Needs 
# symmetry = T - both solutions are scaled to#
# unit variance
# results in procrustes m2



procr <- procrustes (PCA_all_metrics_E1$ind$coord, PCA_all_metrics_E2$ind$coord, symmetric = T)
#choice = c(1,2)  )  # with choice the number of dimensions is chosen
summary (procr)


# Plotting procrustes results##

plot (procr, kind =2)
#Kind 2 plots show the residuals for each sample. 
#This allows identification of samples with the worst fit. 
#The horizontal lines, from bottom to top, are 
#the 25% (dashed), 50% (solid), and 75% (dashed) quantiles of the residuals.

##Test of Significance

#Function protest is a permutational test of the significance of the procrustes result.
#It is based on the correlation from a symmetric Procrustes analysis.


procrtest <- protest(X = PCA_all_metrics_E1$ind$coord , Y = PCA_all_metrics_E2$ind$coord , 
                     scores = "sites", permutations = 99999)



## Procrustes protest raw data  #####
procr_test<- 
  All_metrics_E1_8Sp %>%  select (-Exp, -PlantFamily) %>% 
  #mutate_at (vars (-PlantSpeciesfull), rescale) %>%   ## in case scaled beforehand?? Like in the radar plots , would be measured relative to other plants
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


# Figure: PCA sample wise #####

#All_Metrics_E2_sample <- readRDS ("data/All_metrics_E2_sample.rds")
All_Metrics_E2_sample_df <- 
  All_Metrics_E2_sample %>% 
  relocate(where(is.numeric), .after = where(is.character)) %>%  data.frame(row.names = "sampleID")

##calc of PCA 
PCA_all_metrics_E2_sample <- PCA(All_Metrics_E2_sample_df, quali.sup = c(1:4),   scale.unit = T, graph = F)

#PCA plot Exp 2 -  for each sample #


PCA_arrows_metrics_sample<-
  PCA_all_metrics_E2_sample$var$coord %>%  
  as_tibble (rownames = "metric") %>% 
  dplyr::rename("D1end" = "Dim.1", "D2end"= "Dim.2", "D3end" = "Dim.3")

PCA_metric_data_sample   <- PCA_all_metrics_E2_sample$ind$coord  %>%  as_tibble(rownames = "sampleID") %>% 
  left_join(meta_M1 %>%  select (sampleID, PlantSpeciesfull, PlantFamily) )

saveRDS(PCA_metric_data_sample, "results/PCA_metric_data_sample.rds")

PCA_eig_m_Dim1 <- round (PCA_all_metrics_E2_sample$eig[1,2],1)
PCA_eig_m_Dim1
PCA_eig_m_Dim2 <- round (PCA_all_metrics_E2_sample$eig[2,2],1)
PCA_eig_m_Dim2
PCA_eig_m_Dim3 <- round (PCA_all_metrics_E2_sample$eig[3,2],1)
PCA_eig_m_Dim3
#Plot 


plotPca_dim1_2  <- 
  PCA_metric_data_sample %>%  
  mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
  ggplot (aes (x = Dim.1, y = Dim.2, color = PlantSpeciesfull)) + 
  geom_point (size =3) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  #stat_ellipse(aes (x= Dim.1, y = Dim.2, color = group))  +   # thisis 95%confidence
  geom_segment(data = PCA_arrows_metrics_sample, aes (x=0, xend= D1end*4.5, y = 0, yend = D2end*4.5), 
               arrow = arrow(length = unit(0.3, "picas")), color = "blue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics_sample, aes (x  = D1end*4.5, y = D2end*4.5, label = metric), 
                    color = "blue", inherit.aes = F , force = 0.6) + 
  #geom_text_repel(data = PCA_metric_data_sample, aes (x = Dim.1, y = Dim.2, label = sampleID), 
  #                fontface = "italic", inherit.aes = F) +
  theme_minimal() + 
  xlab(label = "PC1 (61 %)") +
  ylab ("PC3 (17 %)") + 
  guides (color= guide_legend( "Plant species")) +
  theme (legend.position = "bottom") 
#+
 # # scale_color_manual(values=c("Asteraceae"= "#edc948" ,"Cyperaceae" ="#5A6351" ,
 #                              "Fabaceae" = "#f28e2b" , 
 #                              "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
 #                              "Poaceae" = "#7F9A65" )) 


plot1 <- fviz_contrib(PCA_all_metrics_E2_sample, choice = "var", axes = 1, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC1")
plot2 <- fviz_contrib(PCA_all_metrics_E2_sample, choice = "var", axes = 2, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC2")
plot3 <- fviz_contrib(PCA_all_metrics_E2_sample, choice = "var", axes = 3, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC3")






plotPca_dim2_3  <- 
  PCA_metric_data_sample %>%  
  mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
  ggplot (aes (x = Dim.2, y = Dim.3, color = PlantSpeciesfull)) + 
  geom_point (size =3) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  #stat_ellipse(aes (x= Dim.1, y = Dim.2, color = group))  +   # thisis 95%confidence
  geom_segment(data = PCA_arrows_metrics_sample, aes (x=0, xend= D2end*4.5, y = 0, yend = D3end*4.5), 
               arrow = arrow(length = unit(0.3, "picas")), color = "blue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics_sample, aes (x  = D2end*4.5, y = D3end*4.5, label = metric), 
                    color = "blue", inherit.aes = F , force = 0.6) + 
  #geom_text_repel(data = PCA_metric_data_sample, aes (x = Dim.2, y = Dim.3, label = sampleID), 
  #                fontface = "italic", inherit.aes = F) +
  theme_minimal() + 
  xlab(label = "PC2 (17 %)") +
  ylab ("PC3 (13 %)") + 
  guides (color= guide_legend( "Plant species")) +
  theme (legend.position = "bottom") 

ggarrange (ggarrange (plotPca_dim1_2,plotPca_dim2_3, nrow = 1, ncol = 2, labels = c ("A", "B"), common.legend = T),
           ggarrange (plot1, plot2, plot3, nrow = 1, ncol= 3, labels= c( "C", "D", "E")), 
           nrow=2, ncol=1, heights  =  c(2, 1)) 
