## Ch 2publication & Ch 3 data ##

## Plant & Soil Journal ###

#05/12/2022 -- correlations - microbial groups to PCA dim ##

## get PCA dimensions form new PCA ###

#Neither Spearman's rank correlation nor Pearson product-moment correlation are restricted to positive data, 
#because Spearman's just deals with how are values ordered and Pearson's deals with distance of (products of) 
#values to the mean and those work equally fine if values are negative. You can even add an arbitrary real number to any of
#your variables and it won't change any correlation.

## data PLFA ###

#AMF NLFA soil ~  PC   ####
dataNLrootstosoil %>%  # run   dataNLrootstosoil_gen
  select (sampleID, AMF, Dim.1, Dim.2, Dim.3)  %>% 
  dplyr::filter(!is.na (AMF))%>% 
  ungroup() %>% 
  select (!sampleID) %>% 
  pivot_longer(cols =  !AMF, names_to = "Dimx", values_to = "dimValue" ) %>% 
  split (.$Dimx) %>% 
  map_dfr (function(df)   cor_test(df,dimValue , AMF)) %>%  flextable ()

# no statistical correlation below 0.5

# cor = 0.3 with p = 0.09

#### when generalists plants alone, then PC2 corr = 0.4 with p = 0.04

#### check plots !!! 

# PLFA ~ PC ####

##  PLFA microbial groups ~ PC1 ####

PL_FA_conc_Soil %>%  
  as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PC_data_M1roots, by = "PlantSpeciesfull") %>% 
  group_by (Biom2Soil, sampleID) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select (sampleID, Biom2Soil, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>% 
  select (!sampleID)  %>% 
    split (.$Biom2Soil) %>% 
#map(., pivot_longer(.$FA_conc,  names_to = "Dimx", values_to = "DimValues"))
map_dfr (function(df)   cor_test(df,FA_conc , Dim.1), .id = "MG")  %>% 
  filter (MG != "IS") %>% 
  filter (MG != "NAFA")  %>% 
  filter (MG != "AMF") %>% 
  flextable () %>%  
  bold (part = "body", i = ~ p < 0.05, j = "p") %>%  
  add_header_lines (values = "Correlation PC1 ~ PL fatty acid concentration")
  
##  PLFA microbial groups ~ PC2 ####
PL_FA_conc_Soil %>%  
  as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PC_data_M1roots, by = "PlantSpeciesfull") %>% 
  group_by (Biom2Soil, sampleID) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select (sampleID, Biom2Soil, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>% 
  select (!sampleID)  %>% 
  split (.$Biom2Soil) %>% 
  #map(., pivot_longer(.$FA_conc,  names_to = "Dimx", values_to = "DimValues"))
  map_dfr (function(df)   cor_test(df,FA_conc , Dim.2), .id = "MG")  %>% 
  filter (MG != "IS") %>% 
  filter (MG != "NAFA")  %>% 
  filter (MG != "AMF") %>% 
  flextable () %>%  
  bold (part = "body", i = ~ p < 0.05, j = "p") %>%  
  add_header_lines (values = "Correlation PC2 ~ PL fatty acid concentration")


##  PLFA microbial groups ~ PC3 ####
PL_FA_conc_Soil %>%  
  as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PC_data_M1roots, by = "PlantSpeciesfull") %>% 
  group_by (Biom2Soil, sampleID) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select (sampleID, Biom2Soil, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>% 
  select (!sampleID)  %>% 
  split (.$Biom2Soil) %>% 
  #map(., pivot_longer(.$FA_conc,  names_to = "Dimx", values_to = "DimValues"))
  map_dfr (function(df)   cor_test(df,FA_conc , Dim.3), .id = "MG")  %>% 
  filter (MG != "IS") %>% 
  filter (MG != "NAFA")  %>% 
  filter (MG != "AMF") %>% 
  flextable () %>%  
  bold (part = "body", i = ~ p < 0.05, j = "p") %>%  
  add_header_lines (values = "Correlation PC3 ~ PL fatty acid concentration")


# Total PLFA  to PC 1, 2, 3 ####
## total PLFA to PC1 ####
PL_FA_conc_Soil  %>% 
as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PC_data_M1roots, by = "PlantSpeciesfull")  %>% 
  dplyr::filter (Biom2Soil != "IS" & Biom2Soil != "NAFA" & Biom2Soil != "AMF"  & Biom2Soil != "Fungi" )  %>%  
  group_by ( sampleID) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select ( sampleID, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>% 
 # select (!sampleID)  %>% 
 # split (.$Biom2Soil) %>% 
  #map(., pivot_longer(.$FA_conc,  names_to = "Dimx", values_to = "DimValues"))
  cor_test(.,FA_conc , Dim.1)

## total bacterial biomass PLFA to PC1 ####
PL_FA_conc_Soil %>% 
  as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PC_data_M1roots, by = "PlantSpeciesfull")  %>% 
  dplyr::filter (Biom2Soil != "IS" & Biom2Soil != "NAFA" & Biom2Soil != "AMF"  & Biom2Soil != "Fungi" )  %>%  
  group_by ( sampleID) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select (sampleID, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>% 
  ggplot (aes (x = Dim.1 , y = FA_conc)) + geom_point () + geom_smooth(method ="lm")
### check that, because of the differences in total soil PLFA and the extreme values between the PC dimension??

## total PLFA to PC2 ####
PL_FA_conc_Soil %>% 
as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PC_data_M1roots, by = "PlantSpeciesfull")  %>% 
  dplyr::filter (Biom2Soil != "IS" & Biom2Soil != "NAFA" & Biom2Soil != "AMF"  & Biom2Soil != "Fungi" )  %>%  
  group_by ( sampleID) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select (sampleID, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>% 
  # select (!sampleID)  %>% 
  # split (.$Biom2Soil) %>% 
  #map(., pivot_longer(.$FA_conc,  names_to = "Dimx", values_to = "DimValues"))
  cor_test(.,FA_conc , Dim.2)


## total PLFA to PC3 ####
PL_FA_conc_Soil  %>% 
as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PC_data_M1roots, by = "PlantSpeciesfull")  %>% 
  dplyr::filter (Biom2Soil != "IS" & Biom2Soil != "NAFA" & Biom2Soil != "AMF"  & Biom2Soil != "Fungi" )  %>%  
  group_by ( sampleID) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select (sampleID, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>% 
  # select (!sampleID)  %>% 
  # split (.$Biom2Soil) %>% 
  #map(., pivot_longer(.$FA_conc,  names_to = "Dimx", values_to = "DimValues"))
  cor_test(.,FA_conc , Dim.3)





## Gram pos bacteria to PC1 ####
PL_FA_conc_Soil %>% 
  as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PC_data_M1roots, by = "PlantSpeciesfull")  %>% 
  dplyr::filter (Biom2Soil != "IS" & Biom2Soil != "NAFA" & Biom2Soil != "AMF"  & Biom2Soil != "Fungi" )  %>%  
  group_by ( sampleID, Biom2Soil) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select (sampleID, Biom2Soil, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>%  
  filter (Biom2Soil == "gram.pos")  %>% 
  ggplot (aes (x = Dim.1, y = FA_conc)) +
  geom_point () + 
  geom_smooth(method = "lm")


#-----------------------------------------------------##
# using only the generalists???

# GENERALISTS PLANTS  M1 only ####
## new PCA - generalists  #####

# all interaction metrics for generalists only ##
All_metrics_M1_8Sp_gen  <-
All_metrics_M1_8Sp %>%  as_tibble () %>% 
  relocate (PlantSpeciesfull) %>% 
  filter (PlantFamily != "Asteraceae")  %>% 
  data.frame(row.names = "PlantSpeciesfull")


#PCA run with generalist data ##
PCA_all_metrics_M1_generalists <- PCA(All_metrics_M1_8Sp_gen, quali.sup = 8,   scale.unit = T, graph = F)


### PCA Plot data generalists ###

PCA_arrows_metrics_M1_gen<-
  PCA_all_metrics_M1_generalists$var$coord %>%  
  as_tibble (rownames = "metric") %>% 
  dplyr::rename("D1end" = "Dim.1", "D2end"= "Dim.2", "D3end"="Dim.3")

PCA_metric_data_M1_gen   <- PCA_all_metrics_M1_generalists$ind$coord  %>%  
  as_tibble(rownames = "PlantSpeciesfull") %>% 
  left_join(metaM0 %>%  select (PlantSpeciesfull, PlantFamily) %>%  unique ())

PCA_eig_m_Dim1_gen <- round (PCA_all_metrics_M1_generalists$eig[1,2],1)
PCA_eig_m_Dim2_gen <- round (PCA_all_metrics_M1_generalists$eig[2,2],1)

## PCA plot generalists #####


plot_PCA_M1_gen  <-
PCA_metric_data_M1_gen %>%  
  ggplot (aes (x = Dim.1, y = Dim.2, color = PlantFamily)) + 
  geom_point (size =3) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  #stat_ellipse(aes (x= Dim.1, y = Dim.2, color = group))  +   # thisis 95%confidence
  geom_segment(data = PCA_arrows_metrics_M1_gen, aes (x=0, xend= D1end*1, y = 0, yend = D2end*1), 
               arrow = arrow(length = unit(0.3, "picas")), color = "darkblue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics_M1_gen, aes (x  = D1end*1, y = D2end*1, label = metric), 
                    color = "darkblue", inherit.aes = F , force = 0.6) + 
  geom_text_repel(data = PCA_metric_data_M1_gen, aes (x = Dim.1, y = Dim.2, label = PlantSpeciesfull), fontface = "italic", inherit.aes = F) +
  theme_minimal() + 
  xlab(label = "PC1 (45.1%)") +
  ylab ("PC2 (34.5 %)") + 
  guides (color= guide_legend( "Plant family")) +
  theme (legend.position = "bottom") +
  scale_color_manual(values=c("Asteraceae"= "#edc948" ,"Cyperaceae" ="#5A6351" ,
                              "Fabaceae" = "#f28e2b" , 
                              "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
                              "Poaceae" = "#7F9A65" )) 



plot1 <- fviz_contrib(PCA_all_metrics_M1_generalists, choice = "var", axes = 1, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC1")
plot2 <- fviz_contrib(PCA_all_metrics_M1_generalists, choice = "var", axes = 2, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC2")
plot3 <- fviz_contrib(PCA_all_metrics_M1_generalists, choice = "var", axes = 3, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC3")

ggarrange (plot_PCA_M1_gen,
           ggarrange (plot1, plot2, plot3, nrow = 3, ncol= 1, labels= c("B", "C", "D")), 
           nrow=1, ncol=2, widths =  c(2, 1), labels= "A") 





## PLFA NLFA  data table for LM ####
dataNLrootstosoil_gen  <-
   PCA_all_metrics_M1_generalists$ind$coord  %>% 
  as_tibble  (rownames = "PlantSpeciesfull")  %>% 
  arrange (Dim.1) %>%  left_join(dataNLrootstosoil %>%  select (-Dim.1, -Dim.2, -Dim.3, -Dim.4, -Dim.5)) 


#AMF NLFA soil ~  PC   ####
dataNLrootstosoil_gen %>%  # run   dataNLrootstosoil_gen
  select (sampleID, AMF, Dim.1, Dim.2, Dim.3)  %>% 
  dplyr::filter(!is.na (AMF))%>% 
  ungroup() %>% 
  select (!sampleID) %>% 
  pivot_longer(cols =  !AMF, names_to = "Dimx", values_to = "dimValue" ) %>% 
  split (.$Dimx) %>% 
  map (function(df)   cor_test(df,dimValue , AMF))


#### when generalists plants alone, then PC2 corr = 0.4 with p = 0.04

## PC2 is mainly contributed by metrics: PD and unique 

#### check plots !!! 

dataNLrootstosoil_gen %>%  ggplot (aes (x = Dim.2, y = AMF)) + geom_point() + geom_smooth(method = "lm")
## is that pseudoreplication to use PC 1 5 times?

# differences 
# outlier is poa cita 

# PLFA ~ PC ####

##  PLFA microbial groups ~ PC1 ####

PL_FA_conc_Soil %>%  
  as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PCA_metric_data_M1_gen, by = "PlantSpeciesfull") %>% 
  group_by (Biom2Soil, sampleID) %>% 
  filter (!is.na(Dim.1)) %>% 
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select (sampleID, Biom2Soil, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>% 
  select (!sampleID)  %>% 
  split (.$Biom2Soil) %>% 
  #map(., pivot_longer(.$FA_conc,  names_to = "Dimx", values_to = "DimValues"))
  map_dfr (function(df)   cor_test(df,FA_conc , Dim.1), .id = "MG")  %>% 
  filter (MG != "IS") %>% 
  filter (MG != "NAFA")  %>% 
  filter (MG != "AMF") %>% 
  flextable () %>%  
  bold (part = "body", i = ~ p < 0.05, j = "p") %>%  
  add_header_lines (values = "Generalist plants -Correlation PC1 ~ PL fatty acid concentration")

##  PLFA microbial groups ~ PC2 ####
PL_FA_conc_Soil %>%  
  as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PCA_metric_data_M1_gen, by = "PlantSpeciesfull") %>% 
  group_by (Biom2Soil, sampleID) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select (sampleID, Biom2Soil, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>% 
  select (!sampleID)  %>% 
  split (.$Biom2Soil) %>% 
  #map(., pivot_longer(.$FA_conc,  names_to = "Dimx", values_to = "DimValues"))
  map_dfr (function(df)   cor_test(df,FA_conc , Dim.2), .id = "MG")  %>% 
  filter (MG != "IS") %>% 
  filter (MG != "NAFA")  %>% 
  filter (MG != "AMF") %>% 
  flextable () %>%  
  bold (part = "body", i = ~ p < 0.05, j = "p") %>%  
  add_header_lines (values = "Generalist plants Correlation PC2 ~ PL fatty acid concentration")


##  PLFA microbial groups ~ PC3 ####
PL_FA_conc_Soil %>%  
  as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PCA_metric_data_M1_gen, by = "PlantSpeciesfull") %>% 
  group_by (Biom2Soil, sampleID) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select (sampleID, Biom2Soil, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>% 
  select (!sampleID)  %>% 
  split (.$Biom2Soil) %>% 
  #map(., pivot_longer(.$FA_conc,  names_to = "Dimx", values_to = "DimValues"))
  map_dfr (function(df)   cor_test(df,FA_conc , Dim.3), .id = "MG")  %>% 
  filter (MG != "IS") %>% 
  filter (MG != "NAFA")  %>% 
  filter (MG != "AMF") %>% 
  flextable () %>%  
  bold (part = "body", i = ~ p < 0.05, j = "p") %>%  
  add_header_lines (values = "Correlation PC3 ~ PL fatty acid concentration")


# Total PLFA  to PC 1, 2, 3 ####
## total PLFA to PC1 ####
PL_FA_conc_Soil  %>% 
  as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PCA_metric_data_M1_gen, by = "PlantSpeciesfull")  %>% 
  dplyr::filter (Biom2Soil != "IS" & Biom2Soil != "NAFA" & Biom2Soil != "AMF"  & Biom2Soil != "Fungi" )  %>%  
  group_by ( sampleID) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select ( sampleID, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>% 
  # select (!sampleID)  %>% 
  # split (.$Biom2Soil) %>% 
  #map(., pivot_longer(.$FA_conc,  names_to = "Dimx", values_to = "DimValues"))
  cor_test(.,FA_conc , Dim.1)

## total PLFA to PC2 ####
PL_FA_conc_Soil %>% 
  as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PCA_metric_data_M1_gen, by = "PlantSpeciesfull")  %>% 
  dplyr::filter (Biom2Soil != "IS" & Biom2Soil != "NAFA" & Biom2Soil != "AMF"  & Biom2Soil != "Fungi" )  %>%  
  group_by ( sampleID) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select (sampleID, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>% 
  # select (!sampleID)  %>% 
  # split (.$Biom2Soil) %>% 
  #map(., pivot_longer(.$FA_conc,  names_to = "Dimx", values_to = "DimValues"))
  cor_test(.,FA_conc , Dim.2)


## total PLFA to PC3 ####
PL_FA_conc_Soil  %>% 
  as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PCA_metric_data_M1_gen, by = "PlantSpeciesfull")  %>% 
  dplyr::filter (Biom2Soil != "IS" & Biom2Soil != "NAFA" & Biom2Soil != "AMF"  & Biom2Soil != "Fungi" )  %>%  
  group_by ( sampleID) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select (sampleID, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>% 
  # select (!sampleID)  %>% 
  # split (.$Biom2Soil) %>% 
  #map(., pivot_longer(.$FA_conc,  names_to = "Dimx", values_to = "DimValues"))
  cor_test(.,FA_conc , Dim.3)

## total bacterial biomass PLFA to PC1 ####
PL_FA_conc_Soil %>% 
  as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PCA_metric_data_M1_gen, by = "PlantSpeciesfull")  %>% 
  dplyr::filter (Biom2Soil != "IS" & Biom2Soil != "NAFA" & Biom2Soil != "AMF"  & Biom2Soil != "Fungi" )  %>%  
  group_by ( sampleID) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select (sampleID, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>% 
  ggplot (aes (x = Dim.1 , y = FA_conc)) + geom_point () + geom_smooth(method ="lm")
### check that, because of the differences in total soil PLFA and the extreme values between the PC dimension??



## Gram pos bacteria to PC1 ####
PL_FA_conc_Soil %>% 
  as_tibble ()  %>% 
  left_join (meta_M1Wcontrol) %>% 
  select (!pot & !SampleName & !DW_above  & !DW_roots  & !M1_M2_order) %>% 
  left_join (PCA_metric_data_M1_gen, by = "PlantSpeciesfull")  %>% 
  dplyr::filter (Biom2Soil != "IS" & Biom2Soil != "NAFA" & Biom2Soil != "AMF"  & Biom2Soil != "Fungi" )  %>%  
  group_by ( sampleID, Biom2Soil) %>%
  mutate (FA_conc = sum (conc_FA)) %>% 
  ungroup () %>% 
  select (sampleID, Biom2Soil, FA_conc,  Dim.1, Dim.2, Dim.3) %>% 
  unique ()  %>%  
  filter (Biom2Soil == "gram.pos")  %>% 
  ggplot (aes (x = Dim.1, y = FA_conc)) +
  geom_point () + 
  geom_smooth(method = "lm")



