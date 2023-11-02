### Mixed models
#after Zuur 2009
# pages 120 etc, Ch 5


### get the data for all variables, the dependent variable here is total BActeria , measured with PLFA

BacttoAMF  <-   PL_FA_conc_Soil  %>%
  group_by(sampleID, ID,  FvsB) %>%   #
  dplyr::summarise (sum_groups=sum(conc_FA))  %>% #result in microbial groups
  ungroup () %>%
  dplyr::filter (FvsB != "IS")  %>%  #removes  IS
  dplyr::filter (FvsB != "NAFA")  %>%  #removes all NAFA that are of unknown origin or IS
  dplyr::filter (FvsB != "remove") %>%  
  filter (FvsB == "Bacteria") %>% 
  left_join (AMF_conc  %>%  left_join(meta_M1Wcontrol) %>% 
               mutate (fullAMF = AMF* DW_roots) %>% #AMF biomarker amount for total root system
               bind_rows(AMF_conc_Soil %>% 
                           left_join (meta_M1Wcontrol) %>% 
                           mutate (fullAMF = AMF * 1200)) %>% group_by (ID ) %>%  summarize (fullAMF=sum(fullAMF))) %>% 
  mutate (totalBacteria = 1200 * sum_groups) %>%  left_join(meta_M1Wcontrol %>% filter (roots_soil=="roots"), by = "ID")


## use package nlme, follow steps in zuur 2009
# LMER full


totalPLFA %>%  left_join (dataNLrootstosoil %>%  select (pot, AMF, AMFroot), by = "pot")







#LM to test normality and heterosc on the residuals

## AMF = 
lmer1 <- lmer (PLFA ~AMF +   DW_roots* DW_above+(1|PlantFamily/PlaSpe),  data = test)
lmer2 <-  lmer (PLFA ~  AMF  + DW_roots * DW_above + (1|PlaSpe),  data = test)
lmer3 <-  lmer (PLFA ~  AMF  + DW_roots * DW_above + (1|PlantFamily),  data = test)
lmLmer1 <-  gls (PLFA ~  AMF  + DW_roots * DW_above ,  data = test, method = "REML")
lmer4 <-  lmer (PLFA ~ AMF +  DW_above+ (DW_roots|PlaSpe),  data = test)
lmer5 <-  lmer (PLFA ~  AMF  +DW_roots  +  (DW_above|PlantFamily/PlaSpe),  data = test)
lmer6 <-  lmer (PLFA ~ AMF +  DW_above+ (DW_roots|PlantFamily/PlaSpe),  data = test)

# Which model is optimal (random structur?)
anova (lmer1, lmer2, lmer3, lmLmer1, lmer4, lmer5, lmer6)
AIC (lmer1, lmer2, lmer3, lmLmer1)


# now check fixed effects and remove effects that are not ignificant 

summary (lmer3)

# remove AMF
lmerEnd <- lmer (PLFA ~   DW_roots * DW_above + (1|PlantFamily),  data = test)




### only bacteria instead of all PLFA

Bact_to_AMF_PLFA <-   PL_FA_conc_Soil  %>%
  group_by(sampleID, ID,  FvsB) %>%   
  filter(Biom2Soil != "Fungi") %>% #
  dplyr::summarise (sum_groups=sum(conc_FA))  %>% #result in microbial groups
  ungroup () %>%
  dplyr::filter (FvsB != "IS")  %>%  #removes  IS
  dplyr::filter (FvsB != "NAFA")  %>%  #removes all NAFA that are of unknown origin or IS
  dplyr::filter (FvsB != "remove") %>%  
  pivot_wider(names_from = "FvsB", values_from = "sum_groups") %>% 
  left_join(dataNLrootstosoil, by = "") %>% 
  filter (!is.na(Fungi)) 


lmerB1 <- lmer (Bacteria ~Fungi +  DW_roots * DW_above+(1|PlantFamily/PlaSpe),  data = Bact_to_AMF_PLFA , REML = T)
lmerB2<-  lmer (Bacteria ~Fungi  + DW_roots * DW_above + (1|PlaSpe),  data = Bact_to_AMF_PLFA  , REML = T)
lmerB3 <-  lmer  (Bacteria ~Fungi   + DW_roots * DW_above + (1|PlantFamily),  data = Bact_to_AMF_PLFA , REML = T)
lmLmerB1 <-  gls  (Bacteria ~Fungi   + DW_roots * DW_above ,  data = Bact_to_AMF_PLFA , method = "REML")
lmerB4 <-  lmer (Bacteria ~  Fungi +  DW_above+ (DW_roots|PlaSpe),  data = Bact_to_AMF_PLFA )
lmerB5 <-  lmer  (Bacteria ~Fungi   +DW_roots  +  (DW_above|PlantFamily/PlaSpe),  data = Bact_to_AMF_PLFA )
lmerB6 <-  lmer (Bacteria ~ Fungi +  DW_above+ (DW_roots|PlantFamily/PlaSpe),  data = Bact_to_AMF_PLFA )

lmerB7 <- lmer  (Bacteria ~Fungi   + DW_roots * DW_above + (DW_roots * DW_above|PlantFamily),  data = Bact_to_AMF_PLFA , REML = T)


anova (lmerB1, lmerB2, lmerB3, lmLmerB1, lmerB4, lmerB5, lmerB6)

# b3 is best random 


# best model: random effect 1|PlantFamily

# AIC
AIC (lmerBtoAMF1, lmerBtoAMF2, lmerBtoAMF3, lmerBtoAMF4, lmerBtoAMF5, lmerBtoAMF6)

# step- check summary 
summary (lmerB3)


LMM_AMF_PLFA  <- step (lmerB3)





## NL roots and soil

 dataNLrootstosoil  <- 
   AMF_conc %>%
   select (ID, AMF)  %>%  #AMF is the biomarker in µmol per g substrate
   dplyr:: rename ("AMFroot" = "AMF") %>% #renamed AMF for the root samples to calculate root to soil later
   left_join ((AMF_conc_Soil %>%  ungroup () %>% select (ID, AMF)), by = "ID") %>%
   #filter (!is.na (AMFroot)) %>% filter (!is.na (AMF)) %>%  
   left_join (meta) %>%  mutate (totalAMF = AMFroot * DW_roots)

## root AMF biomarker 
mdl3 <-   lm (AMF~  AMFroot , data = dataNLrootstosoil)
lm_PlaFam <-   lm (PLFA~  PlantFamily* PlaSpe %in% PlantFamily *DW_roots*DW_above , data = totalPLFA)# that is plaSpe nested in PlantFamily 
lm_PlaFam <-   lm (PLFA~  fullAMF  + PlantSpeciesfull*DW_roots* DW_above , data = test)# that is plaSpe nested in PlantFamily 

selectedmdl1  <- step (lm_PlaFam)

summary (selectedmdl1)

summary(mdl3)
plot (mdl1, select = c(1))

mdl2 <-   summary (lm (AMFroot~  PlaSpe , data = dataNLrootstosoil))

mdl3 <-   lm (PLFA~  PlantFamily+ DW_roots, data = totalPLFA)# 

lmerAMF1 <- lmer  (totalAMF ~  DW_above  + (1|PlantFamily), data = dataNLrootstosoil)
lmerAMF2 <- lmer  (totalAMF ~  DW_above  + (1|PlantFamily), data = dataNLrootstosoil)
lmerAMF4 <- lmer  (totalAMF ~  DW_above + (1|PlaSpe), data = dataNLrootstosoil)











lmerAMF1 <- lmer  (AMFroot ~  DW_roots*DW_above*PlaSpe + (1+ DW_roots|PlaSpe) + (1+ DW_above|PlaSpe ), data = dataNLrootstosoil)
lmerAMF2 <- lmer  (AMFroot ~  DW_roots*DW_above*PlaSpe + (1+ DW_roots|PlaSpe)  , data = dataNLrootstosoil)
lmerAMF3 <- lmer  (AMFroot ~  DW_roots*DW_above*PlaSpe +   (1+ DW_above|PlaSpe ), data = dataNLrootstosoil)
lmerAMF4 <- lmer  (AMFroot ~  DW_roots*DW_above*PlaSpe + (1|PlaSpe)  , data = dataNLrootstosoil)
lmerAMF3 <- lmer  (AMFroot ~  DW_roots*DW_above*PlaSpe +   (1+ DW_above|PlaSpe ), data = dataNLrootstosoil)

lmerAMF2 <- lmer  (AMFroot ~  DW_roots*PlaSpe + (1+ DW_roots|PlaSpe)  , data = dataNLrootstosoil)
lmerAMF3 <- lmer  (AMFroot ~  DW_above*PlaSpe +   (1+ DW_above|PlaSpe ), data = dataNLrootstosoil)

#max likelihood estimations
anova ( lmerAMF2, lmerAMF4)

##