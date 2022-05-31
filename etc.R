library(here)
library(tidyverse)

load(here("data", "RNAsamples_UGA.Rda"))

RNAsamples.UGA <- RNAsamples.melt %>% 
  filter(Metric == "TPM_g") %>%
  filter(Condition == "CON") %>%
  group_by(Transcript_ID, Stage) %>%
  summarize(Whole_Expression = mean(Expression))

df <- RNAsamples.UGA %>% 
  pivot_wider(names_from = Stage, values_from = Whole_Expression) %>% 
  mutate(type = case_when(
    AF > 5 & AM < 5 & MF < 5 ~ 'Female only',
    AF < 5 & AM > 5 & MF < 5 ~ 'Male only',
    AF > 5 & AM > 5 & MF < 5 ~ 'Adult only',
    AF > 5 & AM < 5 & MF > 5 ~ 'Female and MF',
    AF < 5 & AM > 5 & MF > 5 ~ 'Male and MF',
    AF < 5 & AM < 5 & MF > 5 ~ 'MF only',
    AF < 5 & AM < 5 & MF < 5 ~ 'None',
    AF > 5 & AM > 5 & MF > 5 ~ 'All'
  )) %>% 
  group_by(type) %>% 
  tally()
