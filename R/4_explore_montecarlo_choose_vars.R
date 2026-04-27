library(tidyverse)

fls <- list.files('outputs/montecarlo/', 
                  'Summary.csv', full.names = T, 
                  recursive = T)

dbs <- map(fls, read_csv)


colsadd <- cbind(
  species = rep(gsub('outputs/montecarlo//|/.*', '', fls), each= 50), 
  vartype = rep(c('all', 'bio', 'env'), each=50)
)


alldb <- cbind(colsadd, do.call(rbind, dbs))                 

vardb <- alldb %>% 
  # group_by(species) %>%
  filter(Sprint == 50, 
         # vartype == 'all'
         ) %>% 
  arrange(vartype) 

vardb %>% write_csv('outputs/partial_results/montecarlo_topvars.csv')

####### OLD MONTECARLO
fls <- list.files('outputs/old_montecarlo/', 
                  'Summary.csv', full.names = T, 
                  recursive = T)

dbs <- map(fls, read_csv)


colsadd <- cbind(
  species = rep(gsub('outputs/old_montecarlo//|/.*', '', fls), each= 50), 
  vartype = rep(c('all', 'bio', 'env'), each=50)
)

alldb <- cbind(colsadd, do.call(rbind, dbs))                 

old_vardb <- alldb %>% 
  # group_by(species) %>%
  filter(Sprint == 50, 
         # vartype == 'all'
  ) %>% 
  arrange(vartype) 

rbind(vardb %>% mutate(case = 'new'), old_vardb %>% mutate(case = 'old')) %>% 
  write_csv('outputs/partial_results/compare_montecarlo_topvars.csv')



old_vardb %>% filter(vartype == 'all') %>% select(`Var 1`:last_col()) %>% unlist() %>% 
  table() %>% as.data.frame() %>% arrange(-Freq)

varsIn <- vardb %>% filter(vartype == 'all', 
                 !species  %in% c('Lasiurus varius_summer',
                                  'Lasiurus varius_winter')) %>% 
  select(`Var 1`:last_col()) %>% unlist() %>% 
  table() %>% as.data.frame() %>% arrange(-Freq)

vardb$AUC
var = 'AUC'

alldb %>% 
  group_by(species, vartype) %>% 
  # pivot_longer(
  #   cols=c(AUC, orMTP, AICc)
  # ) %>% 
  # filter(name !='AICc') %>%
  ggplot(aes(Sprint, !!ensym(var), group=vartype, colour=vartype))+
  geom_point(alpha=0.8)+
  geom_smooth()+
  labs(title = var)+
  facet_grid(species ~  vartype,scales='free')


alldb %>% 
  group_by(species, vartype) %>% 
  filter(
    str_detect(species, 'Lasiurus'), 
    Sprint == 50) %>% 
  select(species, vartype, `Var 1`:last_col()) %>% 
  arrange(vartype) %>% 
  View()
