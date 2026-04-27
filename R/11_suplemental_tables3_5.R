library(tidyverse)

sp <- c('Histiotus magellanicus', 
        'Lasiurus varius', 
        'Myotis chiloensis')

csv_tally <- list.files('outputs/montecarlo/', '.csv', 
                        recursive = T, full.names = T)

csv_fls <- csv_tally[ grep('all/Summary.csv', csv_tally)]

csv_in <- csv_fls[c(1,3,4)]
names(csv_in) <- sp

summaries <- imap(csv_in, \(x, n) {
  read_csv(x) %>% mutate(species = n, .before=1)
}) %>% list_rbind()

write_csv(summaries, 'table_s3_montecarlo_variables_and_runs.csv')

all_mods_metrics <- list.files('outputs/models/fitting/', 'block_results.csv', 
                               recursive = T, full.names = T)

names(all_mods_metrics) <- sp

mods_summaries <- imap(all_mods_metrics, \(x, n) {
  read_csv(x) %>% mutate(species = n, .before=1)
}) %>% list_rbind()

write_csv(mods_summaries, 'table_s5_all_mods_tuning_vals.csv')
