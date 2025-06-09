library(readr)
library(dplyr)
library(tidyr)


# Combine DMS data
data_dir <- "substitutions"
combined_dms <- c()
for (fn in list.files(data_dir, full.names=T)) {
  dms <- read_csv(fn, show_col_types=F) %>%
    select(c("mutant", "mutated_sequence", "DMS_score")) %>%
    mutate(name=substring(fn, nchar(data_dir)+2, nchar(fn)-4))
  combined_dms <- rbind(combined_dms, dms)
}

# Clean up and format combined DMS data
df <- combined_dms %>%
  select(-mutated_sequence) %>%
  filter(!grepl(":", mutant, fixed=T),
         !is.na(mutant)) %>%
  separate(mutant, into=c("wildtype", "position", "mutation"), 
           sep=c(1, -1), convert=T) %>%
  group_by(name) %>%
  rename(score=DMS_score) %>%
  arrange(name, position, wildtype, mutation) %>%
  mutate(deleteriousness_rank=rank(-score, ties.method="first")) %>%
  mutate(deleteriousness=(deleteriousness_rank-1)/(length(score)-1))
write_csv(df, "substitutions.csv")
