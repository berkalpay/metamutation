library(readr)
library(dplyr)
library(tidyr)
library(stringr)


df_substitutions <- read_csv("substitutions.csv")
df_structure <- read_csv("structural_features.csv") %>%
  rename(pdb_name=name)
df_dms <- read_csv("dms.csv") %>%
  rename(name=DMS_id) %>%
  mutate(multirange=grepl("|", pdb_file, fixed=T)) %>%
  mutate(pdb_file=ifelse(multirange, paste0(UniProt_ID, ".pdb"), pdb_file),
         pdb_range=ifelse(multirange, paste0(str_split_i(pdb_range, "-", i=1), "-n"), pdb_range)) %>%
  separate_wider_delim(pdb_range, "-", names=c("pdb_start_position", NA)) %>%
  mutate(pdb_start_position=as.numeric(pdb_start_position)) %>%
  separate_wider_delim(pdb_file, ".", names=c("pdb_name", NA)) %>%
  select(name, pdb_name, pdb_start_position)

df_data <- df_substitutions %>%
  inner_join(df_dms, by="name") %>%
  mutate(position=position-pdb_start_position+1) %>%
  inner_join(df_structure, by=c("pdb_name", "position")) %>%
  mutate(position=position+pdb_start_position-1) %>%
  select(-c(pdb_name, pdb_start_position))
write_csv(df_data, "data.csv")
