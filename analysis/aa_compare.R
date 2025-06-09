library(Biostrings)
library(readr)
library(tibble)
library(tidyr)
library(dplyr)


# Yampolsky-Stoltzfus
df_yampolsky <- read_csv("data/matrices/yampolsky.csv") %>%
  rename(wildtype=`...1`) %>%
  pivot_longer(!wildtype, names_to="mutation", values_to="score")
df_yampolsky_imputations <- tibble(wildtype=c("Y", "Y", "Y", "W", "W", "W"),
                                   mutation=c("T", "I", "V", "N", "D", "I"),
                                   imputed_score=c(mean(c(258, 293)), mean(c(258, 293)), mean(c(258, 305)),
                                                   mean(c(142, 258)), mean(c(142, 225)), mean(c(142, 293))))
df_yampolsky_imputed <- df_yampolsky %>%
  full_join(df_yampolsky_imputations, by=c("wildtype", "mutation")) %>%
  mutate(score=ifelse(is.na(score), imputed_score, score)) %>%
  select(wildtype, mutation, score) %>%
  filter(!is.na(score))
write_csv(df_yampolsky_imputed, "data/matrices/yampolsky_imputed.csv")


# BLOSUM
data(BLOSUM62)
blosum <- tibble(data.frame(BLOSUM62))
blosum$wildtype <- names(blosum)
blosum <- blosum %>%
  pivot_longer(!wildtype, names_to="mutation", values_to="score")
write_csv(blosum, "data/matrices/blosum.csv")


# Compare ranks of scores
with(read_csv("results/aa.csv") %>%
       inner_join(df_yampolsky_imputed %>% rename(score.yampolsky=score), by=c("wildtype", "mutation")) %>%
       inner_join(blosum %>% rename(score.blosum=score), by=c("wildtype", "mutation")),
     cor(score, score.yampolsky, method="spearman"))
