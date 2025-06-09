library(glue)
library(stringr)
library(forcats)
library(readr)
library(dplyr)
library(tidyr)
library(rstan)
options(mc.cores = parallel::detectCores())


# Format data
aas <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
aa_to_num <- setNames(1:length(aas), aas)
num_to_aa <- setNames(aas, 1:length(aas))
df <- read_csv("data/data.csv") %>%
  mutate(wildtype=aa_to_num[wildtype],
         mutation=aa_to_num[mutation],
         p=setNames(1:length(unique(name)), unique(name))[name])
format_data <- function(df) {
  P <- length(unique(df$p))
  N <- (df %>% group_by(p) %>% count())$n
  w <- matrix(0, nrow=P, ncol=max(N))
  m <- matrix(0, nrow=P, ncol=max(N))
  surf <- matrix(0, nrow=P, ncol=max(N))
  pii <- 0
  for (pi in unique(df$p)) {
    pii <- pii + 1
    dfp <- df %>% filter(p==pi) %>% arrange(-deleteriousness)
    padding <- rep(0, max(N)-nrow(dfp))
    w[pii,] <- c(dfp$wildtype, padding)
    m[pii,] <- c(dfp$mutation, padding)
    surf[pii,] <- c(dfp$surface_accessibility, padding)
  }
  w <- apply(w, c(1, 2), as.integer)
  m <- apply(m, c(1, 2), as.integer)
  surf <- apply(surf, c(1, 2), as.numeric)
  list(P=P, N=N, w=w, m=m, surf=surf)
}

fit_aa <- function(df) {
  with(format_data(df),
       optimizing(stan_model(file="models/aa.stan"),
                  data=list(P=P, N=N, maxN=max(N), w=w, m=m),
                  iter=10**5, seed=42, verbose=T))
}
fit_full <- function(df) {
  with(format_data(df),
       optimizing(stan_model(file="models/full.stan"),
                  data=list(P=P, N=N, maxN=max(N),
                            w=w, m=m, surf=surf),
                  iter=10**5, seed=42, verbose=T))
}


# Exchangeability model
aa_fit <- fit_aa(df)
unpack_aa <- function(fit) {
  tibble(entry=names(fit$par), score=fit$par) %>%
    mutate(w=as.integer(str_extract(entry, "(?<=\\[)\\d+")),
           m=as.integer(str_extract(entry, "(?<=,)\\d+(?=\\])"))) %>%
    mutate(wildtype=num_to_aa[w],
           mutation=num_to_aa[m]) %>%
    filter(wildtype!=mutation) %>%
    select(wildtype, mutation, score) %>%
    arrange(wildtype, mutation)
}
write_csv(unpack_aa(aa_fit), "results/aa.csv")


# Exchangeability and surface accessibility model
full_fit <- fit_full(df)
unpack_full <- function(fit) {
  tibble(entry=names(fit$par), score=fit$par) %>%
    mutate(w=as.integer(str_extract(entry, "(?<=\\[)\\d+")),
           m=as.integer(str_extract(entry, "(?<=,)\\d+"))) %>%
    mutate(wildtype=num_to_aa[w],
           mutation=num_to_aa[m]) %>%
    filter(wildtype!=mutation,
           grepl("[", entry, fixed=T)) %>%
    separate_wider_delim(entry, "[", names=c("param", NA)) %>%
    pivot_wider(names_from=param, values_from=score) %>%
    select(-c(w, m)) %>%
    arrange(wildtype, mutation)
}
write_csv(unpack_full(full_fit) %>% mutate(sa=full_fit$par[["sa_mu"]]),
          "results/full.csv")


# Split training DMS by type
splits_df <- read_csv("data/dms.csv") %>%
  rename(name=DMS_id) %>%
  group_by(coarse_selection_type) %>%
  summarize(names=list(name))
coarse_selection_types <- splits_df$coarse_selection_type
splits <- splits_df[[2]]
if (!dir.exists("results/split_params_aa/")) dir.create("results/split_params_aa/")
for (i in 1:5) {
  print(i)
  aa_fit <- fit_aa(df %>% filter(name %in% splits[[i]]))
  write_csv(unpack_aa(aa_fit),
            glue("results/split_params_aa/{coarse_selection_types[i]}.csv"))
}


# Cross-validation
df_cv <- df %>%
  inner_join(read_csv("data/dms.csv") %>%
               rename(name=DMS_id) %>%
               select(name, UniProt_ID),
             by="name")
n_splits <- 10
set.seed(42)
P <- length(unique(df_cv$UniProt_ID))
splits <- split(sample(unique(df_cv$UniProt_ID)),
                ceiling(seq_along(sample(1:P))/(P/n_splits)))
pred <- c()
for (i in 1:n_splits) {
  print(i)
  df_train <- df_cv %>% filter(!(UniProt_ID %in% splits[[i]]))
  df_test <- df_cv %>% filter(UniProt_ID %in% splits[[i]])

  print("AA-only")
  scores_aa <- unpack_aa(fit_aa(df_train))
  print("Full")
  full_fit <- fit_full(df_train)
  scores_full <- unpack_full(full_fit) %>%
    mutate(sa=full_fit$par[["sa_mu"]])

  pred_holdout <- df_test %>%
    select(wildtype, position, mutation, name, surface_accessibility) %>%
    mutate(wildtype=num_to_aa[wildtype],
           mutation=num_to_aa[mutation]) %>%
    inner_join(scores_aa, by=c("wildtype", "mutation")) %>%
    mutate(`AA-only`=1-1/(1+exp(-score))) %>%
    select(-score) %>%
    inner_join(scores_full, by=c("wildtype", "mutation")) %>%
    mutate(`Full model`=1-1/(1+exp(-(mu + sa*surface_accessibility)))) %>%
    select(-c(mu, sa, surface_accessibility))
  pred <- rbind(pred, pred_holdout)
}
write_csv(pred, "results/predictions_aa_full.csv")
