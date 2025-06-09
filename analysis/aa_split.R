library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(glue)


aas <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
aa_map <- setNames(1:length(aas), aas)
aa_levels <- read_csv("data/aa_properties.csv") %>%
  arrange(hydropathy) %>%
  pull(aa)

# Compute rank correlations with pan-assay model
coarse_selection_types <- c("Activity", "Binding", "Expression",
                            "OrganismalFitness", "Stability")
df_params <- read_csv("results/aa.csv") %>% mutate(model="full")
for (i in coarse_selection_types) {
  fn <- as.character(glue("results/split_params_aa/{i}.csv"))
  df_params_i <- read_csv(fn) %>% mutate(model=i)
  df_params <- bind_rows(df_params, df_params_i)
}
df_full_params <- filter(df_params, model=="full")
cors <- list()
for (i in coarse_selection_types) {
  corr <- df_full_params %>%
    inner_join(filter(df_params, model==i),
               by=c("wildtype", "mutation")) %>%
    summarize(corr=cor(score.x, score.y, method="spearman")) %>%
    pull(corr)
  cors[[i]] <- corr
}

dfs <- list()
label_names <- c()
split_dir <- "results/split_params_aa/"
for (split_fn in list.files(split_dir)) {
  fn <- paste0(split_dir, split_fn)
  split_name <- strsplit(split_fn, ".csv")[[1]][1]
  corr <- cors[[split_name]]
  split_name <- ifelse(split_name=="OrganismalFitness",
                       "Organismal Fitness", split_name)
  label_names[split_name] <- paste0(split_name,
                                    " (ρ=", as.character(round(corr, 2)), ")")
  dfs[[fn]] <- read_csv(fn) %>%
    inner_join(df_full_params, by=c("wildtype", "mutation")) %>%
    mutate(score_rank_diff=rank(score.x)-rank(score.y)) %>%
    select(wildtype, mutation, score_rank_diff) %>%
    mutate(split_name=split_name)
}
df <- bind_rows(dfs)

ggplot(df %>%
         mutate(wildtype=factor(wildtype, levels=aa_levels),
                mutation=factor(mutation, levels=aa_levels))) +
  geom_tile(aes(fill=score_rank_diff, x=wildtype, y=mutation)) +
  facet_wrap(~split_name, nrow=3, scales="free",
             labeller=as_labeller(label_names)) +
  scale_x_discrete(breaks=names(aa_map)) +
  scale_y_discrete(limits=rev, breaks=names(aa_map)) +
  scale_fill_distiller(limits=c(-380+1, 380-1), breaks=c(-300, -150, 0, 150, 300),
                       palette="RdBu", direction=-1) +
  labs(x="Wildtype", y="Mutation",
       fill="Rank difference from pan-assay exchangeability") +
  guides(fill=guide_colorbar(ticks.colour=NA)) +
  theme_minimal() +
  theme(aspect.ratio=1,
        panel.grid.major=element_blank(),
        axis.text=element_text(color="black"),
        legend.position="bottom", legend.key.height=unit(0.5, "lines"),
        legend.key.width=unit(2, "lines"))
ggsave("figures/aa_split.pdf", width=5.5, height=9, device=cairo_pdf)
