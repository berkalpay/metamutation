library(forcats)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(shadowtext)
library(ggtext)


df_data <- read_csv("data/data.csv")
join_cols <- c("wildtype", "position", "mutation", "name")
df_pred <- df_data %>%
  rename(DMS=score) %>%
  inner_join(read_csv("results/predictions_aa_full.csv"), by=join_cols) %>%
  mutate(`Surf. acc.`=1/(1+exp(-surface_accessibility-0.5))) %>%
  inner_join(read_csv("data/matrices/yampolsky_imputed.csv"), by=c("wildtype", "mutation")) %>%
  rename(`Yampolsky`=score) %>%
  inner_join(read_csv("data/matrices/blosum.csv"), by=c("wildtype", "mutation")) %>%
  rename(BLOSUM62=score)


# Full model
df_dms <- read_csv("data/dms.csv") %>%
  select(DMS_id, taxon, DMS_number_single_mutants, coarse_selection_type) %>%
  rename(name=DMS_id)
df_eval <- df_pred %>%
  group_by(name) %>%
  summarize(spearman=cor(DMS, `Full model`, method="spearman")) %>%
  inner_join(df_dms, by="name")
paste("mean Spearman correlation:", round(mean(df_eval$spearman), 4))
ggplot(df_eval, aes(x=DMS_number_single_mutants, y=spearman, color=taxon)) +
  geom_point(alpha=0.8, shape=21, stroke=0.6) +
  labs(x="DMS size", y="Spearman correlation") +
  theme_classic() +
  theme(legend.position="top", legend.title=element_blank(),
        aspect.ratio=1)
ggsave("figures/dms_correlations.pdf", width=4, height=4)


# Compare to other models
combine_scores <- function(score_col) {
  combined_scores <- c()
  for (fn in list.files("data/scores/", full.names=T)) {
    dms <- read_csv(fn, show_col_types=F) %>%
      mutate(score=get(score_col)) %>%
      select(c("mutant", "score")) %>%
      mutate(name=strsplit(basename(fn), "[.]")[[1]][1])
    combined_scores <- rbind(combined_scores, dms)
  }
  combined_scores %>%
    filter(!grepl(":", mutant, fixed=T),
           !is.na(mutant)) %>%
    separate(mutant, into=c("wildtype", "position", "mutation"),
             sep=c(1, -1), convert=T)
}
join_scores <- function(df_main, df_pred, model_name) {
  df_new <- df_main %>%
    inner_join(df_pred, by=join_cols)
  df_new[model_name] <- df_new$score
  select(df_new, !score)
}

df_pred <- join_scores(df_pred, combine_scores("GEMME"), "GEMME")
df_pred <- join_scores(df_pred, combine_scores("ProSST-2048"), "ProSST")
df_pred <- join_scores(df_pred, combine_scores("TranceptEVE_L"), "TranceptEVE")
df_pred <- join_scores(df_pred, combine_scores("MIF"), "MIF")
print(nrow(df_pred) == nrow(df_data))
write_csv(df_pred, "results/predictions.csv")

cors_df <- c()
score_cols <- c("DMS",
                "Full model", "AA-only", "Surf. acc.",
                "Yampolsky", "BLOSUM62",
                "MIF", "GEMME", "TranceptEVE", "ProSST")
for (score_col1 in score_cols) {
  for (score_col2 in score_cols) {
    m <- mean((df_pred %>%
                 group_by(name) %>%
                 summarize(spearman=cor(get(score_col1), get(score_col2),
                                        method="spearman")))$spearman)
    cors_df <- rbind(cors_df, c(score_col1, score_col2, m))
  }
}
cors_df <- tibble(data.frame(cors_df)) %>%
  mutate(X1=factor(X1, levels=score_cols), X2=factor(X2, levels=score_cols),
         X3=as.numeric(X3))
names(cors_df) <- c("model1", "model2", "spearman")

ggplot(cors_df %>%
         filter(as.integer(model1)<as.integer(model2)),
       aes(model1, model2)) +
  geom_tile(aes(fill=spearman), show.legend=F) +
  geom_shadowtext(aes(label=paste0(".", round(spearman, 2)*100)), size=4.2) +
  scale_x_discrete(labels = paste(c("**DMS**", score_cols[2:length(score_cols)]))) +
  scale_y_discrete(limits=rev, expand=c(0,0)) +
  scale_fill_continuous(limits=c(0, 1), low="white", high="darkgreen") +
  geom_vline(xintercept=0.5, linewidth=1) +
  geom_vline(xintercept=1.505, linewidth=1) +
  annotate("text", x=3.2, y=8.3, angle=-45, label="Physicochemical") +
  annotate("text", x=5.8, y=5.8, angle=-45, label="Matrices") +
  annotate("text", x=8.4, y=3.3, angle=-45, label="Sequence models") +
  theme_minimal() +
  theme(panel.grid.major=element_blank(), aspect.ratio=1,
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_markdown(color="black", angle=45, hjust=1),
        axis.text.y=element_text(color="black"))
ggsave("figures/rank_correlations.pdf", width=5, height=5)
