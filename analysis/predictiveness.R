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
  inner_join(read_csv("data/matrices/yampolsky_imputed.csv"), by=c("wildtype", "mutation")) %>%
  rename(`Yampolsky-Stoltzfus`=score) %>%
  inner_join(read_csv("data/matrices/blosum.csv"), by=c("wildtype", "mutation")) %>%
  rename(BLOSUM62=score)


# Full model
df_dms <- read_csv("data/dms.csv") %>%
  select(DMS_id, DMS_number_single_mutants, coarse_selection_type) %>%
  rename(name=DMS_id)
df_eval <- df_pred %>%
  group_by(name) %>%
  summarize(spearman=cor(DMS, `Full model`, method="spearman")) %>%
  inner_join(df_dms, by="name")
df_eval[df_eval$coarse_selection_type=="OrganismalFitness",]$coarse_selection_type <- "Organismal Fitness"
(df_eval_agg <- df_eval %>%
    group_by(coarse_selection_type) %>%
    summarize(spearman=mean(spearman)) %>%
    mutate(spearman_rounded=round(spearman, 2)))
ggplot(df_eval, aes(x=DMS_number_single_mutants, y=spearman,
                    color=coarse_selection_type)) +
  geom_point(alpha=0.8, shape=21, stroke=0.6, size=1) +
  geom_hline(data=df_eval_agg,
             aes(yintercept=spearman, color=coarse_selection_type),
             linetype="dotted", show.legend=F) +
  labs(x="DMS size", y="Spearman correlation") +
  theme_classic() +
  theme(aspect.ratio=1,
        legend.position="top", legend.title=element_blank(),
        axis.text=element_text(color="black")) +
  guides(color=guide_legend(nrow=2, byrow=T))
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
print(df_pred %>%
        group_by(name) %>%
        summarize(cor=cor(MIF, DMS)) %>%
        summarize(mean_cor=mean(cor)))
df_pred <- df_pred %>%
  select(-MIF) %>%
  rename(Experimental=DMS,
         `Surf. acc.`=surface_accessibility,
         Exchangeability=`AA-only`,
         `Exchangeability\n+ surf. acc.`=`Full model`)

cors_df <- c()
score_cols <- c("Experimental",
                "Exchangeability\n+ surf. acc.", "Surf. acc.", "Exchangeability",
                "Yampolsky-Stoltzfus", "BLOSUM62",
                "GEMME", "TranceptEVE", "ProSST")
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
  scale_y_discrete(limits=rev, expand=c(0,0)) +
  scale_fill_continuous(limits=c(0, 1), low="white", high="darkgreen") +
  geom_vline(xintercept=0.5, linewidth=1) +
  geom_vline(xintercept=1.505, linewidth=1) +
  theme_minimal() +
  coord_fixed(clip="off") +
  theme(aspect.ratio=1,
        panel.grid.major=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_text(color="black", angle=55, hjust=1),
        axis.text.y=element_text(color="black"))
ggsave("figures/rank_correlations.pdf", width=5, height=5)
