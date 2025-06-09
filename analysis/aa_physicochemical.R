library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)


# Read physicochemical properties of each amino acid
aa_properties <- read_csv("data/aa_properties.csv")


# Wildtype-mutant asymmetry
aa_levels <- aa_properties %>%
  arrange(hydropathy) %>%
  pull(aa)
score_asymmetries <- read_csv("results/aa.csv") %>%
  mutate(score=1/(1+exp(-score)),
         wildtype=factor(wildtype, levels=aa_levels),
         mutation=factor(mutation, levels=aa_levels)) %>%
  mutate(aa_pair=ifelse(as.integer(wildtype)<as.integer(mutation),
                        paste(wildtype, mutation),
                        paste(mutation, wildtype))) %>%
  mutate_if(is.character, factor) %>%
  group_by(aa_pair) %>%
  mutate(score_diff=(sum(score)-score)-score) %>%
  arrange(-abs(score_diff))
ggplot(score_asymmetries) +
  geom_tile(aes(wildtype, mutation, fill=-score_diff*10)) +
  scale_y_discrete(limits=rev) +
  scale_fill_distiller(type="div",
                       limits=max(abs(score_asymmetries$score_diff)) * c(-10, 10),
                       breaks=-3:3,
                       palette="RdBu", direction=-1) +
  coord_fixed(ratio=1) +
  labs(x="Wildtype", y="Mutation", fill="Score asymmetry") +
  theme_minimal() +
  theme(legend.position="bottom", panel.grid.major=element_blank(),
        axis.text.x=element_text(color="black"), axis.text.y=element_text(color="black"),
        legend.key.height=unit(0.5, "lines"), legend.key.width=unit(1.5, "lines")) +
  guides(fill=guide_colorbar(ticks.colour=NA))
ggsave("figures/aa_asymmetry.pdf", width=4, height=4)

# Wildtype analysis
aa_inferred <- read_csv("results/aa.csv") %>%
  mutate(score=1/(1+exp(-score))) %>%
  inner_join(aa_properties %>% rename(wildtype=aa), by=c("wildtype"))
print(head(aa_inferred %>% arrange(score), 3))
(aa_inferred %>%
    group_by(wildtype) %>%
    summarize(score=mean(score)) %>%
    arrange(-score) %>%
    head(5))
ggplot(aa_inferred %>%
         group_by(wildtype) %>%
         summarize(volume=mean(volume), score=mean(score)),
       aes(volume, score)) +
  geom_text(aes(label=wildtype)) +
  labs(x="Amino acid volume (cubic angstroms)",
       y="Average deleteriousness score as wildtype") +
  theme_classic() +
  theme(axis.text=element_text(color="black"))
ggsave("figures/aa_volume.pdf", width=4, height=4)
