library(readr)
library(dplyr)
library(ggplot2)
library(shadowtext)


aas <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
aa_map <- setNames(1:length(aas), aas)
aa_levels <- read_csv("data/aa_properties.csv") %>%
  arrange(hydropathy) %>%
  pull(aa)
plot_matrix <- function(df_scores, fill_limits=NULL) {
  df_wildtype_scores <- df_scores %>%
    group_by(wildtype) %>%
    summarize(score=mean(score)) %>%
    mutate(y="0")
  df_mutation_scores <- df_scores %>%
    group_by(mutation) %>%
    summarize(score=mean(score)) %>%
    mutate(x="Z")
  shadow_size <- 0.05
  tile_textsize <- 3
  ggplot(df_scores %>%
           mutate(wildtype=factor(wildtype, levels=aa_levels),
                  mutation=factor(mutation, levels=aa_levels))) +
    geom_tile(aes(fill=score, x=wildtype, y=mutation)) +
    geom_tile(data=df_wildtype_scores %>% mutate(score=0, y="01"),
              aes(x=wildtype, y=y), fill="white", color="white") +
    geom_tile(data=df_mutation_scores %>% mutate(score=0, x="YZ"),
              aes(x=x, y=mutation), fill="white", color="white") +
    geom_tile(data=df_wildtype_scores,
              aes(x=wildtype, fill=score, y=y), color="white") +
    geom_tile(data=df_mutation_scores,
              aes(x=x, fill=score, y=mutation), color="white") +
    geom_shadowtext(aes(x=wildtype, y=mutation, label=round(score*10)),
                    bg.r=shadow_size, size=tile_textsize) +
    geom_shadowtext(data=df_wildtype_scores,
                    aes(x=wildtype, y=y, label=round(score*10)),
                    bg.r=shadow_size, size=tile_textsize) +
    geom_shadowtext(data=df_mutation_scores,
                    aes(x=x, y=mutation, label=round(score*10)),
                    bg.r=shadow_size, size=tile_textsize) +
    scale_x_discrete(breaks=names(aa_map)) +
    scale_y_discrete(limits=rev, breaks=names(aa_map)) +
    scale_fill_gradient(low="white", high="red", limits=fill_limits) +
    coord_fixed(ratio=1) +
    guides(fill="none") +
    labs(x="Wildtype", y="Mutation") +
    theme_minimal() +
    theme(panel.grid.major=element_blank(),
          axis.text=element_text(color="black"))
}

# Exchangeability model
plot_matrix(read_csv("results/aa.csv") %>%
              mutate(score=1/(1+exp(-score))),
            fill_limits=c(0, 1))
ggsave("figures/aa.pdf", width=4, height=4)

# Exchangeability and surface accessibility model
plot_matrix(read_csv("results/full.csv") %>%
              mutate(score=as.integer(cut(mu, 10))/10),
            fill_limits=c(0.1, 1))
ggsave("figures/aa_full.pdf", width=4, height=4)

# Compute correlation between inferred exchangeabilities
print(read_csv("results/full.csv") %>%
  inner_join(read_csv("results/aa.csv"), by=c("wildtype", "mutation")) %>%
  summarize(cor=cor(mu, score, method="spearman")))
