library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)


start_position <- 110
df_region <- read_csv("results/predictions.csv") %>%
  filter(name=="BLAT_ECOLX_Firnberg_2014",
         position>start_position, position<start_position+50)

df_plot <- df_region %>%
  filter(mutation=="H") %>%
  select(wildtype, position,
         deleteriousness,
         surface_accessibility, `AA-only`, `Full model`) %>%
  mutate(`Full model`=rank(-`Full model`)/n(),
         surface_accessibility=rank(-surface_accessibility)/n(),
         `AA-only`=rank(-`AA-only`)/n()) %>%
  rename(Experimental=deleteriousness,
         `Surf. acc.`=surface_accessibility,
         Exchangeability=`AA-only`,
         `Exchangeability + surf. acc.`=`Full model`) %>%
  pivot_longer(!c(wildtype, position), names_to="feature") %>%
  mutate(feature=factor(feature, levels=c("Exchangeability",
                                          "Surf. acc.",
                                          "Exchangeability + surf. acc.",
                                          "Experimental")))

df_labels <- df_region %>%
  group_by(position) %>%
  summarize(wildtype=wildtype[1])
ggplot(df_plot, aes(position, feature)) +
  geom_tile(aes(fill=value), show.legend=F) +
  geom_tile(fill=NA, color="black") +
  scale_x_continuous(breaks=df_labels$position, labels=df_labels$wildtype,
                     expand=c(0, 0)) +
  scale_fill_gradient(low="white", high="red", limits=c(0, 1)) +
  coord_fixed(ratio=1) +
  xlab("Wildtype") +
  theme_classic() +
  theme(axis.text=element_text(color="black"),
        axis.title.x=element_text(size=9),
        axis.title.y=element_blank(), axis.ticks=element_blank())
ggsave("figures/demo.pdf", width=6.5, height=1)
print(with(df_region %>% filter(mutation=="H"),
           cor(deleteriousness, `Full model`, method="spearman")))
