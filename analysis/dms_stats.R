library(readr)
library(dplyr)
library(ggplot2)


df_dms <- read_csv("data/dms.csv") %>%
  select(DMS_id, UniProt_ID, source_organism, taxon,
         DMS_number_single_mutants,
         selection_assay, selection_type, coarse_selection_type) %>%
  mutate(coarse_selection_type=replace(coarse_selection_type,
                                       coarse_selection_type=="OrganismalFitness",
                                       "Organismal Fitness"))

ggplot(df_dms %>%
         mutate(x=rank(-DMS_number_single_mutants)),
       aes(x=x, y=DMS_number_single_mutants, color=coarse_selection_type)) +
  geom_point(size=0.6, alpha=0.7) +
  labs(x="Deep mutational scan", y="Number of single substitutions",
       color="Selection type") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.position="top",
        axis.text=element_text(color="black")) +
  guides(color=guide_legend(nrow=2))
ggsave("figures/dms_sizes.pdf", width=6, height=4)

ggplot(read_csv("data/data.csv") %>%
         group_by(name, position) %>%
         summarize(n_mutants=n()) %>%
         group_by(n_mutants) %>%
         summarize(n=n()) %>%
         mutate(p=n/sum(n)),
       aes(x=n_mutants, y=p)) +
  geom_col() +
  scale_x_continuous(breaks=1:19) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, length.out=5),
                     labels=c("0", "1/4", "1/2", "3/4", "1"),
                     expand=c(0,0)) +
  labs(x="Number of mutations", y="Fraction of sites") +
  theme_classic() +
  theme(axis.text=element_text(color="black"))
ggsave("figures/site_sizes.pdf", width=4, height=3)

aa_levels <- read_csv("data/aa_properties.csv") %>%
  arrange(hydropathy) %>%
  pull(aa)
ggplot(read_csv("data/data.csv") %>%
         group_by(name, position, wildtype) %>%
         summarize() %>%
         group_by(wildtype) %>%
         summarize(n=n()) %>%
         mutate(wildtype=factor(wildtype, levels=aa_levels)),
       aes(wildtype, n)) +
  geom_col() +
  scale_y_continuous(expand=c(0, 0)) +
  labs(x="Wildtype", y="Number of sites") +
  theme_classic() +
  theme(axis.ticks.x=element_blank(), axis.text=element_text(color="black"))
ggsave("figures/aa_sizes.pdf", width=4, height=3)
