library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

format_labels <- function(label) paste0("underline(", label, ")")

aa_levels <- read_csv("data/aa_properties.csv") %>%
  arrange(hydropathy) %>%
  pull(aa)

# Empirical
df <- read_csv("data/data.csv") %>%
  mutate(wildtype=factor(wildtype, levels=aa_levels),
         mutation=factor(mutation, levels=aa_levels))
p <- ggplot(df) +
  facet_wrap(~wildtype, nrow=2,
             labeller=as_labeller(format_labels, default=label_parsed)) +
  geom_histogram(data=(df %>%
                         group_by(name, position, wildtype) %>%
                         summarize(surface_accessibility=median(surface_accessibility))),
                 aes(x=surface_accessibility), bins=25, fill="wheat3", alpha=0.7)
n_max <- max(ggplot_build(p)$data[[1]]$count)
p +
  geom_smooth(aes(x=surface_accessibility, y=deleteriousness*n_max,
                  color=mutation, linetype=mutation),
              linewidth=0.5, se=F) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(sec.axis=sec_axis(~./n_max, name="Rel. deleteriousness",
                                       breaks=c(0, 1/4, 1/2, 3/4, 1),
                                       labels=c("0", "1/4", "1/2", "3/4", "1")),
                     expand=c(0, 0)) +
  scale_linetype_manual(values=rep(c("solid", "dashed", "dotted", "dotdash", "longdash"),
                                   length.out=20)) +
  labs(x="Relative accessible surface area", y="Number of sites",
       color="Mutation", linetype="Mutation") +
  ggtitle(expression(underline("Wildtype"))) +
  theme_bw() +
  theme(legend.position="bottom",
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(-10,-10,-10,-10),
        panel.spacing.x=unit(0.1, "lines"), panel.spacing.y=unit(0, "lines"),
        strip.background=element_blank(), plot.title=element_text(size=11, hjust=0.5),
        axis.ticks.x=element_blank(), axis.text.x=element_blank(),
        axis.title.y.left=element_text(color="wheat3"),
        axis.text.y=element_text(color="black"),
        panel.grid.minor=element_blank()) +
  guides(color=guide_legend(nrow=2), linetype=guide_legend(nrow=2))
ggsave("figures/surface_accessibility_empirical.pdf", width=7, height=3)


# Inferred
axis_df <- c()
for (mu_bin in 1:10) {
  for (surface_accessibility in seq(0, 1, length.out=10^2))
    axis_df <- rbind(axis_df, c(mu_bin, surface_accessibility))
}
axis_df <- tibble(data.frame(axis_df))
names(axis_df) <- c("mu_bin", "surface_accessibility")
full_params <- read_csv("results/full.csv") %>%
  mutate(mu_fill=seq(min(mu), max(mu), length.out=380)) %>%
  mutate(mu_bin=cut(mu_fill, 10)) %>%
  extract(mu_bin, c("mu_bin_left", "mu_bin_right"),
          "\\D(.*),(.*)\\D", convert=T, remove=F) %>%
  mutate(mu_bin_center=(mu_bin_right-mu_bin_left)/2 + mu_bin_left,
         mu_bin=as.integer(mu_bin)) %>%
  group_by(mu_bin) %>%
  summarize(mu_bin_center=mu_bin_center[1],
            sa=sa[1])
ggplot(axis_df %>%
         inner_join(full_params, by=c("mu_bin")) %>%
         mutate(score=1/(1+exp(-(mu_bin_center + sa*surface_accessibility)))),
       aes(x=surface_accessibility, y=score, group=mu_bin)) +
  geom_line(aes(group=mu_bin), color="black", linewidth=1.5) +
  geom_line(aes(color=mu_bin), linewidth=0.75) +
  geom_label(data=full_params %>% mutate(x=(0.4-0.5*(mu_bin-5.5)/10)),
             aes(label=mu_bin, x=x, y=1/(1+exp(-(mu_bin_center + sa*x)))),
             size=2.75, label.size=1/3, label.r=unit(0, "pt")) +
  scale_x_continuous(breaks=c(0, 1/4, 1/2, 3/4, 1), labels=c("0", "1/4", "1/2", "3/4", "1"),
                     expand=c(0, 0)) +
  scale_y_continuous(breaks=c(0.05, 0.15, 0.25), expand=c(0, 0)) +
  scale_color_gradient(low="white", high="red", limits=c(1, 10)) +
  labs(x="Relative accessible surface area", y="Deleteriousness score") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        axis.text=element_text(color="black")) +
  guides(color="none")
ggsave("figures/surface_accessibility_inferred.pdf", width=4, height=4)
