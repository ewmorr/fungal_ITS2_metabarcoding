require(ggplot2)
require(gridExtra)

four54_len_dist = read.table("repo/fungal_ITS2_metabarcoding/454_post_ampliconNoise_len_distribution.txt", header = T)

miseq_len_dist = read.table("repo/fungal_ITS2_metabarcoding/Miseq_merge_len_distribution.txt", header = T)
  
source(ggplot_theme.r)


p1 = ggplot(miseq_len_dist, aes(Length,Count)) +
geom_bar(stat = "identity") +
my_gg_theme +
labs(title = "a", x = "Length (bp)") +
theme(plot.title = element_text(margin = margin(t = 10, b = -20), 
	hjust = 0.01)
)


p2 = ggplot(four54_len_dist, aes(Length,Count)) +
geom_bar(stat = "identity") +
my_gg_theme +
labs(title = "b", x = "Length (bp)") +
theme(plot.title = element_text(margin = margin(t = 10, b = -20), 
	hjust = 0.01)
)


pdf("seq_length_distribution.pdf", width = 12, height = 6)
grid.arrange(p1, p2, nrow = 1)
dev.off()
