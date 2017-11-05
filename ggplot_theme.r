 
my_gg_theme = theme(panel.background = element_rect(fill='white', colour='black'), 
	panel.grid.major=element_blank(), 
	panel.grid.minor= element_blank(), 
	text=element_text(family="sans"), 
	axis.text=element_text(size=15, color="black"),
	axis.ticks = element_line(color = "black"), 
	plot.title = element_text(hjust=0, size=20), 
	axis.title = element_text(size=17), 
	legend.title = element_blank(), 
	legend.text = element_text(size = 19), 
	strip.text = element_text(size = 15),
	axis.title.x = element_text(margin = margin(t= 10)),
	axis.title.y = element_text(margin = margin(r=10))	)