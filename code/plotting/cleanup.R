cleanup=theme_bw() +
  theme(plot.title = element_text(hjust=.5),
        panel.grid.major = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        # axis.line = element_line(color = "black"),
        # axis.text = element_text(face="bold",color="black",size="11"),
        axis.text = element_text(color="black"),
        text=element_text(size=11,face="bold"))
# Cleanup is a plotting script that I use to remove background and clean up the ggplot theme etc
