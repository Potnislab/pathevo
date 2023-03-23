##  MUTATION RATES      

    A <- ggplot(theta2, aes(x = Cultivar, y = tps_mean_depth, fill=factor(Time, level=c('Mid-season','End-season')))) +  
    geom_boxplot() + #geom_point() +
    labs( x = "", y = "Mean Mutation Rate") + ylim(0, 0.0016) +
    scale_fill_manual(name = "Time of Sampling", values=c( "#4169E1", "#FFA500"))+
    theme_bw() + 
    facet_wrap(~Environment) +
    theme(axis.text = element_text(face="bold", size=12), 
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        title =element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=12),
        legend.text=element_text(face = "bold", size=12),
        strip.text.x = element_text(face="bold", size=12)) +
        #ggtitle("") +
        theme(legend.position = "none") 


    B <- ggplot(theta1, aes(x = Cultivar, y = tps_mean_depth, fill=factor(Time, level=c('Mid-season','End-season')))) + 
    geom_boxplot() + #geom_point() +
    labs( x = "", y = "Mean Mutation Rate") + ylim(0, 0.0016) +
    #scale_fill_hue(name = "Time of Sampling") + 
    scale_fill_manual(name = "Time of Sampling", values=c( "#4169E1", "#FFA500"))+
    theme_bw() + 
    facet_wrap(~Environment) +
    #scale_x_discrete(labels=c("Ambient"="Ambient", "Elevated O3"= "Elevated O_3"))+
    theme(axis.text = element_text(face="bold", size=12), 
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        title =element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=12),
        legend.text=element_text(face = "bold", size=12),
        strip.text.x = element_text(face="bold", size=12)) +
    #ggtitle("") +
    theme(legend.position = "bottom"  ) 

    cowplot::plot_grid(A, 
                   B + 
                     theme(axis.text.y = element_blank(),
                           axis.line.y = element_blank(),
                           axis.title.y= element_blank(),
                           axis.ticks.y= element_blank()),
                   nrow = 1,
                   rel_widths = c(2.2, 2),
                   align = 'h', axis = 'tb')

    ggsave("/Users/amanpreetkaur/Downloads/theta.pdf",width = 7, height = 7, units="in", dpi=700)

