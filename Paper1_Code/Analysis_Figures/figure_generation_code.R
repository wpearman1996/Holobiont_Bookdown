#List of figures:
  fig2_envcons
setwd("..//")
########### Figure 2 - environmental conditions #########
library(ggplot2)
library(dplyr)
library(patchwork)
library(readr)
library(ggpubr)
library(tidyverse)
#The following two files are produced up in ../Sims_Main/generating_dataframe.R. We would just modify the code to adjusts the names and the files which are being utilized. For example, general_tm_0_10_data is a scenario where the host genetic contribution to fitness is 1, or microbial contribution is 0. So we would select the files for that.
general_tm_0_10_data<-read_rds("../holobiont_paper1/draft_ms/general_tm_0_10_data_means.RDS") 
general_tm_1_10_data<-read_rds("../holobiont_paper1/draft_ms/general_tm_1_10_data_means.RDS") 
general_tm_0_10_data$G1<-"1"
general_tm_1_10_data$G1<-"0"
general_tm_1_10_data<-rbind(general_tm_0_10_data,general_tm_1_10_data)
general_tm_1_10_data_means<-general_tm_1_10_data
# Split the data by EnvType into a list of data frames
rep2_data<-read_rds("../holobiont_paper1/draft_ms/rerun_1nov_largererun_data_analysed_replicate_2.RDS")
data_list <- rep2_data %>%
  filter(NGens == "1gen", VertCon == 0) %>%
  split(.$EnvType)

# Define the facet labels
facet_labels <- c(
  "highac" = "1a",
  "lowac" = "1b",
  "highac_incmean" = "2a",
  "lowac_incmean" = "2b",
  "highac_incvar" = "3a",
  "lowac_incvar" = "3b"
)

# Create individual plots in a list
plots <- lapply(names(data_list), function(name) {
  data <- data_list[[name]]
  if(name == "highac"){
    p <- ggplot(data, aes(x = Gen, y = EnvCond)) +
      geom_line(linewidth = 0.2) +
      ggtitle(facet_labels[name]) + ylim(-2,2)+
      ylab("Environmental\ncondition    ") +
      theme_classic2(base_size = 15) + scale_x_continuous(breaks = c(0, 750, 1500)) +
      theme(
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(),
        axis.ticks.y = element_line(),
    plot.margin = unit(c(0, 0.3, 0.3, 0), "cm"))
  } else{
  p <- ggplot(data, aes(x = Gen, y = EnvCond)) +
    geom_line(linewidth = 0.2) +
    ggtitle(facet_labels[name]) + ylim(-2,2)+
    theme_classic2(base_size = 15) + scale_x_continuous(breaks = c(0, 750, 1500)) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_line(),
      plot.margin = unit(c(0, 0.3, 0.3, 0), "cm"))
  }
  return(p)
})

fig2_envcons<-plots[[1]] | plots[[4]] | plots[[2]] | plots[[5]] | plots[[3]]| plots[[6]]
ggsave("../holobiont_paper1/draft_ms/fig2_envcons_4nov.jpg",fig2_envcons,dpi=600,width=12,height=3)


legend <- ggplot(general_tm_1_10_data_means) +
  aes(x = Gen, y = HostMicrobeFitness.mean, colour = VertCon) +
  geom_point(shape = "circle", size = 0.5) +
  scale_color_hue(direction = 1, name = "Proportion of vertical inheritance") + 
  theme_classic() + 
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(colour = guide_legend(override.aes = list(size = 5), nrow = 1))

shared_legend <- cowplot::get_plot_component(legend, "guide-box", return_all = TRUE)[[3]]

shared_legend<-cowplot::ggdraw(shared_legend)
########### Figure 3 - Fitness - NGen = 1 ###########
fig_envs <- lapply(names(data_list), function(name) {
  data <- data_list[[name]]
  if(name == "highac"){
    p <- ggplot(data, aes(x = Gen, y = EnvCond)) +
      geom_line(linewidth = 0.2) +
      ylim(-2,2)+
      ylab("Environmental\ncondition    ") +
      theme_classic2(base_size = 15) + scale_x_continuous(breaks = c(0, 750, 1500)) +
      theme(
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(),
        axis.ticks.y = element_line(),
        plot.margin = unit(c(0, 0.3, 0.3, 0), "cm"))
  } else{
    p <- ggplot(data, aes(x = Gen, y = EnvCond)) +
      geom_line(linewidth = 0.2) +
      ylim(-2,2)+
      theme_classic2(base_size = 15) + scale_x_continuous(breaks = c(0, 750, 1500)) +
      theme(
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_line(),
        plot.margin = unit(c(0, 0.3, 0.3, 0), "cm"))
  }
  return(p)
})
fig_envs<-fig_envs[[1]] |  fig_envs[[4]]| fig_envs[[2]]  | fig_envs[[5]]| fig_envs[[3]] | fig_envs[[6]]

theme_classic_WP<-theme_classic2()


multigen_data_means_g1<-general_tm_1_10_data_means
#multigen_data_means_g1$G1<-"0"
#g1_data_means$G1<-"1"
#multigen_data_means_g1<-rbind(multigen_data_means_g1,g1_data_means)

esquisse::esquisser(multigen_data_means_g1 %>%
  filter(NGens %in% "1gen"))

host_fit_lowngen<- multigen_data_means_g1 %>%
  filter(NGens %in% "1gen") %>% filter(G1=="0") %>%
  mutate(EnvType = factor(EnvType, levels = c("highac","lowac", "highac_incmean",  "lowac_incmean","highac_incvar",  "lowac_incvar"))) %>%  # Set facet order
  ggplot() +
  aes(x = Gen, y = HostFitness.mean, colour = VertCon) +
  geom_point(aes(colour = VertCon), shape = "circle", size = 0.2, alpha = 0.1) +  # Plot regular points with color
  geom_point(data = subset(multigen_data_means_g1, G1 == "1"), colour = "black", shape = "circle", size = 0.2, alpha = 0.01, show.legend = FALSE) +  # Add black dots without legend
  geom_smooth(aes(colour = VertCon), span = 0.75) +  # Plot regular smoothed lines
  geom_smooth(data = subset(multigen_data_means_g1, G1 == "1"), colour = "black", span = 0.75, show.legend = FALSE) +  # Add black lines for G1 without legend
  scale_color_hue(direction = 1) +
  theme_classic2(base_size=15) +
  facet_wrap(vars(EnvType), ncol = 6,
             labeller = labeller(EnvType = c(
               "highac" = "g",
               "lowac" = "h",
               "highac_incmean" = "h",
               "lowac_incmean" = "i",
               "highac_incvar" = "j",
               "lowac_incvar" = "k"
             ))) +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  )  +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,linewidth=1.8) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,linewidth=1.8) +
  scale_x_continuous(breaks = c(0, 750, 1500)) + ylim(0.25,1) +
  theme(legend.position = "bottom", legend.direction = "horizontal") + ylab("Host Fitness") +
  xlab("Generation") 


microbe_fit_lowngen <- multigen_data_means_g1 %>%
  filter(NGens %in% "1gen") %>%  filter(G1=="0") %>%
  mutate(EnvType = factor(EnvType, levels = c("highac","lowac", "highac_incmean",  "lowac_incmean","highac_incvar",  "lowac_incvar"))) %>%  # Set facet order
  ggplot() +
  aes(x = Gen, y = HostMicrobeFitness.mean, colour = VertCon) +
  geom_point(aes(colour = VertCon), shape = "circle", size = 0.2, alpha = 0.1) +  # Plot regular points with color
  #geom_point(data = subset(multigen_data_means_g1, G1 == "1"), colour = "black", shape = "circle", size = 0.2, alpha = 0.01, show.legend = FALSE) +  # Add black dots without legend
  geom_smooth(aes(colour = VertCon), span = 0.75) +  # Plot regular smoothed lines
  #geom_smooth(data = subset(multigen_data_means_g1, G1 == "1"), colour = "black", span = 0.75, show.legend = FALSE) +  # Add black lines for G1 without legend
  scale_color_hue(direction = 1) +
  theme_classic2(base_size=15) +
  facet_wrap(vars(EnvType), ncol = 6,
             labeller = labeller(EnvType = c(
               "highac" = "a",
               "lowac" = "b",
               "highac_incmean" = "c",
               "lowac_incmean" = "d",
               "highac_incvar" = "e",
               "lowac_incvar" = "f"
             ))) +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  ) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,linewidth=1.8) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,linewidth=1.8) +
  scale_x_continuous(breaks = c(0, 750, 1500)) + ylim(0.6,1)+
  theme(legend.position = "bottom", legend.direction = "horizontal")+ ylab("Microbe Fitness") +
  xlab("Generation")


fitplots_lowngen<-cowplot::plot_grid(fig_envs,microbe_fit_lowngen + theme(legend.position="none") ,host_fit_lowngen + theme(legend.position="none"),
                                     shared_legend,
                                     nrow = 4,rel_heights = c(2,3,3,1),align = "v",axis="tblr" )

fitplots_lowngen



ggsave("../holobiont_paper1//draft_ms/fig3_lowngen_Fitness_4nov.jpg",fitplots_lowngen,dpi=600,width=12,height=7)





########### Figure 4 - Fitness - NGen = 200 #########


host_fit_highngen<- general_tm_1_10_data_means %>%
  filter(NGens %in% "200gen") %>% filter(G1=="0") %>% 
  mutate(EnvType = factor(EnvType, levels = c("highac","lowac", "highac_incmean",  "lowac_incmean","highac_incvar",  "lowac_incvar"))) %>%  # Set facet order
  ggplot() +
  aes(x = Gen, y = HostFitness.mean, colour = VertCon) +
  geom_point(aes(colour = VertCon), shape = "circle", size = 0.2, alpha = 0.1) +  # Plot regular points with color
  geom_point(data = subset(multigen_data_means_g1, G1 == "1"), colour = "black", shape = "circle", size = 0.2, alpha = 0.01, show.legend = FALSE) +  # Add black dots without legend
  geom_smooth(aes(colour = VertCon), span = 0.75) +  # Plot regular smoothed lines
  geom_smooth(data = subset(multigen_data_means_g1, G1 == "1"), colour = "black", span = 0.75, show.legend = FALSE) +  # Add black lines for G1 without legend
  scale_color_hue(direction = 1) +
  theme_classic2(base_size=15) +
  facet_wrap(vars(EnvType), ncol = 6,
             labeller = labeller(EnvType = c(
               "highac" = "g",
               "lowac" = "h",
               "highac_incmean" = "h",
               "lowac_incmean" = "i",
               "highac_incvar" = "j",
               "lowac_incvar" = "k"
             ))) +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  )  +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,linewidth=1.8) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,linewidth=1.8) +
  scale_x_continuous(breaks = c(0, 750, 1500)) + ylim(0.25,1) +
  theme(legend.position = "bottom", legend.direction = "horizontal") + ylab("Host Fitness") +
  xlab("Generation") 


microbe_fit_highngen <- general_tm_1_10_data_means %>%
  filter(NGens %in% "200gen") %>% filter(G1=="0") %>% 
  mutate(EnvType = factor(EnvType, levels = c("highac","lowac", "highac_incmean",  "lowac_incmean","highac_incvar",  "lowac_incvar"))) %>%  # Set facet order
  ggplot() +
  aes(x = Gen, y = HostMicrobeFitness.mean, colour = VertCon) +
  geom_point(shape = "circle", size = 0.2,alpha=0.2) +
  geom_smooth(span = 0.75) +
  scale_color_hue(direction = 1) +
  theme_classic2(base_size=15) +
  facet_wrap(vars(EnvType), ncol = 6,
             labeller = labeller(EnvType = c(
               "highac" = "a",
               "lowac" = "b",
               "highac_incmean" = "c",
               "lowac_incmean" = "d",
               "highac_incvar" = "e",
               "lowac_incvar" = "f"
             ))) +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  ) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,linewidth=1.8) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,linewidth=1.8) +
  scale_x_continuous(breaks = c(0, 750, 1500)) + ylim(0.6,1)+
  theme(legend.position = "bottom", legend.direction = "horizontal")+ ylab("Microbe Fitness") +
  xlab("Generation")


fitplots_highngen<-cowplot::plot_grid(fig_envs,microbe_fit_highngen + theme(legend.position="none") ,host_fit_highngen + theme(legend.position="none"),
                                      shared_legend,
                                     nrow = 4,rel_heights = c(2,3,3,1),align = "v",axis="tblr" )
fitplots_highngen



ggsave("../holobiont_paper1/draft_ms/fig3_highngen_Fitness_4nov.jpg",fitplots_highngen,dpi=600,width=12,height=7)
 

########### Figure 5 - Diversity Heatmaps
adiv <- general_tm_1_10_data_means %>%
  filter(Gen == 1500) %>%filter(G1=="0") %>% 
  mutate(VertCon_numeric = as.numeric(as.character(VertCon)))  %>%
  mutate(VertCon_rounded = round(VertCon_numeric, 3))  %>%
  mutate(VertCon_rounded_factor = (as.factor(VertCon_rounded)))  %>%
  mutate(EnvType = factor(EnvType, levels = c("highac","lowac", "highac_incmean",  "lowac_incmean","highac_incvar",  "lowac_incvar"))) %>%  # Set facet order
  ggplot() +
  aes(x = GenLab, y = VertCon_rounded_factor, fill = diversity.mean) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  theme_minimal(base_size = 15) +
  facet_wrap(vars(EnvType), ncol = 6,
             labeller = labeller(EnvType = c(
               "highac" = "a",
               "lowac" = "b",
               "highac_incmean" = "c",
               "lowac_incmean" = "d",
               "highac_incvar" = "e",
               "lowac_incvar" = "f"
             ))) +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  )  +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,linewidth=1.8) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,linewidth=1.8) +
  ylab("Proportion of vertically\nacquired microbes") +
  xlab("Number of microbial generations per host generation") +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  ) +
  labs(fill = expression(atop("Alpha Diversity"))) 

library(ggplot2)

microbe_div_highngen <- general_tm_1_10_data_means %>% filter(G1=="0") %>% 
#  mutate(diversity.mean = ifelse(is.na(diversity.mean) | is.nan(diversity.mean), 0, diversity.mean)) %>%
  filter(NGens %in% "200gen") %>%
  mutate(EnvType = factor(EnvType, levels = c("highac","lowac", "highac_incmean",  "lowac_incmean","highac_incvar",  "lowac_incvar"))) %>%  # Set facet order
  ggplot() +
  aes(x = Gen, y = diversity.mean, colour = VertCon) +
  geom_point(shape = "circle", size = 0.2, alpha = 0.2) +
  geom_smooth(span = 0.75) +
  scale_color_hue(direction = 1) +
  theme_classic2(base_size = 15) +
  facet_wrap(vars(EnvType), ncol = 6,
             labeller = labeller(EnvType = c(
               "highac" = "g",
               "lowac" = "h",
               "highac_incmean" = "h",
               "lowac_incmean" = "i",
               "highac_incvar" = "j",
               "lowac_incvar" = "k"
             ))) +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  )  +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,linewidth=1.8) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,linewidth=1.8) +
  scale_x_continuous(breaks = c(0, 750, 1500)) +  ylim(0,1) +
  theme(
    legend.position = "none",  # Remove legend
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"),  # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis tick labels
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  ) +
  ylab(expression(atop("Alpha diversity", paste(T["M"]==200, "\n \n", T["H"]==1500)))) +  # Updated y-axis label with line breaks and subscript notation
  xlab("Generation") 

microbe_div_lowngen <- general_tm_1_10_data_means %>%
  #mutate(diversity.mean = ifelse(is.na(diversity.mean) | is.nan(diversity.mean), 0, diversity.mean)) %>%
  filter(NGens %in% "1gen") %>%filter(G1=="0") %>% 
  mutate(EnvType = factor(EnvType, levels = c("highac","lowac", "highac_incmean",  "lowac_incmean","highac_incvar",  "lowac_incvar"))) %>%  # Set facet order
  ggplot() +
  aes(x = Gen, y = diversity.mean, colour = VertCon) +
  geom_point(shape = "circle", size = 0.2, alpha = 0.2) +
  geom_smooth(span = 0.75) +
  scale_color_hue(direction = 1) +
  theme_classic2(base_size = 15) +
  facet_wrap(vars(EnvType), ncol = 6,
             labeller = labeller(EnvType = c(
               "highac" = "a",
               "lowac" = "b",
               "highac_incmean" = "c",
               "lowac_incmean" = "d",
               "highac_incvar" = "e",
               "lowac_incvar" = "f"
             ))) +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  )  +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,linewidth=1.8) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,linewidth=1.8) +
  scale_x_continuous(breaks = c(0, 750, 1500)) + ylim(0,1) +
  ylab(expression(atop("Alpha diversity", paste(T["M"]==1, "\n \n", T["H"]==1500)))) +  # Updated y-axis label with line breaks and subscript notation
  xlab("Generation") 

# Print the plot
print(microbe_div_lowngen)


div_heatplots<-cowplot::plot_grid(fig_envs,adiv,
                                  nrow = 2,rel_heights = c(3,3),align = "v",axis = "tblr")


div_timeplots<-cowplot::plot_grid(fig_envs,microbe_div_lowngen + theme(legend.position = "none"),microbe_div_highngen,
                                  nrow = 3,rel_heights = c(2,3,3),align = "v",axis = "tblr")


div_timeplots<-cowplot::plot_grid(div_timeplots,shared_legend,nrow=2,rel_heights = c(6,1))
div_timeplots
div_heatplots
ggsave("../holobiont_paper1//draft_ms/fig4_time_diversity_4nov.jpg",div_timeplots,dpi=600,width=12,height=7)
ggsave("../holobiont_paper1//draft_ms/fig4_heatmap_diversity_4nov.jpg",div_heatplots,dpi=600,width=18,height=5)


########### Figure 5 Fitness heat plots ##############
host_fitness_range1 <- range(general_tm_1_10_data_means$HostFitness.mean[general_tm_1_10_data_means$Gen==1500 & general_tm_1_10_data_means$G1=="0"], na.rm = TRUE)
microbe_fitness_range1 <- range(general_tm_1_10_data_means$HostMicrobeFitness.mean[general_tm_1_10_data_means$Gen==1500 & general_tm_1_10_data_means$G1=="0"], na.rm = TRUE)
host_fitness_range2 <- range(effect_VI_data_means$`Host Fitness`[effect_VI_data_means$Gen==1500 ], na.rm = TRUE)
microbe_fitness_range2 <- range(effect_VI_data_means$`Microbe Fitness`[effect_VI_data_means$Gen==1500], na.rm = TRUE)

host_fitness_range<-range(c(host_fitness_range1,host_fitness_range2))
microbe_fitness_range<-range(c(microbe_fitness_range1,microbe_fitness_range2))

microbe_fitness <- general_tm_1_10_data_means %>%
  mutate(EnvType = factor(EnvType, levels = c("highac","lowac", "highac_incmean",  "lowac_incmean","highac_incvar",  "lowac_incvar"))) %>%  # Set facet order
  filter(Gen == 1500) %>%
  mutate(VertCon_numeric = as.numeric(as.character(VertCon)))  %>%
  mutate(VertCon_rounded = round(VertCon_numeric, 3))  %>%
  mutate(VertCon_rounded_factor = (as.factor(VertCon_rounded)))  %>%
  ggplot() +
  aes(x = GenLab, y = VertCon_rounded_factor, fill = HostMicrobeFitness.mean) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = microbe_fitness_range) +
  theme_minimal(base_size = 15) +
  facet_wrap(vars(EnvType), nrow = 1, strip.position = "top",
             labeller = labeller(EnvType = c(
               "highac" = "a",
               "lowac" = "b",
               "highac_incmean" = "c",
               "lowac_incmean" = "d",
               "highac_incvar" = "e",
               "lowac_incvar" = "f"
             ))) +
  ylab("Proportion of vertically\nacquired microbes") +
  xlab("Number of microbial generations per host generation") +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  ) +
  labs(fill = expression(atop("Microbe Fitness"))) 


host_fitness <- general_tm_1_10_data_means %>%
  mutate(EnvType = factor(EnvType, levels = c("highac","lowac", "highac_incmean",  "lowac_incmean","highac_incvar",  "lowac_incvar"))) %>%  # Set facet order
  filter(Gen == 1500) %>%
  mutate(VertCon_numeric = as.numeric(as.character(VertCon)))  %>%
  mutate(VertCon_rounded = round(VertCon_numeric, 3))  %>%
  mutate(VertCon_rounded_factor = (as.factor(VertCon_rounded)))  %>%
  ggplot() +
  aes(x = GenLab, y = VertCon_rounded_factor, fill = HostFitness.mean) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = host_fitness_range) +
  theme_minimal(base_size = 15) +
  facet_wrap(vars(EnvType), nrow = 1, strip.position = "top",
             labeller = labeller(EnvType = c(
               "highac" = "g",
               "lowac" = "h",
               "highac_incmean" = "h",
               "lowac_incmean" = "i",
               "highac_incvar" = "j",
               "lowac_incvar" = "k"
             ))) +
  ylab("Proportion of vertically\nacquired microbes") +
  xlab("Number of microbial generations per host generation") +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  ) +
  labs(fill = expression(atop("Host Fitness"))) 

Fig5_fitness_heatplots<-cowplot::plot_grid(fig_envs,microbe_fitness,host_fitness,
                                  nrow = 3,rel_heights = c(2,3,3),align = "v",axis = "tblr")
Fig5_fitness_heatplots
ggsave("./draft_ms/Fig5_fitness_heatplots4nov.jpg",Fig5_fitness_heatplots,dpi=600,width=18,height=7)


########### Figure 6 Heatplot ratios 





######### Figure 7 - Effective vertical inheritance

library(dplyr)
viests_summary <- viests %>%
  group_by(Var2, InitVert) %>%
  summarize(mean_value = mean(value),
            sd_value = sd(value),
            .groups = 'drop') 
viests_summary$Facet<-"x"

efi_relationship <- ggplot(viests_summary) +
  aes(x = Var2, y = mean_value, group = InitVert, colour = InitVert) +
  geom_line(linewidth = 2) + 
  geom_ribbon(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
              alpha = 0.2, fill = "grey") + 
  scale_color_hue(direction = 1) +
  theme_classic2() + 
  xlab("Microbial Generations within a host") + 
  ylab("Effective Vertical Inheritance")  + 
  facet_wrap(vars(Facet), ncol = 1,
             labeller = labeller(Facet = c(
               "x" = "a"))) +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(), axis.line.y = element_blank(),
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    strip.background = element_blank()
  )  +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=1.8) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=1.8)


# Find the range of HostFitness.mean

efi_hostfit<-effect_VI_data_means %>%
  filter(Gen >= 1500L & Gen <= 1500L) %>%
  mutate(EnvType = factor(EnvType, levels = c("highac","lowac", "highac_incmean",  "lowac_incmean","highac_incvar",  "lowac_incvar"))) %>%  # Set facet order
  ggplot() +
  aes(
    x = Group,
    y = VertCon,
    fill = `Host Fitness`
  ) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = host_fitness_range) +  # Apply the extracted limits
#  scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  theme_minimal() +
  facet_wrap(vars(EnvType), ncol = 6,
             labeller = labeller(EnvType = c(
               "highac" = "h",
               "lowac" = "i",
               "highac_incmean" = "j",
               "lowac_incmean" = "k",
               "highac_incvar" = "l",
               "lowac_incvar" = "m"
             ))) + xlab(NULL) + 
  ylab("Vertical Inheritance") +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    # axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  )   +  scale_y_discrete(position="right")


efi_microbefit<-effect_VI_data_means %>%
  filter(Gen >= 1500L & Gen <= 1500L) %>%
  mutate(EnvType = factor(EnvType, levels = c("highac","lowac", "highac_incmean",  "lowac_incmean","highac_incvar",  "lowac_incvar"))) %>%  # Set facet order
  ggplot() +
  aes(
    x = Group,
    y = VertCon,
    fill = `Microbe Fitness`
  ) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = microbe_fitness_range) +
#  scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  theme_minimal() +
  facet_wrap(vars(EnvType), ncol = 6,
             labeller = labeller(EnvType = c(
               "highac" = "b",
               "lowac" = "e",
               "highac_incmean" = "c",
               "lowac_incmean" = "f",
               "highac_incvar" = "d",
               "lowac_incvar" = "g"
             )))+ xlab(NULL) + 
  ylab("Vertical Inheritance") +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    #axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  )   +  scale_y_discrete(position="right")
fig_envs_evi <- lapply(names(data_list), function(name) {
  data <- data_list[[name]]
  if(name == "highac"){
    p <- ggplot(data, aes(x = Gen, y = EnvCond)) +
      geom_line(linewidth = 0.2) +
      ylim(-2,2)+
      ylab("") +
      theme_classic2(base_size = 15) + scale_x_continuous(breaks = c(0, 750, 1500)) +
      theme(
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(),
        axis.ticks.y = element_line(),
        plot.margin = unit(c(0, 0.3, 0.3, 0), "cm"))
  } else{
    p <- ggplot(data, aes(x = Gen, y = EnvCond)) +
      geom_line(linewidth = 0.2) +
      ylim(-2,2)+
      theme_classic2(base_size = 15) + scale_x_continuous(breaks = c(0, 750, 1500)) +
      theme(
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_line(),
        plot.margin = unit(c(0, 0.3, 0.3, 0), "cm"))
  }
  return(p)
})
fig_envs_evi<-fig_envs_evi[[1]] |  fig_envs_evi[[4]]| fig_envs_evi[[2]]  | fig_envs_evi[[5]]| fig_envs_evi[[3]] | fig_envs_evi[[6]]

efi_heats<-cowplot::plot_grid(fig_envs_evi,efi_microbefit,efi_hostfit,nrow=3,align = "v",rel_heights = c(1,2,2),axis = "tblr")


legendefi <- ggplot(general_tm_1_10_data_means) +
  aes(x = Gen, y = HostMicrobeFitness.mean, colour = VertCon) +
  geom_point(shape = "circle", size = 0.5) +
  scale_color_hue(direction = 1, name = "Proportion of vertical inheritance") + 
  theme_classic() + 
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(colour = guide_legend(override.aes = list(size = 5), nrow = 3))

shared_legend_efi <- cowplot::get_plot_component(legendefi, "guide-box", return_all = TRUE)[[3]]


efi_relationship_legend<-cowplot::plot_grid(efi_relationship + theme(legend.position = "none"),
                                            shared_legend_efi,nrow=2,rel_heights = c(6,1),align = "v")



efi_plot<-cowplot::plot_grid(efi_relationship_legend ,efi_heats,ncol=2,rel_widths = c(2,4))
efi_plot
ggsave("./draft_ms/Fig7_effectivevertinher.jpg",efi_plot,dpi=600,width=15,height=6)





microbe_phenovar_highngen <- general_tm_1_10_data_means %>% filter(G1=="0") %>% 
  filter(NGens %in% "200gen") %>%
  mutate(EnvType = factor(EnvType, levels = c("highac","lowac", "highac_incmean",  "lowac_incmean","highac_incvar",  "lowac_incvar"))) %>%  # Set facet order
  ggplot() +
  aes(x = Gen, y = PhenoVar, colour = VertCon) +
  geom_point(shape = "circle", size = 0.2, alpha = 0.2) +
  geom_smooth(span = 0.75) +
  scale_color_hue(direction = 1) +
  theme_classic2(base_size = 15) +
  facet_wrap(vars(EnvType), ncol = 6,
             labeller = labeller(EnvType = c(
               "highac" = "g",
               "lowac" = "h",
               "highac_incmean" = "h",
               "lowac_incmean" = "i",
               "highac_incvar" = "j",
               "lowac_incvar" = "k"
             ))) +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  )  +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,linewidth=1.8) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,linewidth=1.8) +
  scale_x_continuous(breaks = c(0, 750, 1500)) +  ylim(0,.5) +
  theme(
    legend.position = "none",  # Remove legend
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"),  # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis tick labels
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  ) +
  ylab(expression(atop("Phenotypic Variance", paste(T["M"]==200, "\n \n", T["H"]==1500)))) +  # Updated y-axis label with line breaks and subscript notation
  xlab("Generation") 


microbe_phenovar_lowngen <- general_tm_1_10_data_means %>%filter(G1=="0") %>% 
  filter(NGens %in% "1gen") %>%
  mutate(EnvType = factor(EnvType, levels = c("highac","lowac", "highac_incmean",  "lowac_incmean","highac_incvar",  "lowac_incvar"))) %>%  # Set facet order
  ggplot() +
  aes(x = Gen, y = PhenoVar, colour = VertCon) +
  geom_point(shape = "circle", size = 0.2, alpha = 0.2) +
  geom_smooth(span = 0.75) +
  scale_color_hue(direction = 1) +
  theme_classic2(base_size = 15) +
  facet_wrap(vars(EnvType), ncol = 6,
             labeller = labeller(EnvType = c(
               "highac" = "a",
               "lowac" = "b",
               "highac_incmean" = "c",
               "lowac_incmean" = "d",
               "highac_incvar" = "e",
               "lowac_incvar" = "f"
             ))) +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  )  +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,linewidth=1.8) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,linewidth=1.8) +
  scale_x_continuous(breaks = c(0, 750, 1500)) + ylim(0,0.5) +
  theme(
    legend.position = "none",  # Remove legend
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"),  # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis tick labels
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  ) +
  ylab(expression(atop("Phenotypic Variance", paste(T["M"]==1, "\n \n", T["H"]==1500)))) +  # Updated y-axis label with line breaks and subscript notation
  xlab("Generation") 


phenovar <- general_tm_1_10_data_means %>%filter(G1=="0") %>% 
  filter(Gen == 1500) %>%
  mutate(VertCon_numeric = as.numeric(as.character(VertCon)))  %>%
  mutate(VertCon_rounded = round(VertCon_numeric, 3))  %>%
  mutate(VertCon_rounded_factor = (as.factor(VertCon_rounded)))  %>%
  mutate(EnvType = factor(EnvType, levels = c("highac","lowac", "highac_incmean",  "lowac_incmean","highac_incvar",  "lowac_incvar"))) %>%  # Set facet order
  ggplot() +
  aes(x = GenLab, y = VertCon_rounded_factor, fill = PhenoVar) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  theme_minimal(base_size = 15) +
  facet_wrap(vars(EnvType), ncol = 6,
             labeller = labeller(EnvType = c(
               "highac" = "a",
               "lowac" = "b",
               "highac_incmean" = "c",
               "lowac_incmean" = "d",
               "highac_incvar" = "e",
               "lowac_incvar" = "f"
             ))) +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  )  +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,linewidth=1.8) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,linewidth=1.8) +
  ylab("Proportion of vertically\nacquired microbes") +
  xlab("Number of microbial generations per host generation") +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  ) +
  labs(fill = expression(atop("Phenotypic Variance"))) 

phenovar_timeplots<-cowplot::plot_grid(fig_envs,microbe_phenovar_lowngen + theme(legend.position = "none"),microbe_phenovar_highngen,
                                       nrow = 3,rel_heights = c(2,3,3),align = "v",axis = "tblr")


phenovar_timeplots<-cowplot::plot_grid(phenovar_timeplots,shared_legend,nrow=2,rel_heights = c(6,1))
phenovar<-cowplot::plot_grid(fig_envs,phenovar ,
                                  nrow = 2,rel_heights = c(1,3),align = "v",axis = "tblr")
ggsave("./draft_ms/SuppFig1_phenovar.jpg",phenovar_timeplots,dpi=600,width=12,height=7)

ggsave("./draft_ms/SuppFig2_phenovar.jpg",phenovar,dpi=600,width=20,height=5)
