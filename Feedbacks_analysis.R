#####--------------------------------------------------------
# Load required libraries
library(tidyverse)
library(lme4)
library(cowplot)
library(ggthemr)
library(piecewiseSEM)
library(scales)
library(ggpmisc)
library(ggpubr)

#set ggplot theme
ggthemr('fresh', text_size = 12)

# Set working directory
setwd("~/Project Seagrass Dropbox/Team Project Seagrass/Projects/Academic Papers/Feedbacks/Analysis")

#load data

seed_method<-read.csv("All_seed.csv", header = TRUE, sep=",")
names(seed_method)
summary(seed_method)
str(seed_method)

#recode variables
seed_method$Site <- as.factor(seed_method$Site)
seed_method$Treatment <- as.factor(seed_method$Treatment)
seed_method$Treatment <- factor(seed_method$Treatment, levels = c("Furrows", "Buried bags", "Surface bags"))
seed_method$Algae <- as.numeric(seed_method$Algae)
seed_method <- seed_method %>% mutate(Shoot_count_m2 = Shoot_count * 4)

#mean and SD for shoot count
seed_method_plot <- seed_method %>%
  group_by(Site, Treatment) %>% 
  summarise(Shoot_count = mean(Shoot_count_m2),
            Shoot_count_sd = sd(Shoot_count_m2),
            Leaf = mean(Max_leaf),
            Leaf_sd = sd(Max_leaf),
            Algae_m = mean(Algae),
            Algae_sd = sd(Algae))
seed_method_plot

summary(seed_method_plot)
summary(seed_method$Algae)

#plot of shoot density + shoot length + algae at plot level, with plot standard deviations

a <- ggplot(seed_method_plot, aes(x = Treatment, y = Shoot_count))  +
  geom_point(aes(x = Treatment, y = Shoot_count, group = Site), position=position_dodge(width=0.5), alpha = .8, size = 4) + 
  geom_errorbar(aes(ymin=Shoot_count-Shoot_count_sd, ymax=Shoot_count+Shoot_count_sd, group = Site), alpha = .8, width=.2,
                position=position_dodge(width=0.5)) +
  stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75), colour = "#233B43") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = .3, position = position_dodge(width = .75), colour = "#233B43") +
  labs(y = bquote('Shoot density '("per m"^2)), x = "") +
  theme(legend.position = "bottom") +
  #coord_cartesian(ylim = c(0, 12)) +
  scale_y_continuous(
    labels = label_number(accuracy = 1))
a

b <- ggplot(seed_method_plot, aes(x = Treatment, y = Leaf))  +
  geom_point(aes(group = Site), position=position_dodge(width=0.5), alpha = .8, size = 4) + 
  geom_errorbar(aes(ymin=Leaf-Leaf_sd, ymax=Leaf+Leaf_sd, group = Site), alpha = .8, width=.2,
                position=position_dodge(width=0.5)) +
  stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75), colour = "#233B43") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = .3, position = position_dodge(width = .75), colour = "#233B43") +
  labs(y = "Max leaf length (mm)", x = "") +
  theme(legend.position = "bottom")
  #coord_cartesian(ylim = c(0, 200))
b


c <- ggplot(seed_method_plot, aes(x = Treatment, y = Algae_m))  +
  geom_point(aes(group = Site), position=position_dodge(width=0.5), alpha = .8, size = 4) + 
  geom_errorbar(aes(ymin=Algae_m-Algae_sd, ymax=Algae_m+Algae_sd, group = Site), alpha = .8, width=.2,
                position=position_dodge(width=0.5)) +
  stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75), colour = "#233B43") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = .3, position = position_dodge(width = .75), colour = "#233B43") +
  labs(y = "Macroalgal cover (%)", x = "Planting treatment") +
  theme(legend.position = "bottom")
c

Figure_1 <- plot_grid(a,b,c, nrow = 3, ncol = 1, labels = c('A', 'B', 'C'), align = "v")
Figure_1

#save plot
ggsave("Figure_1.pdf", plot = Figure_1, device = "pdf",
       width = 20,
       height = 25,
       units = c("cm"),
       dpi = 300)

#plot of correlation between shoot count and algae
Figure_2 <- ggplot(seed_method_plot, aes(y = Shoot_count, x = Algae_m)) +
  geom_point(size=8, aes(colour = Treatment)) +
  #geom_smooth(method=glm, fullrange=TRUE) +
  labs(y = bquote('Shoot density '("per m"^2)), x = "Macroalgal cover (%)", colour = "Planting treatment") +
  coord_cartesian(ylim = c(-5, 15), xlim = c(0, 60)) +
  stat_poly_line() +
  #stat_poly_eq() 
  stat_correlation(use_label(c("R", "R2", "P", "n")), small.r = TRUE, small.p = TRUE)
 
Figure_2

ggsave("Figure_2.pdf",plot = Figure_2, device = "pdf",
       width = 15,
       height = 15,
       units = c("cm"),
       dpi = 300)

#modeled effect of treatment on shoot count
seed_poisson <- glmer(Shoot_count_m2  ~ Treatment*Algae + (1|Site),  data = subset(seed_method, Shoot_count_m2 > 0), family=poisson(link = "log"))
summary(seed_poisson)
rsquared(seed_poisson)
seed_poisson

drop1(seed_poisson)

#binomial part
seed_binomial <- glmer((Shoot_count_m2!=0) ~ Treatment + Algae + (1|Site),  data = seed_method, family=binomial(link="log"))
summary(seed_binomial)
rsquared(seed_binomial)
drop1(seed_binomial)


# predicted probability of mf
seed_method$prob_seed <- predict(seed_binomial, type = "response")

# predicted mean mf for sites with fish
seed_method_bin <- subset(seed_method, Shoot_count > 0)
seed_method_bin$seed_growth <- predict(seed_poisson, type = "response")

# combine data frames
seed_method2 <- left_join(seed_method, select(seed_method_bin, c(ID, seed_growth)), by = "ID")

# replace NA with 0
seed_method2$seed_growth <- ifelse(is.na(seed_method2$seed_growth), 0, seed_method2$seed_growth)

# calculate predicted mf for all
seed_method2 <- seed_method2 %>% mutate(seed_pred = prob_seed * seed_growth)
seed_method2$seed_pred <- as.integer(seed_method2$seed_pred)

seed_method2 <- seed_method2 %>% mutate(prob_seed_perc = prob_seed * 100)

p1 <- ggplot(seed_method2, aes(x = Treatment, y = prob_seed_perc))  +
  geom_point(aes(colour=Treatment), position = position_jitter(width = 0.2), alpha = .4, size = 4) + 
  stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75), colour = "#233B43") +
  stat_summary(fun.data = mean_cl_boot, na.rm =TRUE, geom = "errorbar", width = .3, position = position_dodge(width = .75), colour = "#233B43") +
  #coord_cartesian(ylim = c(0, 0.8)) +
  labs(y = "Likelihood of emergence (%)", x = "") +
  #theme_bw(base_size = 12) +
  theme(legend.position = "none")
p1

p2<-ggplot(seed_method2, aes(x = Algae, y = prob_seed_perc))  +
  geom_point(aes(colour=Treatment), position = position_jitter(width = 0.2), alpha = .4, size = 4) +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), linetype = 4, colour = "#233B43") +
  #stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75)) +
  #stat_summary(fun.data = mean_se, na.rm =TRUE, geom = "errorbar", width = .3, position = position_dodge(width = .75)) +
  #coord_cartesian(ylim = c(-0.5, 3.5)) +
  labs(y = "Likelihood of emergence (%)", x = "", color = "Planting treatment") +
  theme(legend.position = "none")+
  facet_grid(~Treatment)

p2

p3 <- ggplot(subset(seed_method2, Max_leaf > 0), aes(x = Treatment, y = seed_growth))  +
  geom_point(aes(colour=Treatment),position = position_jitter(width = 0.2), alpha = .4, size = 4) + 
  stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75), colour = "#233B43") +
  stat_summary(fun.data = mean_cl_boot, na.rm =TRUE, geom = "errorbar", width = .3, position = position_dodge(width = .75), colour = "#233B43") +
  #coord_cartesian(ylim = c(0, 0.8)) +
  labs(y = bquote('Shoot density '("per m"^2)), x = "") +
  #theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 40))

p3

p4<-ggplot(subset(seed_method2, Max_leaf > 0), aes(x = Algae, y = seed_growth))  +
  geom_point(aes(colour=Treatment), position = position_jitter(width = 0.2), alpha = .4, size = 4) +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), linetype = 4, colour = "#233B43") +
  #stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75)) +
  #stat_summary(fun.data = mean_se, na.rm =TRUE, geom = "errorbar", width = .3, position = position_dodge(width = .75)) +
  #coord_cartesian(ylim = c(-0.5, 3.5)) +
  labs(y = bquote('Shoot density '("per m"^2)), x = "", color = "Planting treatment") +
  coord_cartesian(ylim = c(0, 40)) +
  theme(legend.position = "none") + 
  facet_grid(~Treatment)

p4

p5 <- ggplot(seed_method2, aes(x = Treatment, y = seed_pred))  +
  geom_point(aes(colour=Treatment), position = position_jitter(width = 0.2), alpha = .4, size = 4) + 
  stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75), colour = "#233B43") +
  stat_summary(fun.data = mean_cl_boot, na.rm =TRUE, geom = "errorbar", width = .3, position = position_dodge(width = .75), colour = "#233B43") +
  #coord_cartesian(ylim = c(-0.5, 3.5)) +
  labs(y = "Seagrass emergence success", x = "Planting treatment") +
  #theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 12)) +
  scale_y_continuous(
    labels = label_number(accuracy = 1))

p5

p6 <- ggplot(seed_method2, aes(x = Algae, y = seed_pred))  +
  geom_point(aes(colour=Treatment), position = position_jitter(width = 0.2), alpha = .4, size = 4) +
  stat_smooth(method = "glm", formula = y ~ poly(x, 2), linetype = 4, colour = "#233B43") +
  #stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75)) +
  #stat_summary(fun.data = mean_se, na.rm =TRUE, geom = "errorbar", width = .3, position = position_dodge(width = .75)) +
  #coord_cartesian(ylim = c(-0.5, 3.5)) +
  labs(y = "Seagrass emergence success", x = "Macroalgal cover (%)", color = "Planting treatment") +
  coord_cartesian(ylim = c(0, 12)) +
  scale_y_continuous(
    labels = label_number(accuracy = 1)) +
  theme(legend.position = "bottom") + 
  facet_grid(~Treatment)
p6
p_legend<-get_legend(p6)
p6 <- p6 + theme(legend.position = "none")

Figure_3 <- plot_grid(plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, nrow=3, labels = c("A", "B", "C", "D", "E", "F"), align = "v", axis = "l"), p_legend, nrow =2, rel_heights = c(1,0.05))
Figure_3

ggsave("Figure_3_alternative.tiff",plot = Figure_3, device = "tiff",
       width = 20,
       height = 25,
       units = c("cm"),
       dpi = 300)

####leaf length####
#modeled effect of treatment on leaf length
seed_poisson <- glmer(Max_leaf  ~ Treatment + Algae + (1|Site),  data = subset(seed_method, Shoot_count > 0), family = gaussian(link = "log"))
summary(seed_poisson)
rsquared(seed_poisson)
seed_poisson

drop1(seed_poisson)

#binomial part
seed_binomial <- glmer((Shoot_count_m2!=0) ~ Treatment + Algae + (1|Site),  data = seed_method, family=binomial(link = "log"))
summary(seed_binomial)
rsquared(seed_binomial)
drop1(seed_binomial)


# predicted probability of mf
seed_method$prob_seed <- predict(seed_binomial, type = "response")

# predicted mean mf for sites with fish
seed_method_bin <- subset(seed_method, Shoot_count > 0)
seed_method_bin$seed_growth <- predict(seed_poisson, type = "response")

# combine data frames
seed_method2 <- left_join(seed_method, select(seed_method_bin, c(ID, seed_growth)), by = "ID")

# replace NA with 0
seed_method2$seed_growth <- ifelse(is.na(seed_method2$seed_growth), 0, seed_method2$seed_growth)

# calculate predicted mf for all
seed_method2 <- seed_method2 %>% mutate(seed_pred = prob_seed * seed_growth)

seed_method2 <- seed_method2 %>% mutate(prob_seed_perc = prob_seed * 100)

seed_method2 <- seed_method2 %>% mutate(seed_pred_large = prob_seed_perc * seed_growth)

p1 <- ggplot(seed_method2, aes(x = Treatment, y = prob_seed_perc))  +
  geom_point(aes(colour=Treatment), position = position_jitter(width = 0.2), alpha = .4, size = 4) + 
  stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75), colour = "#233B43") +
  stat_summary(fun.data = mean_cl_boot, na.rm =TRUE, geom = "errorbar", width = .3, position = position_dodge(width = .75), colour = "#233B43") +
  #coord_cartesian(ylim = c(0, 0.8)) +
  labs(y = "Likelihood of emergence (%)", x = "") +
  #theme_bw(base_size = 12) +
  theme(legend.position = "none")
p1

p2<-ggplot(seed_method2, aes(x = Algae, y = prob_seed_perc))  +
  geom_point(aes(colour=Treatment), position = position_jitter(width = 0.2), alpha = .4, size = 4) +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), linetype = 4, colour = "#233B43") +
  #stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75)) +
  #stat_summary(fun.data = mean_se, na.rm =TRUE, geom = "errorbar", width = .3, position = position_dodge(width = .75)) +
  #coord_cartesian(ylim = c(-0.5, 3.5)) +
  labs(y = "Likelihood of emergence (%)", x = "") +
  theme(legend.position = "none")
p2

p3 <- ggplot(subset(seed_method2, Max_leaf > 0), aes(x = Treatment, y = seed_growth))  +
  geom_point(aes(colour=Treatment), position = position_jitter(width = 0.2), alpha = .4, size = 4) + 
  stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75), colour = "#233B43") +
  stat_summary(fun.data = mean_cl_boot, na.rm =TRUE, geom = "errorbar", width = .3, position = position_dodge(width = .75), colour = "#233B43") +
  #coord_cartesian(ylim = c(0, 0.8)) +
  labs(y = "Max leaf length (mm)", x = "") +
  #theme_bw(base_size = 12) +
  theme(legend.position = "none")

p3

p4<-ggplot(subset(seed_method2, Max_leaf > 0), aes(x = Algae, y = seed_growth))  +
  geom_point(aes(colour=Treatment), position = position_jitter(width = 0.2), alpha = .4, size = 4) +
  geom_smooth(method = "glm", formula = y ~ poly(x, 2), linetype = 4, colour = "#233B43") +
  #stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75)) +
  #stat_summary(fun.data = mean_se, na.rm =TRUE, geom = "errorbar", width = .3, position = position_dodge(width = .75)) +
  #coord_cartesian(ylim = c(-0.5, 3.5)) +
  labs(y = "Max leaf length (mm)", x = "") +
  theme(legend.position = "none")

p4

p5 <- ggplot(seed_method2, aes(x = Treatment, y = seed_pred))  +
  geom_point(aes(colour=Treatment), position = position_jitter(width = 0.2), alpha = .4, size = 4) + 
  stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75), colour = "#233B43") +
  stat_summary(fun.data = mean_cl_boot, na.rm =TRUE, geom = "errorbar", width = .3, position = position_dodge(width = .75), colour = "#233B43") +
  #coord_cartesian(ylim = c(-0.5, 3.5)) +
  labs(y = "Seagrass growth", x = "Planting treatment") +
  #theme_bw(base_size = 12) +
  theme(legend.position = "none")
p5

p6 <- ggplot(seed_method2, aes(x = Algae, y = seed_pred))  +
  geom_point(aes(colour=Treatment), position = position_jitter(width = 0.2), alpha = .4, size = 4) +
  stat_smooth(method = "glm", formula = y ~ poly(x, 2), linetype = 4, colour = "#233B43") +
  #stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75)) +
  #stat_summary(fun.data = mean_se, na.rm =TRUE, geom = "errorbar", width = .3, position = position_dodge(width = .75)) +
  #coord_cartesian(ylim = c(-0.5, 3.5)) +
  labs(y = "Seagrass growth", x = "Macroalgal cover (%)") +
  theme(legend.position = "bottom")
p6

p6
p_legend<-get_legend(p6)
p6 <- p6 + theme(legend.position = "none")

Figure_4 <- plot_grid(plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, nrow=3, labels = c("A", "B", "C", "D", "E", "F"), align = "v", axis = "l"), p_legend, nrow =2, rel_heights = c(1,0.05))
Figure_4

ggsave("Figure_4.tiff",plot = Figure_4, device = "tiff",
       width = 20,
       height = 25,
       units = c("cm"),
       dpi = 300)


####Supplimentary plots####

#plot of shood density + shoot length + algae - quadrats

Sa <- ggplot(seed_method, aes(x = Treatment, y = Shoot_count_m2))  +
  geom_jitter(alpha = .8, size = 4, width = 0.2) + 
  stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75), colour = "#233B43") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .3, position = position_dodge(width = .75), colour = "#233B43") +
  labs(y = bquote('Shoot density '("per m"^2)), x = "") +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 60)) +
  scale_y_continuous(labels = label_number(accuracy = 1))
Sa

Sb <- ggplot(seed_method, aes(x = Treatment, y = Max_leaf))  +
  geom_jitter(alpha = .8, size = 4, width = 0.2) + 
  stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75), colour = "#233B43") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .3, position = position_dodge(width = .75), colour = "#233B43") +
  labs(y = "Max leaf length (mm)", x = "") +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 700))
Sb
Sc <- ggplot(seed_method, aes(x = Treatment, y = Algae))  +
  geom_jitter(alpha = .8, size = 4, width = 0.2) +  
  stat_summary(fun = mean, na.rm = TRUE, geom = "point", size = 8, position = position_dodge(width = .75), colour = "#233B43") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .3, position = position_dodge(width = .75), colour = "#233B43") +
  labs(y = "Macroalgal cover (%)", x = "Planting treatment") +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 100))
Sc
Supp_Figure_1 <- plot_grid(Sa,Sb,Sc, nrow = 3, ncol = 1, labels = c('A', 'B', 'C'), align = "v")
Supp_Figure_1


ggsave("Figure_S1.tiff",plot = Supp_Figure_1, device = "tiff",
       width = 20,
       height = 25,
       units = c("cm"),
       dpi = 300)

