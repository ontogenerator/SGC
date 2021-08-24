
# chase flights detection algorithm
options(digits.secs = 3) # change options to see milliseconds in DateTime 

library(tidyverse)
library(lubridate)
library(ggpubr)
library(broom.mixed)
library(MCMCglmm)
library(PlayerRatings)
library(ggbeeswarm)
library(ggpattern)
folder <- "analysis/data/"

# File too large to upload to github, so needs to be generated first with script in load.R
allnights <- read.csv2(file = paste0(folder, "EventData.csv"),
                                header = TRUE, sep = ";", dec = ".", na.strings = "NA") 

dur_epsilon <- 200 # 200 ms as a definition of a the upper threshold for the duration of a TR detection caused by chasing event

# are LS ever preceded by Reader events?
# very rarely and never looks like a chase event!!
# allnights_n <- allnights_n %>% 
#   group_by(night) %>% 
#   mutate(postreader = ifelse(lag(loc) == loc & str_detect(lag(unitLabel), "Reader") & str_detect(unitLabel, "LS"), TRUE, FALSE))

allnights_n <- allnights %>% 
  group_by(night) %>% 
  mutate(DateTime = as.POSIXct(DateTime, tz = "UTC"),
         n_event = 1:n(),
         chase = case_when(
           n_event < 3 ~ FALSE, # first two events cannot be chases by definition
           is.na(IdLabel) ~ FALSE, # not a chase if corrupted RFID number detected that cannot be assigned to any ind
           eventDuration >= dur_epsilon ~ FALSE, # not a chase if the duration is higher than the epsilon
           !str_detect(unitLabel, "Re") ~ FALSE, # not a chase if not a Reader event
           !str_detect(lag(unitLabel), "LS") ~ FALSE, # not a chase if previous event not an IR sensor event
           loc != lag(loc) | loc != lag(loc, 2) ~ FALSE, # not a chase if previous two events at different flower
           !str_detect(lag(unitLabel, 2), "Co") ~ FALSE,  # not a chase if next-to-last event prior was not a CondMod
           TRUE ~ TRUE # otherwise chase visit
           ),
         chased = ifelse(chase == TRUE, lag(IdLabel, 2), NA), 
         # dur_chase_seq = ifelse(chase == TRUE, 
         #                        lag(eventDuration, 2)/1000 -
         #                          as.numeric(DateTime - lag(DateTime, 2)), NA),
         patch = case_when(
           loc < 6 ~ 1L,
           loc > 5 ~ 2L,
           TRUE ~ NA_integer_),
         group = factor(group, levels = c("mixed1", "mixed2", "mixed3", "mixed4", "6m", "6f")),
         sex = factor(sex, levels = c("m", "f"))) %>%
  select(chase, chased, everything())

phase_durs <- allnights_n %>% 
  filter(n_event == 1 | SystemMsg == "switch") %>% 
  group_by(night) %>% 
  summarise(hours_ph1 = round(max(DateTime) - min(DateTime), 1),
            hours_ph2 = 12 - hours_ph1) %>% 
  pivot_longer(-night, names_to = "phase", values_to = "phase_dur", names_prefix = "hours_ph") %>% 
  mutate(phase = as.integer(phase))

vis_summaries <- allnights_n %>% 
  filter(eventDuration > 300, choice == TRUE, !is.na(IdLabel)) %>% 
  group_by(night, group, group_night, IdLabel, phase, cond, patch, sex) %>% 
  summarise(n = n(), vol_consumed = sum(vol) / 1000) %>% 
  left_join(phase_durs) %>% 
  mutate(vis_hour = n/as.numeric(phase_dur),
         ml_hour = vol_consumed/as.numeric(phase_dur))

rev_count <- vis_summaries %>% 
  filter(cond == "test") %>% 
  ungroup() %>% 
  count(group, group_night) %>% 
  group_by(group) %>% 
  mutate(n = n():1) %>% 
  rename(rev_night = n) # reverse count the nights to find the last two nights (different for each group)

vis_summaries <- vis_summaries %>% 
  ungroup() %>% 
  left_join(rev_count)

write.table(vis_summaries, file = paste0(folder, "VisitSummary.csv"), sep = ";", row.names = FALSE)

chase_summ <- allnights_n %>% 
  filter(!is.na(IdLabel), !is.na(cond)) %>% 
  group_by(group_night, cond, IdLabel, group, sex, phase) %>% 
  summarise(n_detections = n(),
            n_chases = sum(chase),
            prop_chases = n_chases/n_detections)

chased_n <- allnights_n %>% 
  ungroup() %>%
  filter(!is.na(chased)) %>% 
  count(cond, group_night, chased, group, phase) %>%
  rename(IdLabel = chased,
         n_chased = n) 

bat_info <- read.csv2(file = paste0(folder, "metadata/conditions.csv"),
                                          header = TRUE, sep = ";", dec = ".", na.strings = "NA") %>%
  filter(group_night == 5) %>% 
  group_by(IdLabel) %>% 
  summarise(weight = mean(weight),
            sex = sex)

# Glicko rating analysis, based on last two nights only

only_chases <- allnights_n %>%
  ungroup() %>% 
  filter(chase == TRUE) %>% 
  select(cond, group_night, group, IdLabel, chased) %>% 
  mutate(outcome = 1) %>% 
  left_join(rev_count) 

gl <- only_chases %>% 
  filter(rev_night < 3, cond == "test") %>% 
  select(-rev_night, -cond) %>% 
  nest_by(group) %>% 
  mutate(glicko = list(glicko(data))) %>% 
  select(-data) %>% 
  mutate(IdLabel = list(pluck(glicko, 1, "Player")),
         glicko_rating = list(pluck(glicko, 1, "Rating")))

rankings <- gl %>% 
  select(-glicko) %>% 
  unnest(c(IdLabel, glicko_rating)) %>% 
  mutate(glicko_rank = 1:n())

chase_summ <- chase_summ %>% 
  left_join(chased_n) %>%
  left_join(bat_info) %>% 
  left_join(rankings) %>% 
  replace_na(list(n_chased = 0)) %>% 
  mutate(prop_chased = n_chased / n_detections,
         phase = factor(phase), 
         sex = factor(sex))

write.table(chase_summ, file = paste0(folder, "ChaseSummary.csv"), sep = ";", row.names = FALSE)

# function for plotting autocorrelations from https://github.com/tmalsburg/MCMCglmm-intro
plot.acfs <- function(x) {
  n <- dim(x)[2]
  par(mfrow = c(ceiling(n/2), 2), mar = c(3, 2, 3, 0))
  for (i in 1:n) {
    acf(x[, i], lag.max = 100, main = colnames(x)[i])
    grid()
  }
}

mcmc_to_tibble <- function(object) { # extract the model summary into a tibble 
  nF <- object$Fixed$nfl
  nL <- object$Fixed$nll
  
  solutions <- cbind(colMeans(object$Sol[, 1:nF, drop = FALSE]), 
                     coda::HPDinterval(object$Sol[, 1:nF, drop = FALSE]), 
                     effectiveSize(object$Sol[, 1:nF, drop = FALSE]), 2 * 
                       pmax(0.5/dim(object$Sol)[1], 
                            pmin(colSums(object$Sol[, 1:nF, drop = FALSE] > 0)/dim(object$Sol)[1],
                                 1 - colSums(object$Sol[, 1:nF, drop = FALSE] > 0)/dim(object$Sol)[1])))
  if (nL > 0) {
    solutions <- rbind(solutions, cbind(colMeans(object$Lambda), 
                                        coda::HPDinterval(object$Lambda), effectiveSize(object$Lambda), 
                                        2 * pmax(0.5/dim(object$Lambda)[1],
                                                 pmin(colSums(object$Lambda > 0)/dim(object$Lambda)[1],
                                                      1 - colSums(object$Lambda > 0)/dim(object$Sol)[1]))))
  }
  colnames(solutions) <- c("estimate", "conf.low", 
                           "conf.high", "eff.samp", "pMCMC")
  solutions %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "term") %>% 
    mutate(across(starts_with("est") | starts_with("conf"), ~round(., 2)))
}

prior.1 <- list(R = list(V = 1, nu = 0.002),
                G = list(G1 = list(V = diag(1), nu = 1, 
                                   alpha.mu = c(0), alpha.V = diag(1)*10000),
                         G2 = list(V = diag(1), nu = 1, 
                                   alpha.mu = c(0), alpha.V = diag(1)*10000)))

#prior.1a<-list(R = list(V=1, nu = 0.002), G = list(G1 = list(V = diag(2), nu = 0.002), 
#G2 = list(V = 1, nu = 0.002)))

set.seed(42)
mcmc_chases <- MCMCglmm(cbind(n_chases, n_detections - n_chases) ~ sex * phase + weight, 
                     random = ~IdLabel + group,
                     data = as.data.frame(chase_summ), family = "multinomial2", pr = TRUE,
                     prior = prior.1, verbose = FALSE, saveX = TRUE, saveZ = TRUE,
                     nitt = 200E3, thin = 100, burnin = 100E3)
summary(mcmc_chases)

geweke.plot(mcmc_chases$VCV)
plot.acfs(mcmc_chases$VCV)
plot(mcmc_chases)

est_chases <- mcmc_to_tibble(mcmc_chases) %>% 
  mutate(Model = "Chasing")

set.seed(123)
mcmc_chased <- MCMCglmm(cbind(n_chased, n_detections - n_chased) ~ sex * phase + weight, 
                        random = ~IdLabel + group,
                        data = as.data.frame(chase_summ), family = "multinomial2", pr = TRUE,
                        prior = prior.1, verbose = FALSE, saveX = TRUE, saveZ = TRUE,
                        nitt = 200E3, thin = 100, burnin = 100E3)
summary(mcmc_chased)
geweke.plot(mcmc_chased$VCV)
plot(mcmc_chased)


est_chased <- mcmc_to_tibble(mcmc_chased) %>% 
  mutate(Model = "Being chased")
est_chases <- est_chases %>% 
  bind_rows(est_chased) %>% 
  mutate(CI = paste0("(", conf.low, ", ", conf.high, ")"),
         term = str_replace(term, "sexf", "sex (female)"),
         term = str_replace(term, "phase2", "condition (distributed)"))

write.table(est_chases, file = paste0(folder, "mcmcChases.csv"), sep = ";", row.names = FALSE)

# 
vol_sd_time <- vis_summaries %>%
  filter(group_night > 0) %>%
  group_by(group, IdLabel, phase, sex, group_night) %>%
  summarise(vol_hr = sum(vol_consumed)/sum(as.numeric(phase_dur))) %>%
  group_by(group, group_night, phase, sex) %>%
  summarise(sd_consumed = sd(vol_hr)) %>%
  mutate(group2 = paste(group, sex),
         phase = factor(phase)) %>% 
  ungroup() %>% 
  mutate(night = scale(group_night, center = TRUE, scale = FALSE))

write.table(vol_sd_time, file = paste0(folder, "IntakeTime.csv"), sep = ";", row.names = FALSE)


prior.2 <- list(R = list(V = 1, nu = 0.002),
                G = list(G1 = list(V = diag(2), nu = 0.002)))

#prior.1a<-list(R = list(V=1, nu = 0.002), G = list(G1 = list(V = diag(2), nu = 0.002), 
#G2 = list(V = 1, nu = 0.002)))

set.seed(815) 
mcmc_vol_time <- MCMCglmm(sd_consumed ~ sex + phase + night +
                            sex*phase*night, 
                        random = ~us(1 + night):group,
                        data = as.data.frame(vol_sd_time), family = "gaussian", pr = TRUE,
                        prior = prior.2, verbose = FALSE, saveX = TRUE, saveZ = TRUE,
                        nitt = 200E4, thin = 1000, burnin = 100E3)
summary(mcmc_vol_time)

geweke.plot(mcmc_vol_time$VCV)
plot(mcmc_vol_time)

est_vol_time <- mcmc_to_tibble(mcmc_vol_time) %>% 
  mutate(CI = paste0("(", conf.low, ", ", conf.high, ")"),
         term = str_replace(term, "sexf", "sex (female)"),
         term = str_replace(term, "phase2", "condition (distributed)"),
         term = str_replace(term, "group_night", "night"))

write.table(est_vol_time, file = paste0(folder, "mcmcIntakeTime.csv"), sep = ";", row.names = FALSE)


mode_HPD <- function(mcmc) {
  round(c(posterior.mode(mcmc), HPDinterval(mcmc)), 2) %>% set_names(c("estimate", "conf.low", "conf.high"))
}
# mp2 <- int + phase2
# fp2 <- int + phase2 + sexf + sexfphase2
# mslp2 <- night + phase2night
# fslp2 <- night + phase2night + sexfnight + sexfnightphase2

sex_phase2 <- mcmc_vol_time$Sol[, "sexf"] + mcmc_vol_time$Sol[, "sexf:phase2"] # checks out
night_phase2 <- mcmc_vol_time$Sol[, "night"] + mcmc_vol_time$Sol[, "phase2:night"] # checks out
females_phase2 <- mcmc_vol_time$Sol[, "phase2"] + mcmc_vol_time$Sol[, "sexf:phase2"] # checks out
sex_night_phase2 <- -(mcmc_vol_time$Sol[, "sexf:phase2:night"] +  mcmc_vol_time$Sol[, "sexf:night"]) # checks out
females_nightph1 <- mcmc_vol_time$Sol[, "night"] + mcmc_vol_time$Sol[, "sexf:night"] # checks out
females_nightph2 <- mcmc_vol_time$Sol[, "night"] + mcmc_vol_time$Sol[, "sexf:night"] + 
  mcmc_vol_time$Sol[, "phase2:night"] +  mcmc_vol_time$Sol[, "sexf:phase2:night"] # double check!
females_night_contrast <- females_nightph2 - females_nightph1

intake_interactions <- mode_HPD(females_phase2) %>% # for females, sign. difference between resource conditions
  rbind(mode_HPD(night_phase2)) %>%  # slope of night for males equals zero in phase 2 (distributed resource condition)
  rbind(mode_HPD(sex_phase2)) %>% # no difference between sexes in phase 2 (distributed resource condition)
  rbind(mode_HPD(sex_night_phase2)) %>% # no difference in slopes between sexes in phase 2 either
  rbind(mode_HPD(females_nightph1)) %>% # no effect of night in females in clumped
  rbind(mode_HPD(females_nightph2)) %>% # nor in distributed condition
  rbind(mode_HPD(females_night_contrast)) %>% # no difference in slope for females between conditions
  as_tibble() %>%
  mutate(term = c("females_phase2", "night_phase2", "sex_phase2",
                  "sex_night_phase2", "females_night", "females_night_phase2",
                  "females_night_contrast")) %>% 
  select(term, everything())

predicted <- predict(mcmc_vol_time, type = "response", posterior = "mode") %>% as.vector()
predicted <- tibble(predicted = predicted)
vol_sd_time %>% 
  # filter(!str_detect(group, "6")) %>% 
  bind_cols(predicted) %>% 
  ggplot(aes(group_night, predicted, color = group)) +
  # ggplot(aes(group_night, sd_consumed, color = group)) +
  geom_point() +
  geom_smooth(aes(group = group2, linetype = sex), method = "lm") +
  facet_grid(. ~ phase, labeller = labeller(phase = c(`1` = "Clumped resource condition",
                                                      `2` = "Distributed resource condition"))) +
  scale_color_viridis_d(option = "turbo", labels = c("Mixed Group 1", "Mixed Group 2",
                                                     "Mixed Group 3", "Mixed Group 4",
                                                     "Males-only Group", "Females-only Group"),
                        direction = -1, name = NULL) +
  scale_linetype_manual(labels = c("Male", "Female"), name = NULL, values=c(2, 1)) +
  scale_x_continuous(breaks = 1:9) +
  ylim(c(0, 1)) +
  theme_serif() +
  labs(x = "Experimental night", y = "Standard deviation of\nnectar intake")

write.table(intake_interactions, file = paste0(folder, "IntakeInteractions.csv"), sep = ";", row.names = FALSE)

onlyrewarded <- allnights_n %>% 
  left_join(statuses) %>% # statuses calculation in Rmd file!
  ungroup() %>% 
  filter(vol > 0, cond == "test", !is.na(IdLabel)) %>% 
  select(group_night, IdLabel, status, group, DateTime, loc, phase)

write.table(onlyrewarded, 
            file = paste0(folder,"RewardedVisits.csv"), sep = ";", row.names = FALSE)

### Figures for slides in presentations:

cons_clumped <- cons_group1 %>%
  filter(phase == 1) %>%
  ggplot(aes(x = group_night, y = vol_consumed, shape = Individual, group = Individual)) +
  labs(x = "Night", title = "Clumped condition - one patch",
       y = "Nectar intake [mL]") +
  ylim(c(0, 18)) +
  scale_shape_manual(values = c(21, 22, 24, 21, 22, 24), guide = FALSE) +
  scale_color_manual(values = c(rgb(red = 46, green = 117, blue = 182, max = 255),
                                rgb(red = 184, green = 13, blue = 72, max = 255)),
                     name = "", labels = c("Male", "Female")) +
  scale_fill_manual(name = "", values = c(rgb(red = 46, green = 117, blue = 182, max = 255),
                                          rgb(red = 184, green = 13, blue = 72, max = 255)), guide = FALSE) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7), limits = c(0.8, 7)) +
  geom_line(aes(color = sex)) +
  geom_point(aes(color = sex, fill = sex), size = 3) +
  theme_serif() +
  theme(legend.position = c(0.2, 0.98))

cons_distr <- cons_group1 %>% 
  filter(phase == 2, group_night != 4) %>%
  ggplot(aes(x = group_night, y = vol_consumed)) +
  labs(x = "Night", title = "Distributed condition - two patches",
       y = "") +
  scale_shape_manual(values = c(21, 22, 24, 21, 22, 24), guide = FALSE) +
  scale_color_manual(values = c(rgb(red = 46, green = 117, blue = 182, max = 255),
                                rgb(red = 184, green = 13, blue = 72, max = 255)), guide = FALSE) +
  scale_fill_manual(name = "", values = c(rgb(red = 46, green = 117, blue = 182, max = 255),
                                          rgb(red = 184, green = 13, blue = 72, max = 255)), guide = FALSE) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7), limits = c(0.8, 7)) +
  geom_line(aes(group = Individual, color = sex)) +
  geom_point(aes(shape = Individual, color = sex, fill = sex), size = 3) +
  ylim(c(0, 18)) +
  theme_serif()

ggarrange(cons_clumped, cons_distr,
          font.label = list(size = 18, family = "serif"))

fig_chases <- chase_summ %>% 
  filter(phase == 1) %>% 
  group_by(sex, group, IdLabel) %>% 
  summarise(prop_chases = mean(prop_chases)) %>% 
  ggplot() +
  geom_beeswarm(aes(group, prop_chases, color = sex),
              alpha = 0.7, cex = 2, size = 4) +
  # geom_boxplot(aes(group, prop_chases, fill = sex),
  #              width = 0.3, size = 0.6, alpha = 0.7,
  #              position = position_dodge(width = 1)) +
  labs(title = "Chasing", x = "", y = "Proportion chasing") +
  theme_serif() +
  scale_color_manual(values = c(rgb(red = 46, green = 117, blue = 182, max = 255),
                                rgb(red = 184, green = 13, blue = 72, max = 255)),
                     name = "", labels = c("Male", "Female")) +
  # scale_fill_manual("", values = c("grey24", "white"), guide = FALSE) +
  scale_x_discrete(labels = c("Mixed \n Group1", "Mixed \n Group2", "Mixed \n Group3",
                              "Mixed \n Group4", "Male \n Group", "Female \n Group")) +
  scale_y_continuous(limits = c(0, 0.016), labels = c(0, "", 0.01, "")) +
  theme(axis.ticks.x = element_line(size = 1.2),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 18),
        legend.position = c(0.25, 0.9))

fig_chased <- chase_summ %>%  
  filter(phase == 1) %>% 
  group_by(sex, group, IdLabel) %>% 
  summarise(prop_chased = mean(prop_chased)) %>% 
  ggplot() +
  geom_beeswarm(aes(group, prop_chased, color = sex),
                alpha = 0.7, cex = 2, size = 4) +
  labs(title = "Being chased", x = "", y = "Proportion chased") +
  theme_serif() +
  scale_color_manual(values = c(rgb(red = 46, green = 117, blue = 182, max = 255),
                                rgb(red = 184, green = 13, blue = 72, max = 255)),
                     guide = FALSE) +
  scale_x_discrete(labels = c("Mixed \n Group1", "Mixed \n Group2", "Mixed \n Group3",
                              "Mixed \n Group4", "Male \n Group", "Female \n Group")) +
  scale_y_continuous(limits = c(0, 0.016), labels = c(0, "", 0.01, "")) +
  theme(axis.ticks.x = element_line(size = 1.2),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 18))

ggarrange(fig_chases, fig_chased,
          font.label = list(size = 18, family = "serif"))


ggarrange(chase_males, chase_females,
          font.label = list(size = 18, family = "serif"),
          common.legend = TRUE, legend = "none")


stat.test2 <- cons_contrasts %>%
  group_by(partition, status) %>%
  t_test(vol_hr ~ phase) %>% 
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>% 
  add_xy_position(x = "phase", group = "status", fun = "min") 


ggboxplot(cons_contrasts, x = "phase", y = "vol_hr",
          color = "status", line.color = "gray", line.size = 0.4,
          palette = c("dodgerblue3","grey50","purple"),
          add = "jitter") +
  facet_grid(. ~ partition) +
  # stat_pvalue_manual(
  #   stat.test,  label = "p.adj.signif", tip.length = 0.02, hide.ns = TRUE
  # ) +
  stat_pvalue_manual(
    stat.test,  label = "p.adj", tip.length = 0.02, hide.ns = TRUE
  ) +
   stat_pvalue_manual(
    stat.test2,  label = "p.adj", tip.length = 0.02, hide.ns = TRUE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = "", y = bquote("Nectar intake ["~mL~h^-1*"]"), color = "") +
  theme_serif()



chase_males <- chase_nectar %>%
  filter(sex == "m") %>%
  ggplot(aes(prop_chases, vol_hr)) +
  geom_text(aes(label = glicko_rank, color = group), show.legend  = FALSE) +
  geom_point(aes(color = group), alpha = 0) +
  scale_color_viridis_d(option = "turbo", name = "Groups", 
                        labels = c("Mixed Group 1", "Mixed Group 2",
                                   "Mixed Group 3", "Mixed Group 4",
                                   "Group 6 Males", "Group 6 Females"),
                        drop = FALSE, direction = -1) +
  # geom_point(aes(prop_chases, vol_hr, shape = group), size = 3) +
  stat_ellipse(data = . %>%  filter(vol_hr > 0.75, prop_chases > 0.003),
               aes(prop_chases, vol_hr), level = 0.89, linetype = 2) +
  labs(x = "Proportion of chasing", 
       title = "Clumped resource condition \n Males",
       y = expression(atop("Nectar intake ["~mL~h^-1*"]", "(Over last two nights)"))) +
  # geom_vline(xintercept = 0.003, linetype = 3) +
  # geom_hline(yintercept = 0.75, linetype = 3) +
  ylim(0, 2) + 
  xlim(0, 0.013) +
  theme_serif() +
  theme(legend.text = element_text(size = 18, family = "serif"))  +
  guides(color = guide_legend(override.aes = list(alpha = 1)))


chase_females <- chase_nectar %>%
  filter(sex == "f") %>% 
  ggplot(aes(prop_chases, vol_hr)) +
  geom_text(aes(label = glicko_rank, color = group), show.legend  = FALSE) +
  geom_point(aes(color = group), alpha = 0) +
  scale_color_viridis_d(option = "turbo", name = "Groups", 
                        labels = c("Mixed Group 1", "Mixed Group 2",
                                   "Mixed Group 3", "Mixed Group 4",
                                   "Group 6 Males", "Group 6 Females"),
                        drop = FALSE, direction = -1) +
  labs(x = "Proportion of chasing", 
       title = "Clumped resource condition \n Females",
       y = bquote(atop("Nectar intake ["~mL~h^-1*"]", "(Over last two nights)"))) +
  # geom_vline(xintercept = 0.003, linetype = 3) +
  # geom_hline(yintercept = 0.75, linetype = 3) +
  ylim(0, 2) + 
  xlim(0, 0.013) +
  theme_serif() +
  theme(legend.text = element_text(size = 18, family = "serif")) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))


ggarrange(chase_males, chase_females, labels = c("A", "B"),
          font.label = list(size = 18, family = "serif"),
          common.legend = TRUE, legend = "right")

######## Analyses made to answer questions from reviewers

# defending on nights when switch happened earlier/later
allnights_n %>%
  ungroup() %>% 
  filter(cond == "test", chase == TRUE) %>% 
  mutate(hour = factor(hour(DateTime),
                       levels = c(15:23, 0:14))) %>% 
  count(group, night, group_night, hour) %>% 
  left_join(filter(phase_durs, phase == 1)) %>% 
  mutate(phase_dur = as.numeric(round(phase_dur)) + 1) %>% 
  ggplot(aes(hour, n, color = group, group = group)) +
  geom_point() +
  geom_line() +
  ylab("total_chases") +
  scale_x_discrete(breaks = c(16, 18, 20, 22, 0, 2, 4)) +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  facet_grid(group ~ group_night, labeller = label_both) +
  geom_vline(aes(xintercept = phase_dur), linetype = 2) +
  guides(color = FALSE)


allnights_n %>%
  ungroup() %>% 
  filter(cond == "test", chase == TRUE) %>% 
  left_join(statuses) %>% 
  filter(status == "Dominant males") %>% 
  mutate(hour = factor(hour(DateTime),
                       levels = c(15:23, 0:14))) %>% 
  count(group, night, group_night, hour) %>% 
  left_join(filter(phase_durs, phase == 1)) %>% 
  mutate(phase_dur = as.numeric(round(phase_dur)) + 1) %>% 
  ggplot(aes(hour, n, color = group, group = group)) +
  geom_point() +
  geom_line() +
  ylab("chases_by_dominant") +
  scale_x_discrete(breaks = c(16, 18, 20, 22, 0, 2, 4)) +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  facet_grid(group ~ group_night, labeller = label_both) +
  geom_vline(aes(xintercept = phase_dur), linetype = 2) +
  guides(color = FALSE)

allnights_n %>%
  ungroup() %>% 
  filter(group_night == 0, chase == TRUE) %>% 
  count(group) %>% 
  ggplot(aes(group, n, fill = group)) +
  geom_col() +
  scale_fill_viridis_d(option = "turbo", direction = -1, guide = FALSE) +
  scale_x_discrete(labels = c("Mixed \n Group1", "Mixed \n Group2", "Mixed \n Group3",
                              "Mixed \n Group4", "Male \n Group", "Female \n Group")) +
  labs(x = "", y = "Total chases in training") +
  theme_serif()

lastnight_f <- read.csv2(file = paste0(folder, "raw/group6f/", "Day10_AllActive_Group2-15.03.26.csv"),
                         sep = ";", dec = ",", header = TRUE,
                         fileEncoding = "UTF-16LE", as.is = TRUE, row.names = NULL) %>% 
  filter(!str_detect(DateTime, "#")) %>% 
  arrange(DateTime) %>% #sort chronologically
  mutate(DateTime = as.numeric(str_replace(DateTime, ",", ".")),
         DateTime = as.POSIXct(DateTime * (60 * 60 * 24),
                               origin = "1899-12-30", tz = "UTC"))

sex_info <- bat_info %>% 
  select(-weight) %>% 
  rename(chased = IdLabel,
         sex_chased = sex)

only_chases %>% 
  left_join(bat_info) %>% 
  left_join(sex_info) %>% 
  filter(!str_detect(group, "6")) %>% 
  group_by(IdLabel, sex, group) %>% 
  summarise(prop_chased_males = mean(sex_chased == "m"),
            prop_chased_females = mean(sex_chased == "f")) %>% 
  ggplot(aes(group, prop_chased_males, color = sex)) +
  geom_jitter(width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_serif() +
  scale_color_viridis_d(direction = -1)


only_chases %>% 
  filter(group_night > 0) %>% 
  left_join(bat_info) %>% 
  left_join(sex_info) %>% 
  left_join(statuses) %>% 
  count(group_night, IdLabel, sex, status, group) %>%
  ggplot(aes(group_night, n, color = status, group = IdLabel)) +
  geom_jitter(width = 0.2) +
  scale_x_continuous(breaks = 1:9) +
  ylab("Total number of chases") +
  geom_line() +
  facet_wrap(~group) +
  theme_serif() +
  scale_color_viridis_d(direction = -1)

allnights_n %>% 
  filter(group_night > 0) %>% 
  count(IdLabel, group_night, group) %>% 
  left_join(bat_info) %>% 
  left_join(statuses) %>% 
  ggplot(aes(group_night, n, color = status, group = IdLabel)) +
  geom_jitter(width = 0.2) +
  scale_x_continuous(breaks = 1:9) +
  ylab("Total number of detections") +
  geom_line() +
  facet_wrap(~group) +
  theme_serif() +
  scale_color_viridis_d(direction = -1)



#being chased vs nectar intake

chase_nectar <- chase_summ %>% 
  filter(phase == 1, cond == "test") %>% 
  group_by(group, sex, IdLabel, glicko_rank, glicko_rating) %>% 
  summarise(prop_chases = sum(n_chases)/sum(n_detections),
            prop_chased = sum(n_chased)/sum(n_detections),
            n_chases = sum(n_chases),
            n_chased = sum(n_chased)) %>% 
  left_join(intake_lastnights_clumped)


chased_males <- chase_nectar %>%
  filter(sex == "m") %>% 
  ggplot() +
  geom_point(aes(prop_chased, vol_hr, shape = group), size = 3) +
  stat_ellipse(data = . %>%  filter(vol_hr > 0.75, prop_chases > 0.003),
               aes(prop_chases, vol_hr), level = 0.89, linetype = 2) +
  scale_shape_manual("Groups", labels = c("Mixed Group 1", "Mixed Group 2",
                                          "Mixed Group 3", "Mixed Group 4",
                                          "Males-only Group", "Females-only Group"),
                     values = c(1, 2, 16, 17, 8, 4), drop = FALSE) +
  labs(x = "Proportion of being chased", 
       title = "Clumped resource condition \n Males",
       y = expression(atop("Nectar intake ["~mL~h^-1*"]", "(Over last two nights)"))) +
  # geom_vline(xintercept = 0.003, linetype = 3) +
  # geom_hline(yintercept = 0.75, linetype = 3) +
  ylim(0, 2) + 
  xlim(0, 0.013) +
  theme_serif() +
  theme(legend.text = element_text(size = 18, family = "serif"))



chased_females <- chase_nectar %>%
  filter(sex == "f") %>% 
  ggplot(aes(prop_chased, vol_hr, shape = group)) +
  geom_point(size = 3) +
  scale_shape_manual("Groups", labels = c("Mixed Group 1", "Mixed Group 2",
                                          "Mixed Group 3", "Mixed Group 4",
                                          "Males-only Group", "Females-only Group"),
                     values = c(1, 2, 16, 17, 8, 4), drop = FALSE) +
  labs(x = "Proportion of being chased", 
       title = "Clumped resource condition \n Females",
       y = bquote(atop("Nectar intake ["~mL~h^-1*"]", "(Over last two nights)"))) +
  # geom_vline(xintercept = 0.003, linetype = 3) +
  # geom_hline(yintercept = 0.75, linetype = 3) +
  ylim(0, 2) + 
  xlim(0, 0.013) +
  theme_serif() +
  theme(legend.text = element_text(size = 18, family = "serif"))

ggarrange(chased_males, chased_females, labels = c("A", "B"),
          font.label = list(size = 18, family = "serif"),
          common.legend = TRUE, legend = "right")


chased_males <- chase_nectar %>%
  filter(sex == "m") %>% 
  ggplot() +
  geom_point(aes(n_chased, vol_hr, shape = group), size = 3) +
  stat_ellipse(data = . %>%  filter(vol_hr > 0.75, prop_chases > 0.003),
               aes(prop_chases, vol_hr), level = 0.89, linetype = 2) +
  scale_shape_manual("Groups", labels = c("Mixed Group 1", "Mixed Group 2",
                                          "Mixed Group 3", "Mixed Group 4",
                                          "Males-only Group", "Females-only Group"),
                     values = c(1, 2, 16, 17, 8, 4), drop = FALSE) +
  labs(x = "# being chased", 
       title = "Clumped resource condition \n Males",
       y = expression(atop("Nectar intake ["~mL~h^-1*"]", "(Over last two nights)"))) +
  # geom_vline(xintercept = 0.003, linetype = 3) +
  # geom_hline(yintercept = 0.75, linetype = 3) +
  ylim(0, 2) + 
  # xlim(0, 0.013) +
  theme_serif() +
  theme(legend.text = element_text(size = 18, family = "serif"))



chased_females <- chase_nectar %>%
  filter(sex == "f") %>% 
  ggplot(aes(n_chased, vol_hr, shape = group)) +
  geom_point(size = 3) +
  scale_shape_manual("Groups", labels = c("Mixed Group 1", "Mixed Group 2",
                                          "Mixed Group 3", "Mixed Group 4",
                                          "Males-only Group", "Females-only Group"),
                     values = c(1, 2, 16, 17, 8, 4), drop = FALSE) +
  labs(x = "# being chased", 
       title = "Clumped resource condition \n Females",
       y = bquote(atop("Nectar intake ["~mL~h^-1*"]", "(Over last two nights)"))) +
  # geom_vline(xintercept = 0.003, linetype = 3) +
  # geom_hline(yintercept = 0.75, linetype = 3) +
  ylim(0, 2) + 
  # xlim(0, 0.013) +
  theme_serif() +
  theme(legend.text = element_text(size = 18, family = "serif"))

ggarrange(chased_males, chased_females, labels = c("A", "B"),
          font.label = list(size = 18, family = "serif"),
          common.legend = TRUE, legend = "right")


# total n chases vs nectar intake


chases_males <- chase_nectar %>%
  filter(sex == "m") %>% 
  ggplot() +
  geom_point(aes(n_chases, vol_hr, shape = group), size = 3) +
  stat_ellipse(data = . %>%  filter(vol_hr > 0.75, prop_chases > 0.003),
               aes(prop_chases, vol_hr), level = 0.89, linetype = 2) +
  scale_shape_manual("Groups", labels = c("Mixed Group 1", "Mixed Group 2",
                                          "Mixed Group 3", "Mixed Group 4",
                                          "Males-only Group", "Females-only Group"),
                     values = c(1, 2, 16, 17, 8, 4), drop = FALSE) +
  labs(x = "# chasing", 
       title = "Clumped resource condition \n Males",
       y = expression(atop("Nectar intake ["~mL~h^-1*"]", "(Over last two nights)"))) +
  # geom_vline(xintercept = 0.003, linetype = 3) +
  # geom_hline(yintercept = 0.75, linetype = 3) +
  ylim(0, 2) + 
  # xlim(0, 0.013) +
  theme_serif() +
  theme(legend.text = element_text(size = 18, family = "serif"))



chases_females <- chase_nectar %>%
  filter(sex == "f") %>% 
  ggplot(aes(n_chases, vol_hr, shape = group)) +
  geom_point(size = 3) +
  scale_shape_manual("Groups", labels = c("Mixed Group 1", "Mixed Group 2",
                                          "Mixed Group 3", "Mixed Group 4",
                                          "Males-only Group", "Females-only Group"),
                     values = c(1, 2, 16, 17, 8, 4), drop = FALSE) +
  labs(x = "# chasing", 
       title = "Clumped resource condition \n Females",
       y = bquote(atop("Nectar intake ["~mL~h^-1*"]", "(Over last two nights)"))) +
  # geom_vline(xintercept = 0.003, linetype = 3) +
  # geom_hline(yintercept = 0.75, linetype = 3) +
  ylim(0, 2) + 
  # xlim(0, 0.013) +
  theme_serif() +
  theme(legend.text = element_text(size = 18, family = "serif"))

ggarrange(chases_males, chases_females, labels = c("A", "B"),
          font.label = list(size = 18, family = "serif"),
          common.legend = TRUE, legend = "right")
## first vs. second night

p_crit <- 0.1

cons_contrasts <- vis_summaries %>%
  filter(group_night < 3) %>% 
  group_by(group, sex, IdLabel) %>% 
  mutate(partition = case_when(
    group_night == 1 ~ "First night",
    group_night == 2 ~ "Second night",
    TRUE ~ NA_character_
  ),
  partition = factor(partition, levels = c("First night", "Second night")),
  phase = factor(phase, labels = c("Clumped resource\ncondition",
                                   "Distributed resource\ncondition")))

cons_contrasts <- cons_contrasts %>% 
  ungroup() %>%
  filter(!is.na(partition)) %>% 
  group_by(group, IdLabel, phase, partition, sex) %>%
  summarise(vol_hr = sum(vol_consumed)/sum(as.numeric(phase_dur))) %>% 
  left_join(statuses) %>% 
  arrange(partition, IdLabel)

stat.test <- cons_contrasts %>%
  group_by(status, phase) %>%
  t_test(vol_hr ~ partition) %>% 
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>% 
  add_xy_position(x = "partition", group = "status", fun = "max") %>% 
  filter(p.adj < p_crit) %>% 
  mutate(p.adj = ifelse(p.adj < 0.001, "<0.001", round(p.adj, 3)))

ggboxplot(cons_contrasts, x = "partition", y = "vol_hr",
          color = "status", line.color = "gray", line.size = 0.4,
          palette = c("dodgerblue3","grey50","purple"),
          add = "jitter") +
  facet_grid(. ~ phase) +
  # stat_pvalue_manual(
  #   stat.test,  label = "p.adj.signif", tip.length = 0.02, hide.ns = TRUE
  # ) +
  stat_pvalue_manual(
    stat.test,  label = "p.adj", tip.length = 0.02
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = "", y = bquote("Nectar intake ["~mL~h^-1*"]"), color = "") +
  theme_serif()


mean_intake_all <- vis_summaries %>% 
  filter(cond == "test", phase == 1) %>% 
  group_by(group_night, IdLabel) %>% 
  summarise(vol_hr = sum(vol_consumed)/sum(phase_dur)) %>% 
  group_by(IdLabel) %>% 
  summarise(mean_vol_hr = mean(vol_hr)) %>% 
  summarise(mean_intake = round(mean(mean_vol_hr), 2),
            sd = round(sd(mean_vol_hr), 2))


### overall Glicko (not only last two nights)

gl <- only_chases %>% 
  filter(cond == "test") %>% 
  select(-rev_night, -cond) %>% 
  nest_by(group) %>% 
  mutate(glicko = list(glicko(data))) %>% 
  select(-data) %>% 
  mutate(IdLabel = list(pluck(glicko, 1, "Player")),
         glicko_rating = list(pluck(glicko, 1, "Rating")))

rankings <- gl %>% 
  select(-glicko) %>% 
  unnest(c(IdLabel, glicko_rating)) %>% 
  mutate(glicko_rank = 1:n())

chase_summ_all <- chase_summ %>% 
  left_join(chased_n) %>%
  left_join(bat_info) %>% 
  left_join(rankings) %>% 
  replace_na(list(n_chased = 0)) %>% 
  mutate(prop_chased = n_chased / n_detections,
         phase = factor(phase), 
         sex = factor(sex))
glicko_all <- chase_summ_all %>% 
  # select(glicko_rank) %>% 
  rename(glicko_all = glicko_rank) %>% 
  ungroup() %>% 
  select(glicko_all)

chase_summ %>% 
  left_join(rankings) %>% 
  left_join(statuses) %>% 
  bind_cols(glicko_all) %>%
  ggplot(aes(glicko_all, glicko_rank, color = status)) +
  geom_point() +
  facet_wrap(~group) +
  geom_abline()

# show dom female on boxplot, why not same as other plot?

dom_female <- chase_nectar %>% 
  ungroup() %>% 
  mutate(status = factor(case_when(
    prop_chases > 0.003 & vol_hr > 0.75 ~ "Dominant",
    TRUE ~ "Subordinate"
  ), levels = c("Dominant", "Subordinate"))) %>%
  select(IdLabel, status, sex) %>% 
  filter(sex == "f", status == "Dominant")

dom_f_intake <- cons_contrasts %>% 
  filter(IdLabel == dom_female$IdLabel)

ggboxplot(cons_contrasts, x = "phase", y = "vol_hr",
          color = "status", line.color = "gray", line.size = 0.4,
          palette = c("dodgerblue3","grey50","purple"),
          add = "jitter") +
  facet_grid(. ~ partition) +
  geom_point(data = dom_f_intake, color = "red") +
  xlab("")

onlychoices <- allnights_n %>% 
  left_join(statuses) %>% # statuses calculation in Rmd file!
  ungroup() %>% 
  filter(cond == "test", !is.na(IdLabel), choice == TRUE) %>% 
  select(group_night, IdLabel, status, group, DateTime, loc, phase)

flower_distribution(onlychoices, "6m", phases = "both")

# not much gained by these alternative visualizations