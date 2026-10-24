
# Load R packages  --------------------------------------------------------

library(lme4)          
library(DHARMa)        
library(performance)   
library(emmeans)       
library(broom.mixed)   
library(dplyr)
library(ggplot2)
library(bbmle)         
library(ggnewscale)
library(devEMF)        
library(ggeffects)     
library(lattice)       


# Load data into R --------------------------------------------------------

setwd("C:/Users/samue/Desktop/Research/Honours/Chapter_2_Relatedness_paper/Kinference")

infile  <- "kinference_all.csv"
outdir  <- "C:/Users/samue/Desktop"

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

all_pairs2 <- read.csv(infile, check.names = FALSE)

# Inspect

head(all_pairs2)
glimpse(all_pairs2)

# Get metadata 

meta <- read.csv("metadata.csv")

glimpse(meta)

meta1 <- meta %>%
  dplyr::select(id, year_caught, Site_ID) %>%
  dplyr::rename(
    year_caught_1 = year_caught,
    site_ID_1 = Site_ID
  )

meta2 <- meta %>%
  dplyr::select(id, year_caught, Site_ID) %>%
  dplyr::rename(
    year_caught_2 = year_caught,
    site_ID_2 = Site_ID
  )

# Tidy things up by joining our individuals and getting things ready --------

all_pairs2 <- all_pairs2 %>%
  left_join(meta1, by = c("id_1" = "id")) %>%
  left_join(meta2, by = c("id_2" = "id"))

all_pairs2$pw_yr <- ifelse(all_pairs2$year_caught_1 == all_pairs2$year_caught_2, 1, 0)

all_pairs2$pw_loc <- ifelse(all_pairs2$site_ID_1 == all_pairs2$site_ID_2, 1, 0)

# Factors for pairwise indicators (expect 0/1 coding)
all_pairs2$pw_yr  <- as.factor(all_pairs2$pw_yr)
all_pairs2$pw_loc <- as.factor(all_pairs2$pw_loc)
all_pairs2$rel <- as.factor(all_pairs2$rel)

# Quick tables#pw_yr Quick tables
table(all_pairs2$rel,   all_pairs2$pw_yr)
table(all_pairs2$pw_loc, all_pairs2$pw_yr)

names(all_pairs2)[names(all_pairs2) == ""] <- "empty_col"

## make fsp and hsp just "related"

all_pairs2 <- all_pairs2 %>%
  mutate(related = as.factor(rel %in% c("full-sibs", "half-sibs")))

# Run the models: year, location, and null --------------------------------

## Model 1: year-only

year_model <- glmer(related ~ pw_yr + (1 | id_1) + (1 | id_2),
                    data   = all_pairs2,
                    family = binomial,
                    control = glmerControl(optimizer = "bobyqa",
                       optCtrl = list(maxfun = 2e5)))

summary(year_model)
exp(0.727)
broom::tidy(year_model)

## Model 2: location-only

loc_model <- glmer(
  related ~ pw_loc + (1 | id_1) + (1 | id_2),
  data   = all_pairs2,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 2e5)))

summary(loc_model)
exp(-3.29385)
broom::tidy(loc_model)

## Model 3: Null model 

null_model <- glmer(
  related ~ 1 + (1 | id_1) + (1 | id_2),
  data   = all_pairs2,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 2e5))
)

summary(null_model)
exp(14.032490)
broom::tidy(null_model)

# Look at model outputs ---------------------------------------------------

## Model selection (AICc table)

AICcmodavg::aictab(list(year_model, 
            loc_model, 
            null_model), 
            modnames = c( "Cohort","Location", "Null"))  # AICc, LL

# 3) Diagnostics ----
set.seed(42)
simres_year <- simulateResiduals(year_model, n = 1000)
plot(simres_year)

set.seed(42)
simres_loc  <- simulateResiduals(loc_model,  n = 1000)
plot(simres_loc)

performance::check_overdispersion(year_model)
performance::check_singularity(year_model)

## 4) Effect sizes (odds ratios) ----

year_OR <- broom.mixed::tidy(
  year_model, effects = "fixed", conf.int = TRUE, conf.method = "Wald"
) %>%
  mutate(
    OR      = exp(estimate),
    OR_low  = exp(conf.low),
    OR_high = exp(conf.high)
  ) %>%
  select(term, estimate, std.error, OR, OR_low, OR_high, p.value)

year_OR

loc_OR <- broom.mixed::tidy(
  loc_model, effects = "fixed", conf.int = TRUE, conf.method = "Wald"
) %>%
  mutate(
    OR      = exp(estimate),
    OR_low  = exp(conf.low),
    OR_high = exp(conf.high)
  ) %>%
  select(term, estimate, std.error, OR, OR_low, OR_high, p.value)

loc_OR

null_OR <- broom.mixed::tidy(
  null_model, effects = "fixed", conf.int = TRUE, conf.method = "Wald"
) %>%
  mutate(
    OR      = exp(estimate),
    OR_low  = exp(conf.low),
    OR_high = exp(conf.high)
  ) %>%
  select(term, estimate, std.error, OR, OR_low, OR_high, p.value)

null_OR

## 5) Marginal means (year) ----

emm_resp <- emmeans::emmeans(year_model, ~ pw_yr, type = "response")
emm_df   <- as.data.frame(emm_resp)  # cols: pw_yr, prob, SE, df, asymp.LCL, asymp.UCL
emm_df

## OR via contrasts on link scale
emm_link <- emmeans::emmeans(year_model, ~ pw_yr, type = "link")
ct <- pairs(emm_link, reverse = TRUE) %>% summary(infer = TRUE, adjust = "none")

ct_num <- ct %>%
  mutate(
    estimate = as.numeric(estimate),
    asymp.LCL = as.numeric(asymp.LCL),
    asymp.UCL = as.numeric(asymp.UCL)
  )

or_label <- sprintf(
  "OR (Same vs Different) = %.2f [%.2f, %.2f]",
  exp(ct_num$estimate[1]),
  exp(ct_num$asymp.LCL[1]),
  exp(ct_num$asymp.UCL[1])
)

or_label

## Predicted probabilities for each row (year model)

all_pairs2$predicted <- predict(year_model, type = "response")

all_pairs2$predicted <- predict(
  year_model,
  newdata = all_pairs2,
  type = "response",
  allow.new.levels = TRUE
)

## Nice labels for plotting (ensure 0/1 are present)
all_pairs2 <- data.frame(all_pairs2) %>%
  dplyr::mutate(
    pw_yr_f = factor(pw_yr,  levels = c("0","1"), labels = c("Different year","Same year")),
    pw_loc_f = factor(pw_loc, levels = c("0","1"), labels = c("Different location","Same location"))
  )

formula(year_model)

## 6) Plot: year model ----
year_model_plot <- ggplot() +
  geom_jitter(
    data = all_pairs2,
    aes(x = pw_yr, y = predicted, fill = pw_yr),
    width = 0.2, height = 0.025,
    alpha = 0.8, colour = "black", shape = 21, size = 4
  ) +
  scale_fill_manual(values = c("grey", "black")) +
  geom_errorbar(
    data = emm_df,
    aes(x = pw_yr, ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.1, linewidth = 1.25, color = "red1"
  ) +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = emm_df,
    aes(x = pw_yr, y = prob, fill = pw_yr),
    size = 6, shape = 21, colour = "black"
  ) +
  geom_line(
    data = emm_df,
    aes(x = as.numeric(pw_yr), y = prob, group = 1),
    color = "red", linewidth = 1.25
  ) +
  scale_x_discrete(drop = FALSE, labels = c("0" = "Different", "1" = "Same")) +
  scale_y_continuous(limits = c(0, 0.078)) +
  labs(x = "Capture Year", y = "Predicted P(Related)") +
  theme_bw() +
  scale_colour_manual(values = c("grey", "black")) +
  scale_fill_manual(values = c("red", "red")) +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

print(year_model_plot)

devEMF::emf(file.path(outdir, "GLMM_year_model.emf"), width = 10, height = 6.5)
print(year_model_plot)
dev.off()

#### Location model results & plot ####

## Effect sizes (location)
loc_tab <- broom.mixed::tidy(
  loc_model, effects = "fixed", conf.int = TRUE, conf.method = "Wald"
) %>%
  mutate(
    OR      = exp(estimate),
    OR_low  = exp(conf.low),
    OR_high = exp(conf.high)
  ) %>%
  select(term, estimate, std.error, OR, OR_low, OR_high, p.value)

loc_tab

## Marginal means (location)
loc_emm <- emmeans::emmeans(loc_model, ~ pw_loc, type = "response")
loc_emm_df <- as.data.frame(loc_emm)

## Contrasts on link scale
loc_emm_link <- emmeans::emmeans(loc_model, ~ pw_loc, type = "link")
loc_ct <- pairs(loc_emm_link, reverse = TRUE) %>% summary(infer = TRUE, adjust = "none")

loc_ct_num <- loc_ct %>%
  mutate(
    estimate  = as.numeric(estimate),
    asymp.LCL = as.numeric(asymp.LCL),
    asymp.UCL = as.numeric(asymp.UCL)
  )

loc_or_label <- sprintf(
  "OR (Same vs Different) = %.2f [%.2f, %.2f]",
  exp(loc_ct_num$estimate[1]),
  exp(loc_ct_num$asymp.LCL[1]),
  exp(loc_ct_num$asymp.UCL[1])
)

loc_or_label

## Row-level predictions (location)
all_pairs2$predicted_loc <- predict(loc_model, type = "response")

all_pairs2$predicted_loc <- predict(
  loc_model,
  newdata = all_pairs2,
  type = "response",
  allow.new.levels = TRUE
)

all_pairs2$pw_loc <- factor(all_pairs2$pw_loc,
                            levels = c(0, 1),
                            labels = c("Different", "Same"))

loc_emm_df$pw_loc <- factor(loc_emm_df$pw_loc,
                            levels = c(0, 1),
                            labels = c("Different", "Same"))

## Plot: location model
loc_model_plot <- ggplot() +
  geom_jitter(
    data = all_pairs2,
    aes(x = pw_loc, y = predicted_loc, fill = as.factor(pw_loc)),
    width = 0.2, height = 0.025,
    alpha = 0.8, colour = "black", shape = 21, size = 4
  ) +
  scale_fill_manual(values = c("grey", "black")) +
  geom_errorbar(
    data = loc_emm_df,
    aes(x = pw_loc, ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.1, linewidth = 1.25, color = "red1"
  ) +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = loc_emm_df,
    aes(x = pw_loc, y = prob, fill = as.factor(pw_loc)),
    size = 6, shape = 21, colour = "black"
  ) +
  geom_line(
    data = loc_emm_df,
    aes(x = as.numeric(pw_loc), y = prob, group = 1),
    color = "red", linewidth = 1.25
  ) +
  scale_x_discrete(drop = FALSE, labels = c("0" = "Different", "1" = "Same")) +
  scale_y_continuous(limits = c(0, 0.08)) +
  labs(x = "Capture Location", y = "Predicted P(Related)") +
  theme_bw() +
  scale_colour_manual(values = c("grey", "black")) +
  scale_fill_manual(values = c("red", "red")) +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

print(loc_model_plot)

devEMF::emf(file.path(outdir, "GLMM_location_model.emf"), width = 10, height = 6.5)
print(loc_model_plot)
dev.off()

## 7) Side-by-side figure ----
library(ggpubr)

glmm_combined <- ggarrange(
  year_model_plot, 
  loc_model_plot,
  labels = c("A", "B"),
  ncol   = 2, nrow = 1
)
plot(glmm_combined)

devEMF::emf(file.path(outdir, "GLMM_models.emf"), width = 10, height = 6)
print(glmm_combined)
dev.off()

## 8) Additional visual (year model, jitter + smoother) ----
glmm_plot <- ggplot(all_pairs2, aes(x = pw_yr, y = predicted, color = pw_yr)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 4) +
  stat_smooth(method = "glm", method.args = list(family = binomial), se = FALSE) +
  labs(x = "Same Cohort (0 = Different, 1 = Same)",
       y = "Predicted Probability of Kinship",
       color = "Same Cohort") +
  theme_bw() +
  scale_colour_manual(values = c("grey", "black")) +
  theme(legend.position = "top")

print(glmm_plot)

devEMF::emf(file.path(outdir, "kinship_glmm.emf"), width = 10, height = 8)
print(glmm_plot)
dev.off()

## 9) Interaction model (explicit) ----
glmer_model_interaction <- glmer(
  Relatives ~ pw_loc * pw_yr + (1 | id_1) + (1 | id_2),
  data = all_pairs2, 
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
summary(glmer_model_interaction)

## Compare all three (drop the undefined 'glmm_model2')
AICctab(glmm_model, year_model, loc_model, glmer_model_interaction)

## 10) ggpredict & odds figure ----
preds <- ggpredict(glmm_model, terms = c("pw_yr"),
                   bias_correction = TRUE, type = "fixed", back.transform = FALSE)
glimpse(preds)

preds$odds      <- exp(preds$predicted)
preds$odds_low  <- exp(preds$conf.low)
preds$odds_high <- exp(preds$conf.high)

glmm_plot2 <- ggplot(preds, aes(x = x, y = odds, color = group)) +
  geom_point(size = 4.5) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = odds_low, ymax = odds_high, fill = group),
              alpha = 0.2, color = NA) +
  scale_y_log10() +
  labs(
    x = "Within Cohort (0 = Different, 1 = Same)",
    y = "Odds of Being Relatives (log10 scale)",
    color = "Within Cohort",
    fill  = "Within Cohort"
  ) +
  scale_colour_manual(values = c("grey", "black")) +
  scale_fill_manual(values   = c("grey", "black")) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        legend.position = "top")

print(glmm_plot2)

devEMF::emf(file.path(outdir, "kinship_glmm_2.emf"), width = 10, height = 8)
print(glmm_plot2)
dev.off()

