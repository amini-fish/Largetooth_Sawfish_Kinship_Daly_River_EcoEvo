#### Simplified GLS Modelling ####

library(nlme)
library(corMLPE)
library(bbmle)
library(ggeffects)
library(tidyverse)
library(paletteer)
library(performance)
library(emmeans)
library(lme4)
library(MuMIn)
library(rcompanion)
library(ggsignif)
library(glmmTMB)
library(dplyr)
library(ggpubr)

## remove samples from 2011 and/or those from the main channel of the Daly River

data_3 <- data2 %>%
  dplyr::filter(!billabong_both %in% c("Daly River:Daly River", "Daly River:Wunullen", "Dalagarr:Daly River", "Daly River:Tyumalagun")) 

data_3$billabong_both <- droplevels(data_3$billabong_both)
data_3$year_caught_both <- droplevels(data_3$year_caught_both)

levels(data_3$billabong_both)
levels(data_3$year_caught_both)

## Rerun the sample gls models and re-examine

location_MLPE_model <- gls(log(relatedness+0.001) ~ billabong_both, 
            data = data_3, 
            correlation = corMLPE(form = ~ID1 + ID2), 
            method = "ML")

year_MLPE_model <- gls(log(relatedness+0.001) ~ year_caught_both, 
            data = data_3, 
            correlation = corMLPE(form = ~ID1 + ID2), 
            method = "ML")

null_MLPE_model <- gls(log(relatedness+0.001) ~ 1, 
              data = data_3, 
              correlation = corMLPE(form = ~ID1 + ID2), 
              method = "ML")

summary(location_MLPE_model)
summary(year_MLPE_model)
summary(null_MLPE_model)

anova(location_MLPE_model)
anova(year_MLPE_model)
anova(null_MLPE_model)

AICctab(location_MLPE_model, year_MLPE_model, null_MLPE_model)


table(data_3$year_caught_both)
table(data_3$billabong_both)

## General summary statistics #

aictab(list(location_MLPE_model, year_MLPE_model, null_MLPE_model), modnames = c("Location", "Cohort", "Null"))  # AICc, LL
sapply(list(location_MLPE_model, year_MLPE_model, null_MLPE_model), function(m) coef(m$modelStruct$corStruct, unconstrained = FALSE))  # rho
anova(m_1)$`p-value`[2]   # exact p for location, if you want more than "< 0.0001"

#### Model Assessments ####

loc_emm <- emmeans(location_MLPE_model, ~ billabong_both); loc_emm
pairs(loc_emm, simple = "each")

df_loc <- data.frame(regrid(loc_emm)); df_loc
range(df_loc$response)

summary(loc_emm, type = "response") 

year_emm <- emmeans(year_MLPE_model, ~ year_caught_both); year_emm
df <- data.frame(regrid(year_emm)); df
range(df$response)

year_emm_df <- as.data.frame(year_emm)
year_emm_df$pred_relatedness <- exp(year_emm_df$emmean) - 0.001
year_emm_df$lower_CL <- exp(year_emm_df$lower.CL) - 0.001
year_emm_df$upper_CL <- exp(year_emm_df$upper.CL) - 0.001

year_emm_df

emmeans(year_MLPE_model, ~ year_caught_both, 
        type = "response")


## 2) Predictions from the model (back-transformed), coerced to same factor order

gp <- ggpredict(location_MLPE_model, terms = "billabong_both", back.transform = TRUE) %>%
  as.data.frame()

gp

location_MLPE_model

# Estimated marginal means (on the log scale you modeled)
emm <- emmeans(location_MLPE_model, ~ billabong_both)

# Tukey-adjusted pairwise comparisons
pairs(emm, adjust = "tukey")
confint(emm, adjust = "tukey")

# 3. Convert to a data frame and back-transform to the relatedness scale
emm_df <- as.data.frame(emm)
emm_df$pred_relatedness <- exp(emm_df$emmean) - 0.001
emm_df$lower_CL <- exp(emm_df$lower.CL) - 0.001
emm_df$upper_CL <- exp(emm_df$upper.CL) - 0.001

emm_df

# Pairwise contrasts from your GLS
pair_contrasts <- pairs(emmeans(location_MLPE_model, ~ billabong_both), adjust = "tukey")
pair_df <- as.data.frame(pair_contrasts)

head(pair_df)

pair_df <- within(pair_df, {
  ratios        <- exp(estimate)
  lower_ratio   <- exp(estimate - 1.96 * SE)
  upper_ratio   <- exp(estimate + 1.96 * SE)
  change_percent <- (ratios - 1) * 100
})

pair_df

### Cap year 

# Estimated marginal means (on the log scale you modeled)
emm2 <- emmeans(year_MLPE_model, ~ year_caught_both)

# Tukey-adjusted pairwise comparisons

pairs(emm2, adjust = "tukey")
confint(emm2, adjust = "tukey")

# 3. Convert to a data frame and back-transform to the relatedness scale

emm_df2 <- as.data.frame(emm2)
emm_df2$pred_relatedness <- exp(emm_df2$emmean) - 0.001
emm_df2$lower_CL <- exp(emm_df2$lower.CL) - 0.001
emm_df2$upper_CL <- exp(emm_df2$upper.CL) - 0.001

emm_df2

# Pairwise contrasts from your GLS
pair_contrasts2 <- pairs(emmeans(year_MLPE_model, ~ year_caught_both), adjust = "tukey")
pair_df2 <- as.data.frame(pair_contrasts2)

head(pair_df2)

pair_df2 <- within(pair_df2, {
  ratios        <- exp(estimate)
  lower_ratio   <- exp(estimate - 1.96 * SE)
  upper_ratio   <- exp(estimate + 1.96 * SE)
  change_percent <- (ratios - 1) * 100
})

pair_df2 %>%
  filter(p.value < 0.1)


## 3) Styling (shared)
pal <- scale_fill_brewer(palette = "RdYlBu")
pal_col <- scale_color_brewer(palette = "RdYlBu")
base_theme <- theme_bw() + 
  theme(
    axis.title  = element_text(size = 14),
    axis.text   = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

## 4) Panel A — model effect (mean ± 95% CI) with raw points (optional jitter)
p1 <- ggplot(gp, aes(x, predicted, fill = x), 
             col = "black") +
  geom_pointrange(aes(ymin = conf.low, 
                      ymax = conf.high), 
                  size =1.5, 
                  shape = 21,
                  alpha = 0.9) +
  pal +
  labs(x = "Billabong Pair", y = "Predicted Rxy") +
  theme_bw()

print(p1)

###################
location_plot <- ggplot() +
  geom_jitter(data = data_3, 
              aes(x = billabong_both, y = relatedness, fill = billabong_both), 
              alpha = 0.4,
              size = 3,
              shape = 21, 
              colour = "black") +
  geom_errorbar(data = gp,
                aes(x = x, ymin = conf.low, ymax = conf.high),
                width = 0.1, linewidth = 1.25, color = "black") +
  geom_point(data = gp,
             aes(x = x, y = predicted, fill = x),
             size = 6, shape = 21, colour = "black") +
  
  pal +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, max(gp$conf.high, na.rm = TRUE) * 1.05)) +
  labs(x = "Pairwise capture location", y = "Predicted Rxy") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 14), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


print(location_plot)
###################

m2_2 <- gls(log(relatedness+0.001) ~ year_caught_both, 
            data = data_3, 
            correlation = corMLPE(form = ~ID1 + ID2), 
            method = "REML", 
            control = list(singular.ok = T))

gp_yr <- ggpredict(m2_2, terms = "year_caught_both", back.transform = TRUE) %>%
  as.data.frame() 

gp_yr <- arrange(gp_yr, desc(predicted))
gp_yr


## 4) Panel A — model effect (mean ± 95% CI) with raw points (optional jitter)
p2 <- ggplot(gp_yr, 
             aes( x = reorder(x, -predicted), predicted, fill = x), 
             col = "black") +
  geom_pointrange(aes(ymin = conf.low, 
                      ymax = conf.high), 
                  size =1.5, 
                  shape = 21,
                  alpha = 0.9) +
  pal +
  labs(x = "Year Pair", y = "Predicted Rxy") +
  base_theme 

print(p2)

## OR 

year_plot <- ggplot() +
  geom_jitter(data = data_3, 
              aes(x = year_caught_both, y = relatedness, fill = year_caught_both), 
              alpha = 0.5,
              size = 3,
              shape = 21, 
              colour = "black") +
  geom_errorbar(data = gp_yr,
                aes(x = x, ymin = conf.low, ymax = conf.high),
                width = 0.1, linewidth = 1.25, color = "black") +
  geom_point(data = gp_yr,
             aes(x = x, y = predicted, fill = x),
             size = 6, shape = 21, colour = "black") +
  
  pal +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, max(gp_yr$conf.high, na.rm = TRUE) * 1.05)) +
  labs(x = "Pairwise capture year", y = "Predicted Rxy") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 14), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

new_combined <- ggarrange(location_plot, 
                          year_plot, 
                          ncol = 1, 
                          nrow = 2, 
                          labels = c("A", "B"),
                          align = "v")

print(new_combined)

devEMF::emf("C:/Users/samue/Desktop/corMLPE_model_output.emf", width = 10, height = 10)  
print(new_combined)
dev.off()

anova(location_MLPE_model)
anova(m2_2)

summary(location_MLPE_model)
summary(year_MLPE_model)


anova(m1)
anova(m2)
summary(m1)
summary(m2)

# Model 1_1 residuals

plot(location_MLPE_model, resid(., type="normalized") ~ fitted(.), abline = c(0,0))


# Extract fitted and normalized residuals
df_res <- data.frame(
  fitted = fitted(location_MLPE_model),
  resid  = resid(location_MLPE_model, type = "normalized")
)

# ggplot version of: plot(location_MLPE_model, resid(., type="normalized") ~ fitted(.), abline = c(0,0))

anova(location_MLPE_model)
plot_res1a <- ggplot(df_res, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  labs(
    x = "Fitted values",
    y = "Normalised residuals",
    title = "Residuals vs Fitted values (GLS–MLPE model)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 14), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_res1b <- ggplot(data_3, aes(y=resid(location_MLPE_model, type="normalized"), 
                                 x= billabong_both)) + 
  geom_boxplot(fill = "grey") + 
  xlab("Pairwise capture location") + ylab("Normalised residuals") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 14), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

res1_combined <- ggarrange(plot_res1a, 
                           plot_res1b, 
                           ncol = 1, 
                           nrow = 2, 
                           labels = c("A", "B"),
                           align = "v")

print(res1_combined)

devEMF::emf("C:/Users/samue/Desktop/model1_resids.emf", width = 10, height = 8)  
print(res1_combined)
dev.off()

nagelkerke(location_MLPE_model) ## Note that the mdoel was refitted with ML not REML 

#### Model 2 

# Extract fitted and normalized residuals
df_res2 <- data.frame(
  fitted = fitted(year_MLPE_model),
  resid  = resid(year_MLPE_model, type = "normalized")
)

# ggplot version of: plot(location_MLPE_model, resid(., type="normalized") ~ fitted(.), abline = c(0,0))
plot_res2a <- ggplot(df_res2, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  labs(
    x = "Fitted values",
    y = "Normalised residuals",
    title = "Residuals vs Fitted values (GLS–MLPE model)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 14), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_res2b <- ggplot(data_3, aes(y=resid(year_MLPE_model, type="normalized"), 
                                 x= year_caught_both)) + 
  geom_boxplot(fill = "grey") + 
  xlab("Pairwise capture year") + ylab("Normalised residuals") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 14), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

res2_combined <- ggarrange(plot_res2a, 
                           plot_res2b, 
                           ncol = 1, 
                           nrow = 2, 
                           labels = c("A", "B"),
                           align = "v")

print(res2_combined)

devEMF::emf("C:/Users/samue/Desktop/model2_resids.emf", width = 10, height = 8)  
print(res2_combined)
dev.off()

