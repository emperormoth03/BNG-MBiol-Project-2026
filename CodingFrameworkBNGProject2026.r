############ AI Disclosure ----
#I used generative AI tools (ChatGPT) occasionally to assist with 
#debugging, as well as troubleshooting statistical and data
#visualisation workflows. All analytical decisions and 
#final code were reviewed and produced by myself. 
############ 1.1 Packages ---- 
library(readxl)       #To load in excel files
library(tidyverse)    #Contains ggplot2, dplyr, tidyr
library(vegan)        #Calculate Shannon diversity and Mantel
library(glmmTMB)      #GLMMs 
library(lmerTest)     #LMMs (Carabid Shannon diversity)
library(ggeffects)    #Generate model predictions
library(janitor)      #Clean data
library(DHARMa)       #Model diagnostics
library(performance)  #Model diagnostics 
library(flextable)    #Create flextables for diss
library(officer)      #Flextable captions
library(geosphere)    #Mantel test

############ 1.2 GLMM Diagnostics function ----
#Diagnostics function
check_glmm <- function(model) {
  
  #Simulated residuals
  sim_res <- simulateResiduals(model, n = 1000)
  
  #Basic plots
  plot(sim_res)                      #QQ + residual vs predicted
  testDispersion(sim_res)            #over/underdispersion
  testZeroInflation(sim_res)         #zero inflation
  testUniformity(sim_res)            #overall distribution
  
  #Residuals vs predictors
  plotResiduals(sim_res, form = model.frame(model)$d_c_sc)
  plotResiduals(sim_res, form = model.frame(model)$sampling_round)
  
  #Model performance checks
  print(check_overdispersion(model))
  print(check_singularity(model))
  print(check_collinearity(model))
  
  #AIC for quick reference
  print(AIC(model))
}


############ 1.3 GLMM Summary Table function ----
#Create appendix tables for GLMMs
appendix_glmm_flextable <- function(model, digits = 3, model_name = deparse(substitute(model))) {
  sm <- summary(model)
  
  #Extract fixed effects with failsafe
  coef_mat <- if (!is.null(sm$coefficients$cond)) {
    sm$coefficients$cond
  } else if (!is.null(sm$coefficients)) {
    sm$coefficients
  } else {
    stop("Could not find fixed-effect coefficients in model summary.")
  }
  
  fixed_df <- as.data.frame(coef_mat) %>%
    rownames_to_column("Term") %>%
    mutate(
      CI_low = Estimate - 1.96 * `Std. Error`,
      CI_high = Estimate + 1.96 * `Std. Error`,
      p_value = ifelse(`Pr(>|z|)` < 0.001, "<0.001", sprintf(paste0("%.", digits, "f"), `Pr(>|z|)`)),
      Estimate = round(Estimate, digits),
      `Std. Error` = round(`Std. Error`, digits),
      `z value` = round(`z value`, digits),
      CI_low = round(CI_low, digits),
      CI_high = round(CI_high, digits)
    ) %>%
    select(Term, Estimate, `Std. Error`, `z value`, p_value, CI_low, CI_high) %>%
    rename(`p value` = p_value, `95% CI low` = CI_low, `95% CI high` = CI_high)
  
  fixed_df <- fixed_df %>%
    mutate(
      Term = recode(Term,
                    "(Intercept)" = "Intercept",
                    "d_c_sc" = "Distinctiveness*Condition Score",
                    "sampling_roundTwo" = "Sampling round (2 vs 1)",
                    "sampling_roundThree" = "Sampling round (3 vs 1)",
                    "d_c_sc:sampling_roundTwo" = "Distinctiveness*Condition Score × round 2",
                    "d_c_sc:sampling_roundThree" = "Distinctiveness*Condition Score × round 3"))
  
  #Model summary info
  model_info <- tibble(
    Metric = c("Family", "Link", "N", "AIC", "BIC", "logLik"),
    Value = c(
      family(model)$family,
      family(model)$link,
      nobs(model),
      round(AIC(model), digits),
      round(BIC(model), digits),
      round(as.numeric(logLik(model)), digits)
    )
  )
  
  #Random effects summary if available
  re_df <- tryCatch({
    vc <- as.data.frame(VarCorr(model))
    if (nrow(vc) == 0) {
      NULL
    } else {
      vc %>%
        mutate(
          Effect = ifelse(is.na(var1) | var1 == "(Intercept)", paste0(grp, " (Intercept)"), paste0(grp, ": ", var1)),
          Variance = round(vcov, digits),
          `Std. Dev.` = round(sdcor, digits)
        ) %>%
        select(Effect, Variance, `Std. Dev.`)
    }
  }, error = function(e) NULL)
  
  ft_fixed <- flextable(fixed_df) %>%
    set_header_labels(
      Term = "Term",
      Estimate = "Estimate",
      `Std. Error` = "Std. Error",
      `z value` = "z-value",
      `p value` = "p-value",
      `95% CI low` = "95% CI low",
      `95% CI high` = "95% CI high"
    ) %>%
    theme_booktabs() %>%
    autofit() %>%
    bold(part = "header") %>%
    align(align = "center", part = "all") %>%
    align(j = "Term", align = "left", part = "all") %>%
    set_caption(as_paragraph(as_chunk(paste0("Fixed effects: ", model_name))))
  
  ft_info <- flextable(model_info) %>%
    theme_booktabs() %>%
    autofit() %>%
    bold(part = "header") %>%
    align(align = "center", part = "all") %>%
    align(j = "Metric", align = "left", part = "all") %>%
    set_caption(as_paragraph(as_chunk(paste0("Model summary: ", model_name))))
  
  ft_re <- if (!is.null(re_df)) {
    flextable(re_df) %>%
      theme_booktabs() %>%
      autofit() %>%
      bold(part = "header") %>%
      align(align = "center", part = "all") %>%
      align(j = "Effect", align = "left", part = "all") %>%
      set_caption(as_paragraph(as_chunk(paste0("Random effects: ", model_name))))
  } else {
    NULL
  }
  
  list(
    model_summary = ft_info,
    fixed_effects = ft_fixed,
    random_effects = ft_re
  )
}

#Create table, enter desired model 
tabs <- appendix_glmm_flextable(m_richness_hedge_fams_int)

tabs$model_summary
tabs$fixed_effects
tabs$random_effects


#Function below required for Shannon model as lmer used
tabs <- appendix_mixed_flextable(m_carab_shannon)

tabs$model_summary
tabs$fixed_effects

appendix_mixed_flextable <- function(model, digits = 3, model_name = deparse(substitute(model))) {
  sm <- summary(model)
  
  coef_obj <- sm$coefficients
  coef_mat <- if (is.matrix(coef_obj)) {
    coef_obj
  } else if (is.list(coef_obj) && !is.null(coef_obj$cond) && is.matrix(coef_obj$cond)) {
    coef_obj$cond
  } else {
    stop("Could not find a coefficient matrix in the model summary.")
  }
  
  coef_df <- as.data.frame(coef_mat) %>%
    rownames_to_column("Term")
  
  stat_col <- intersect(c("t value", "z value"), names(coef_df))[1]
  p_col <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), names(coef_df))[1]
  
  if (is.na(stat_col) || is.na(p_col)) {
    stop("Could not identify statistic or p-value columns.")
  }
  
  fixed_df <- coef_df %>%
    mutate(
      CI_low = Estimate - 1.96 * `Std. Error`,
      CI_high = Estimate + 1.96 * `Std. Error`
    ) %>%
    mutate(
      Estimate = round(Estimate, digits),
      `Std. Error` = round(`Std. Error`, digits),
      !!stat_col := round(.data[[stat_col]], digits),
      CI_low = round(CI_low, digits),
      CI_high = round(CI_high, digits),
      `p-value` = ifelse(.data[[p_col]] < 0.001, "<0.001", sprintf(paste0("%.", digits, "f"), .data[[p_col]]))
    ) %>%
    select(Term, Estimate, `Std. Error`, all_of(stat_col), `p-value`, CI_low, CI_high) %>%
    rename(
      Statistic = all_of(stat_col),
      `95% CI low` = CI_low,
      `95% CI high` = CI_high
    ) %>%
    mutate(
      Term = recode(Term,
                    "(Intercept)" = "Intercept",
                    "d_c_sc" = "Distinctiveness*Condition Score",
                    "sampling_roundTwo" = "Sampling round (2 vs 1)",
                    "sampling_roundThree" = "Sampling round (3 vs 1)"
      )
    )
  
  model_info <- tibble(
    Metric = c("Family", "Link", "N", "AIC", "BIC", "logLik"),
    Value = c(
      family(model)$family,
      family(model)$link,
      nobs(model),
      round(AIC(model), digits),
      round(BIC(model), digits),
      round(as.numeric(logLik(model)), digits)
    )
  )
  
  ft_fixed <- flextable(fixed_df) %>%
    set_header_labels(
      Term = "Term",
      Estimate = "Estimate",
      `Std. Error` = "Std. Error",
      Statistic = ifelse(stat_col == "t value", "t-value", "z-value"),
      `p-value` = "p-value",
      `95% CI low` = "95% CI low",
      `95% CI high` = "95% CI high"
    ) %>%
    theme_booktabs() %>%
    autofit() %>%
    bold(part = "header") %>%
    align(align = "center", part = "all") %>%
    align(j = "Term", align = "left", part = "all") %>%
    set_caption(as_paragraph(as_chunk(paste0("Fixed effects: ", model_name))))
  
  ft_info <- flextable(model_info) %>%
    theme_booktabs() %>%
    autofit() %>%
    bold(part = "header") %>%
    align(align = "center", part = "all") %>%
    align(j = "Metric", align = "left", part = "all") %>%
    set_caption(as_paragraph(as_chunk(paste0("Model summary: ", model_name))))
  
  list(
    model_summary = ft_info,
    fixed_effects = ft_fixed
  )
}
############ 1.4 Hedgerows Flextable ----
#Read Hedgerows sheet
hedgerows <- read_excel("MasterDataNew.xlsx", sheet = "Hedgerows")

#Select and rename columns
pitfall_table <- hedgerows %>%
  select(
    Hedgerow,
    Latitude,
    Longitude,
    `Habitat Type`,
    Distinctiveness,
    Condition,
    `D*C`) %>%
  rename(
    `Hedgerow` = Hedgerow,
    `Habitat type` = `Habitat Type`,
    `Distinctiveness score` = Distinctiveness,
    `Condition score` = Condition,
    `Distinctiveness*Condition` = `D*C`)

#Create flextable
ft_pitfalls <- flextable(pitfall_table) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = c("Hedgerow", "Habitat type"), align = "left", part = "all") %>%
  autofit() %>%
  width(j = c("Habitat type"), width = 5) %>%
  width(j = c("Hedgerow"), width = 0.9) %>%
  width(j = c("Latitude", "Longitude", "Distinctiveness score", "Condition score", "Distinctiveness*Condition"), width = 1.1) %>%
  set_caption("Pitfall locations, habitat types, and score values")

ft_pitfalls

############ 2.1 Load and Prepare Datasets ---- 
PitfallOrders <- read_excel(
  "MasterDataNew.xlsx",
  sheet = "PitfallOrders") %>%
  clean_names()

PitfallFamilies <- read_excel(
  "MasterDataNew.xlsx",
  sheet = "PitfallFamilies") %>%
  clean_names()

#NOTE: PitfallOrders and PitfallFamilies exclude ants, this can be added back in manually if desired

HedgeOrders <- read_excel(
  "MasterDataNew.xlsx",
  sheet = "HedgeOrders") %>%
  clean_names()

HedgeFamilies <- read_excel(
  "MasterDataNew.xlsx",
  sheet = "HedgeFamilies") %>%
  clean_names()

#Identify non-taxon columns
meta_cols <- c(
  "pitfall",
  "sampling_round",
  "habitat_type",
  "distinctiveness",
  "condition",
  "d_c",
  "total_abundance",
  "taxon_richness",
  "adjacent_land_use",
  "width",
  "height", 
  "aspect", 
  "direction",
  "latitude",
  "longitude")

#Derive land_use_bin and clean data
prep_pitfall <- function(df, meta_cols) {
  
  taxon_cols <- setdiff(names(df), meta_cols)
  
  df_clean <- df %>%
    mutate(
      pitfall = factor(pitfall),
      sampling_round   = factor(sampling_round, levels = c("One", "Two", "Three"))) %>%
    mutate(
      across(all_of(taxon_cols), ~ as.numeric(as.character(.))),
      across(all_of(taxon_cols), ~ replace_na(., 0)))
  
  df_clean %>%
    rowwise() %>%
    ungroup() %>%
    filter(
      pitfall != "MPH10",  #Ant colony pitfall to remove
      total_abundance > 0) %>%
    mutate(
      land_use_bin = case_when(
        adjacent_land_use %in% c("Wheat Field","Barley Field","Oat Field") ~ "agricultural_field",
        adjacent_land_use %in% c("Species-rich Meadow","Wild Grassland","Mowed Grass Field","Meadow then Barley","Path Verge") ~ "non_agricultural_grasslands",
        TRUE ~ NA_character_),
      land_use_bin = factor(land_use_bin, levels = c("agricultural_field","non_agricultural_grasslands")),
      d_c_sc = as.numeric(scale(d_c)))
}

PitfallOrdersReady <- prep_pitfall(PitfallOrders, meta_cols)
PitfallFamiliesReady  <- prep_pitfall(PitfallFamilies, meta_cols)
HedgeOrdersReady <- prep_pitfall(HedgeOrders, meta_cols)
HedgeFamiliesReady  <- prep_pitfall(HedgeFamilies, meta_cols)

#Stage requires carabid_cols to be read first (found in 2.5 Abundance Plots)
CarabidData <- read_excel(
  "MasterDataNew.xlsx",
  sheet = "Carabids") %>%
  clean_names()

#Remove lost pitfalls
  CarabidData <- CarabidData %>%
  filter(
    !(pitfall == "MPH10"),
    !(pitfall == "MPH23" & sampling_round == "One"),
    !(pitfall == "MPH9"  & sampling_round == "Two"),
    !(pitfall == "MPH20" & sampling_round == "Three"))

#Identify carabid species columns
#Carabid Shannon diversity 
CarabidData <- CarabidData %>%
  rowwise() %>%
  mutate(
    carab_abundance = sum(c_across(all_of(carabid_cols))),
    carab_shannon   = ifelse(
      carab_abundance > 0,
      diversity(c_across(all_of(carabid_cols)), index = "shannon"),
      NA_real_)) %>%
  ungroup()

#Scale BNG score
CarabidData <- CarabidData %>%
  mutate(
    d_c_sc = as.numeric(scale(d_c)))

CarabidDataReady <- CarabidData

#Create land_use_bin
CarabidDataReady <- CarabidDataReady %>%
  mutate(
    land_use_bin = case_when(
      adjacent_land_use %in% c("Wheat Field","Barley Field","Oat Field") ~ "agricultural_field",
      adjacent_land_use %in% c("Species-rich Meadow","Wild Grassland","Mowed Grass Field","Meadow then Barley","Path Verge") ~ "non_agricultural_grasslands",
      TRUE ~ NA_character_
    ),
    land_use_bin = factor(land_use_bin, levels = c("agricultural_field","non_agricultural_grasslands"))
  )

############ 2.2 Pitfall Abundance GLMMs ---- 
#GLMMs for abundance - order and family level (order level used in diss, as contains all individuals)
m_abund_orders <- glmmTMB(
  total_abundance ~ d_c_sc + sampling_round + (1 | pitfall),
  family = nbinom2,
  data = PitfallOrdersReady)

m_abund_orders_int <- glmmTMB(
  total_abundance ~ d_c_sc * sampling_round + (1 | pitfall),
  family = nbinom2,
  data = PitfallOrdersReady)

m_abund_fams <- glmmTMB(
  total_abundance ~ d_c_sc + sampling_round + (1 | pitfall),
  family = nbinom2,
  data = PitfallFamiliesReady)

m_abund_fams_int <- glmmTMB(
  total_abundance ~ d_c_sc * sampling_round + (1 | pitfall),
  family = nbinom2,
  data = PitfallFamiliesReady)

summary(m_abund_orders)
summary(m_abund_orders_int)
summary(m_abund_fams)
summary(m_abund_fams_int)

#Interaction model passes diagnostics but additive doesn't
check_glmm(m_abund_orders_int)

anova(m_abund_orders, m_abund_orders_int, test = "Chisq")

#Interaction model has lower AIC score for order (all counts)
AIC(m_abund_orders)
AIC(m_abund_orders_int)
AIC(m_abund_fams)
AIC(m_abund_fams_int)

#Sensitivity Analysis - adjacent land use
m_abund_orders_int2 <- glmmTMB(
  total_abundanceance ~ d_c_sc*sampling_round + land_use_bin + (1 | pitfall),
  family = nbinom2,
  data = PitfallOrdersReady)

#Incorporating adjacent land has no sig dif with ANOVA
summary(m_abund_orders_int2)
AIC(m_abund_orders_int,m_abund_orders_int2)
anova(m_abund_orders_int,m_abund_orders_int2)

#Function to extract coefficients and calculate IRR
extract_irr <- function(model, conf.level = 0.95, exclude_intercept = TRUE) {
  sm <- summary(model)
  
  #Fixed-effect table from glmmTMB summary
  coefs <- as.data.frame(sm$coefficients$cond) %>%
    rownames_to_column("term") %>%
    as_tibble()
  
  #Standardise column names
  names(coefs) <- sub("^Estimate$", "beta", names(coefs))
  names(coefs) <- sub("^Std\\. Error$", "se", names(coefs))
  names(coefs) <- sub("^Pr\\(>\\|z\\|\\)$", "p.value", names(coefs))
  
  z <- qnorm(1 - (1 - conf.level) / 2)
  
  out <- coefs %>%
    mutate(
      IRR = exp(beta),
      conf.low = exp(beta - z * se),
      conf.high = exp(beta + z * se)) %>%
    select(term, IRR, conf.low, conf.high, beta, se, p.value)
  
  if (exclude_intercept) {
    out <- out %>% filter(term != "(Intercept)")
  }
  
  out
}
extract_irr(m_abund_orders_int)





#PLOTS 
#Settings
point_alpha <- 0.8
point_size  <- 2.2
jitter_width <- 0.3
plot_width <- 8
plot_height <- 5
dpi <- 600

#Helper: make prediction data (raw d_c) and predict response-scale + CI from a glmmTMB model
make_pred_df_rawx <- function(model, data, raw_x = "d_c", scaled_x = "d_c_sc", group_var = "sampling_round",
                              n = 200, re.form = NA, level = 0.95) {
  #Ensure group factor levels
  data <- data %>% mutate(!!group_var := as.factor(.data[[group_var]]))
  #Mean/sd used to create scaled predictor in the data (assumes d_c_sc = scale(d_c))
  mean_x <- mean(data[[raw_x]], na.rm = TRUE)
  sd_x   <- sd(data[[raw_x]], na.rm = TRUE)
  #Prediction grid on raw scale
  newdat <- expand.grid(
    d_c = seq(min(data[[raw_x]], na.rm = TRUE),
              max(data[[raw_x]], na.rm = TRUE),
              length.out = n),
    sampling_round = levels(as.factor(data[[group_var]])))
  #Create scaled predictor matching model input
  newdat[[scaled_x]] <- (newdat$d_c - mean_x) / sd_x
  #Ensure sampling_round is same class
  newdat[[group_var]] <- factor(newdat[[group_var]], levels = levels(as.factor(data[[group_var]])))
  #Predict on link scale (for glmmTMB with log link) and attempt to get se
  pred_link <- tryCatch(predict(model, newdata = newdat, type = "link", se.fit = TRUE, re.form = re.form),
                        error = function(e) predict(model, newdata = newdat, type = "link", re.form = re.form))
  #Extract fit and se.fit 
  if(is.list(pred_link) && !is.null(pred_link$fit)) {
    fit_link <- as.numeric(pred_link$fit)
    se_link  <- if(!is.null(pred_link$se.fit)) as.numeric(pred_link$se.fit) else NA_real_
  } else if(is.numeric(pred_link)) {
    #no se provided
    fit_link <- as.numeric(pred_link)
    se_link  <- rep(NA_real_, length(fit_link))
  } else {
    stop("Unexpected predict() return structure. Inspect predict output.")
  }
  #Back-transform to response scale (log link -> exp)
  z <- qnorm(1 - (1 - level) / 2)
  newdat$predicted <- exp(fit_link)
  newdat$conf.low  <- ifelse(is.na(se_link), NA, exp(fit_link - z * se_link))
  newdat$conf.high <- ifelse(is.na(se_link), NA, exp(fit_link + z * se_link))
  #rename raw column to standard name for plotting
  newdat <- newdat %>% rename(!!raw_x := d_c)
  return(newdat)
}

#Orders plot (raw d_c on x axis)
pred_orders_df <- make_pred_df_rawx(m_abund_orders_int, PitfallOrdersReady, raw_x = "d_c", scaled_x = "d_c_sc", group_var = "sampling_round")

#Determine max D×C for consistent scaling
max_dc_orders <- max(PitfallOrdersReady$d_c, na.rm = TRUE)
max_dc_fams   <- max(PitfallFamiliesReady$d_c, na.rm = TRUE)

max_dc <- max(max_dc_orders, max_dc_fams)

#Breaks in multiples of 2
dc_breaks <- seq(0, ceiling(max_dc / 2) * 2, by = 2)

sampling_roundcols <- c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")

pitfall_abundance_plot <- ggplot() +
  geom_jitter(data = PitfallOrdersReady, aes(x = d_c, y = total_abundance, color = sampling_round, shape = sampling_round),
              width = jitter_width, alpha = point_alpha, size = point_size) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_line(data = pred_orders_df, aes(x = d_c, y = predicted, colour = sampling_round, linetype = sampling_round), linewidth = 1) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_ribbon(data = pred_orders_df, aes(x = d_c, ymin = conf.low, ymax = conf.high, fill = sampling_round), alpha = 0.16, colour = NA) +
  scale_fill_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  #  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(limits = c(3, 25), breaks = dc_breaks, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 260), expand = c(0, 0)) +
  labs(x = "Distinctiveness*Condition Score",
       y = "Invertebrate Abundance",
#       title = "Pitfalls Invertebrate Abundance vs Distinctiveness*Condition Score",
       colour = "sampling_round", fill = "sampling_round", shape = "sampling_round", linetype = "sampling_round") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(legend.key.height = unit(0.5, 'cm'), #Change legend key height
        legend.key.width = unit(1.4, 'cm'), #Change legend key width
        legend.title = element_text(size=12), #Change legend title font size
        legend.text = element_text(size=10), #Change legend text font size
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))
pitfall_abundance_plot

ggsave(
  "Figure_Abundance_vs_BNG.png",
  plot = pitfall_abundance_plot,
  width = 6.5,
  height = 5,
  dpi = 600)


############ 2.3 Hedge Beating Abundance GLMMs ----
#GLMMs for abundance - orders and families (order level used in diss as contains all counts)
m_abund_orders_hedge <- glmmTMB(
  total_abundance ~ d_c_sc + sampling_round + (1 | pitfall),
  family = nbinom2,
  data = HedgeOrdersReady)

m_abund_orders_hedge_int <- glmmTMB(
  total_abundance ~ d_c_sc*sampling_round + (1 | pitfall),
  family = nbinom2,
  data = HedgeOrdersReady)

m_abund_fams_hedge <- glmmTMB(
  total_abundance ~ d_c_sc + sampling_round + (1 | pitfall),
  family = nbinom2,
  data = HedgeFamiliesReady)

m_abund_fams_hedge_int <- glmmTMB(
  total_abundance ~ d_c_sc*sampling_round + (1 | pitfall),
  family = nbinom2,
  data = HedgeFamiliesReady)

summary(m_abund_orders_hedge)
summary(m_abund_fams_hedge)
summary(m_abund_orders_hedge_int)
summary(m_abund_fams_hedge_int)

#Order level model for reporting passes all tests
check_glmm(m_abund_orders_hedge)

anova(m_abund_orders_hedge_int, m_abund_orders_hedge, test = "Chisq")

#AIC scores very similar, slightly lower for the additive models
AIC(m_abund_fams_hedge)
AIC(m_abund_fams_hedge_int)
AIC(m_abund_orders_hedge)
AIC(m_abund_orders_hedge_int)

extract_irr(m_abund_orders_hedge)

#Sensitivity Analysis - adjacent land use
m_abund_orders_hedge_int2 <- glmmTMB(
  total_abundance ~ d_c_sc*sampling_round + land_use_bin + (1 | pitfall),
  family = nbinom2,
  data = HedgeOrdersReady)

#Incorporating adjacent land use doesn't decrease AIC, no sig dif with ANOVA
summary(m_abund_orders_hedge_int2)
AIC(m_abund_orders_hedge_int,m_abund_orders_hedge_int2)
anova(m_abund_orders_hedge_int,m_abund_orders_hedge_int2)


#Plotting
pred_orders_df <- make_pred_df_rawx(m_abund_orders_hedge, HedgeOrdersReady, raw_x = "d_c", scaled_x = "d_c_sc", group_var = "sampling_round")

#determine max D×C for consistent scaling
max_dc <- max(HedgeOrdersReady$d_c, na.rm = TRUE)

#breaks in multiples of 2
dc_breaks <- seq(0, ceiling(max_dc / 2) * 2, by = 2)

hedge_abundance_plot <- ggplot() +
  geom_jitter(data = HedgeOrdersReady, aes(x = d_c, y = total_abundance, colour = sampling_round, shape = sampling_round),
              width = jitter_width, alpha = point_alpha, size = point_size) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_line(data = pred_orders_df, aes(x = d_c, y = predicted, colour = sampling_round, linetype = sampling_round), linewidth = 1) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_ribbon(data = pred_orders_df, aes(x = d_c, ymin = conf.low, ymax = conf.high, fill = sampling_round), alpha = 0.16, colour = NA) +
  scale_fill_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  #  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(limits = c(3, 25), breaks = dc_breaks, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
  labs(x = "Distinctiveness*Condition Score",
       y = "Invertebrate Abundance",
#       title = "Hedge Beating Invertebrate Abundance vs Distinctiveness*Condition Score",
       colour = "sampling_round", fill = "sampling_round", shape = "sampling_round", linetype = "sampling_round") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(1.4, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10), #change legend text font size
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))
hedge_abundance_plot 

ggsave(
  "Figure_Hedge_Abundance_vs_BNG.png",
  plot = hedge_abundance_plot,
  width = 6.5,
  height = 5,
  dpi = 600)

############ 2.4 Carabid Abundance GLMMs ----
#Abundance (NB GLMM)
m_carab_abund <- glmmTMB(
  total_abundance ~ d_c_sc + sampling_round + (1 | pitfall),
  family = nbinom2,
  data   = CarabidDataReady)

m_carab_abund2 <- glmmTMB(
  total_abundance ~ d_c_sc + sampling_round + land_use_bin + (1 | pitfall),
  family = nbinom2,
  data   = CarabidDataReady)

m_carab_abund_int <- glmmTMB(
  total_abundance ~ d_c_sc * sampling_round + (1 | pitfall),
  family = nbinom2,
  data   = CarabidDataReady)

#Both models pass diagnostic checks
summary(m_carab_abund)
summary(m_carab_abund_int)
extract_irr(m_carab_abund)

check_glmm(m_carab_abund)

#Additive has lower AIC score (528 compared to 531 for interaction)
AIC(m_carab_abund)
AIC(m_carab_abund_int)

#Sensitivity Analysis - adjacent land use
m_carab_abund2 <- glmmTMB(
  total_abundance ~ d_c_sc + sampling_round + land_use_bin + (1 | pitfall),
  family = nbinom2,
  data = CarabidDataReady)

#Incorporating adjacent land use decreases AIC, no sig dif with ANOVA
AIC(m_carab_abund,m_carab_abund2)
anova(m_carab_abund,m_carab_abund2)



#Recover scaling info for d_c (raw) from the same data used to fit the model
dc_mean <- mean(CarabidDataReady$d_c, na.rm = TRUE)
dc_sd   <- sd(CarabidDataReady$d_c, na.rm = TRUE)

#Determine breaks (multiples of 2) up to a sensible max
max_dc <- max(CarabidDataReady$d_c, na.rm = TRUE)
dc_breaks <- seq(0, ceiling(max_dc / 2) * 2, by = 2)

#Get model predictions across the observed range of d_c_sc and by sampling_round
pred <- ggpredict(m_carab_abund, terms = c("d_c_sc [all]", "sampling_round"))
pred_df <- as.data.frame(pred) %>%
  rename(d_c_sc = x, sampling_round = group) %>%
  mutate(d_c = d_c_sc * dc_sd + dc_mean)

#Original points (use the same dataframe used to fit the model)
plot_df <- CarabidDataReady

#Choose x/y limits
x_min <- min(plot_df$d_c, na.rm = TRUE); x_max <- max(plot_df$d_c, na.rm = TRUE)
y_min <- min(plot_df$taxon_richness, na.rm = TRUE); y_max <- max(plot_df$taxon_richness, na.rm = TRUE)

carabid_abundance_plot <- ggplot() +
  geom_jitter(data = CarabidDataReady, aes(x = d_c, y = total_abundance, color = sampling_round, shape = sampling_round),
              width = jitter_width, alpha = point_alpha, size = point_size) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_line(data = pred_df, aes(x = d_c, y = predicted, colour = sampling_round, linetype = sampling_round), linewidth = 1) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_ribbon(data = pred_df, aes(x = d_c, ymin = conf.low, ymax = conf.high, fill = sampling_round), alpha = 0.2, colour = NA) +
  scale_fill_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  scale_shape_manual(values = c("One" = 16, "Two" = 17, "Three" = 15)) +
  scale_x_continuous(limits = c(3,25), breaks = dc_breaks, expand = c(0,0)) +
  scale_y_continuous(limits = c(-1, 30), expand = c(0,0)) +
  labs(x = "Distinctiveness*Condition Score",
       y = "Carabid Abundance",
       #title = "Carabids - Total Abundance vs Distinctiveness*Condition Score",
       colour = "sampling_round", shape = "sampling_round",fill = "sampling_round", linetype = "sampling_round") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(1.4, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10), #change legend text font size
        axis.title = element_text(size = 15),  # axis labels (x and y titles)
        axis.text  = element_text(size = 12)) 
carabid_abundance_plot

ggsave(
  "Figure_Carab_Abundance_vs_BNG.png",
  plot = carabid_abundance_plot,
  width = 6.5,
  height = 5,
  dpi = 600)

############ 2.5 Abundance Tables ----
#Function to make taxon abundance bar plots
make_rank_bar_plot <- function(data, taxon_cols, title_text) {

  #Summary (top 15 for barplot)
  summary_df <- data %>%
    mutate(total = rowSums(across(all_of(taxon_cols)), na.rm = TRUE)) %>%
    pivot_longer(all_of(taxon_cols), names_to = "taxon", values_to = "count") %>%
    mutate(taxon = tools::toTitleCase(gsub("_", " ", taxon))) %>%
    #mutate(
    #  taxon = gsub("_", " ", taxon),             #Un-hash to make carabid plot
    #  taxon = stringr::str_to_sentence(taxon),   #Genus capitalised and everything italicised
    #  taxon = sub(
    #    "^(\\w+)\\s+(\\w+)",
    #    "\\1 \\L\\2",
    #    taxon,
    #    perl = TRUE
    #  )
    #) %>%
    group_by(taxon) %>%
    summarise(
      mean_pct = mean(count / total * 100, na.rm = TRUE),
      mean_count = mean(count, na.rm = TRUE),
      freq_pct   = mean(count > 0, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    arrange(desc(mean_pct))
  
  top_df <- summary_df %>% slice_head(n = 15)
  
  #Rank-abundance (ALL taxa)
  rank_df <- summary_df %>%
    mutate(rank = row_number()) %>% 
    slice_head(n = 30)
  
  #Rank-abundance plot - not included in diss
  p_rank <- ggplot(rank_df, aes(x = rank, y = mean_count)) +
    geom_line() +
    geom_point() +
    labs(
      x = "Taxon rank",
      y = "Mean % abundance",
      title = paste0(title_text, " – Rank-abundance")) +
    theme_bw()
  
  #Barplot (top 15)
  p_bar <- ggplot(top_df, aes(x = mean_count, y = reorder(taxon, mean_count))) +
    geom_col(fill = "#59A14F") + 
    labs(
      x = "Mean count per sample",
      y = NULL,
      title = paste0(title_text, " – Top 15 taxa")) +
    theme_bw()
  
  #Combine if desired
  #p_rank / p_bar
  p_bar +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 16))
    #axis.text.y = element_text(face = "italic"))
}

taxon_cols_pitfall_orders <- c(
  "acari",
  "araneae",
  "chilopoda",
  "clitellata",
  "dermaptera",
  "gastropoda",
  "isopoda",
  "julida",
  "lepidoptera",
  "opiliones",
  "polydesmida",
  "pseudoscorpions",
  "psocoptera",
  "entomobryomorpha",
  "symphypleona",
  "hemiptera",
  "hymenoptera",
  "diptera",
  "coleoptera")

taxon_cols_pitfall_families <- c(
  "byrrhidae",
  "chrysomelidae",
  "coccinellidae",
  "corylophidae",
  "cryptophagidae",
  "curculionidae",
  "latridiidae",
  "leiodidae",
  "lygaeidae",
  "scarabaeoidea",
  "silphidae",
  "staphylinidae",
  "tenebrionidae",
  "carabidae",
  "aphidae",
  "berytidae",
  "cicadellidae",
  "cydnidae",
  "hemerobiidae",
  "issidae",
  "lygaeidae",
  "nabidae",
  "thyreocoridae",
  "tingidae",
  "apidae",
  "vespidae",
  "scelionidae",
  "diapriidae",
  "anthomyiidae",
  "calliphoridae",
  "chironomidae",
  "empididae",
  "hybotidae",
  "phoridae",
  "sarcophagidae",
  "sepsidae",
  "sciaridae",
  "sphaeroceridae",
  "forficulidae",
  "julidae",
  "polydesmidae")

taxon_cols_hedge_orders <- c(
  "acari",
  "araneae",
  "dermaptera",
  "gastropoda",
  "isopoda",
  "lepidoptera",
  "myriapoda",
  "opiliones",
  "psocoptera",
  "tettigoniidae",
  "entomobryomorpha",
  "hemiptera",
  "hymenoptera",
  "coleoptera",
  "diptera")

taxon_cols_hedge_families <- c(
  "aphidae",
  "anthocoridae",
  "berytidae",
  "cercopoidae",
  "cicadellidae",
  "coreidae",
  "pentatomidae",
  "hemerobiidae",
  "heterogastridae",
  "miridae",
  "nabidae",
  "psyllidae",
  "reduviidae",
  "tingidae",
  "formicidae",
  "parasitica",
  "symphyta",
  "anthicidae",
  "buprestidae",
  "cantharidae",
  "carabidae",
  "coccinellidae",
  "curculionidae",
  "chrysomelidae",
  "cryptophagidae",
  "phalacridae",
  "ptinidae",
  "latridiidae",
  "leiodidae",
  "melyridae",
  "mordellidae",
  "nitidulidae",
  "staphylinidae",
  "tenebrionidae",
  "anthomyiidae",
  "calliphoridae",
  "chironomidae",
  "chloropidae",
  "cecidomyiidae",
  "drosophilidae",
  "ephydridae",
  "hybotidae",
  "lauxaniidae",
  "lonchopteridae",
  "simuliidae",
  "stratiomyidae",
  "tachinidae",
  "tipulidae",
  "dermaptera",
  "tettigoniidae")

carabid_cols <- c(
  "amara_bifrons",
  "amara_communis",
  "amara_familiaris",
  "amara_ovata",
  "amara_similata",
  "anchomenus_dorsalis",
  "badister_bullatus",
  "bembidion_guttula",
  "bembidion_lampros",
  "bembidion_quadrimaculatum",
  "brachinus_crepitans",
  "calathus_fuscipes",
  "calathus_rotundicollis",
  "clivinia_fossor",
  "curtonotus_aulicus",
  "harpalus_affinis",
  "harpalus_latus",
  "harpalus_rufipes",
  "harpalus_tardus",
  "microlestes_minutulus",
  "nebria_brevicollis",
  "notiophilus_biguttatus",
  "notiophilus_germinyi",
  "notiophilus_palustris",
  "notiophilus_rufipes",
  "ophonus_rufibarbis",
  "oxypselaphus_obscurus",
  "poecilus_cupreus",
  "poecilus_versicolor",
  "pterostichus_madidus",
  "pterostichus_melanarius",
  "pterostichus_niger",
  "pterostichus_strenuus",
  "syntomus_foveatus",
  "syntomus_truncatellus")

pordbar=make_rank_bar_plot(PitfallOrdersReady, taxon_cols_pitfall_orders, "Pitfall assemblage (order level)")
pfambar=make_rank_bar_plot(PitfallFamiliesReady, taxon_cols_pitfall_families, "Pitfall assemblage (family level)")
hordbar=make_rank_bar_plot(HedgeOrdersReady, taxon_cols_hedge_orders, "Hedge beating assemblage (order level)")
hfambar=make_rank_bar_plot(HedgeFamiliesReady, taxon_cols_hedge_families, "Hedge beating assemblage (family level)")
carabar=make_rank_bar_plot(CarabidDataReady, carabid_cols, "Carabid assemblage (species level)")

pordbar
pfambar
hordbar
hfambar
carabar

pitfall_col <- "#4C78A8"   #blue
hedge_col   <- "#59A14F"   #green


############ 3.1 Pitfall Taxon Richness GLMMs ----
#GLMMs for pitfall taxon richness with compois distribution
m_richness_pitfall_orders <- glmmTMB(taxon_richness ~ d_c_sc + sampling_round + (1 | pitfall),
                       family = compois,
                       data = PitfallOrdersReady)

m_richness_pitfall_orders_int <- glmmTMB(taxon_richness ~ d_c_sc * sampling_round + (1 | pitfall),
                           family = compois,
                           data = PitfallOrdersReady)

m_richness_pitfall_fams <- glmmTMB(taxon_richness ~ d_c_sc + sampling_round + (1 | pitfall),
                     family = compois,
                     data = PitfallFamiliesReady)

m_richness_pitfall_fams_int <- glmmTMB(taxon_richness ~ d_c_sc * sampling_round + (1 | pitfall),
                         family = compois,
                         data = PitfallFamiliesReady)

#Additive better for both
AIC(m_richness_pitfall_orders)
AIC(m_richness_pitfall_orders_int)
AIC(m_richness_pitfall_fams)
AIC(m_richness_pitfall_fams_int)

#Both models pass diagnostic checks
summary(m_richness_pitfall_orders)
summary(m_richness_pitfall_fams)
check_glmm(m_richness_pitfall_orders)
check_glmm(m_richness_pitfall_fams)

extract_irr(m_richness_pitfall_fams)

#Sensitivity Analysis - adjacent land use
m_richness_pitfall_fams2 <- glmmTMB(taxon_richness ~ d_c_sc + sampling_round + land_use_bin + (1 | pitfall),
                      family = compois,
                      data = PitfallFamiliesReady)

#Incorporating adjacent land use doesn't decrease AIC, no sig dif with ANOVA
summary(m_richness_pitfall_fams2)
AIC(m_richness_pitfall_fams,m_richness_pitfall_fams2)
anova(m_richness_pitfall_fams,m_richness_pitfall_fams2)




#Family taxon richness
dc_mean <- mean(PitfallFamiliesReady$d_c, na.rm = TRUE)
dc_sd   <- sd(PitfallFamiliesReady$d_c, na.rm = TRUE)

pred <- ggpredict(m_richness_pitfall_fams, terms = c("d_c_sc [all]", "sampling_round"))
pred_df <- as.data.frame(pred)

#Convert scaled predictor back to d_c
pred_df <- pred_df %>%
  rename(d_c_sc = x, sampling_round = group) %>%
  mutate(d_c = d_c_sc * dc_sd + dc_mean)

#Determine breaks (multiples of 2)
max_dc <- max(PitfallFamiliesReady$d_c, na.rm = TRUE)
dc_breaks <- seq(0, ceiling(max_dc / 2) * 2, by = 2)

PitfallFamilyTaxonRich <- ggplot() +
  geom_jitter(data = PitfallFamiliesReady, aes(x = d_c, y = taxon_richness, color = sampling_round, shape = sampling_round),
              width = jitter_width, alpha = point_alpha, size = point_size) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_line(data = pred_df, aes(x = d_c, y = predicted, colour = sampling_round, linetype = sampling_round), linewidth = 1) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_ribbon(data = pred_df, aes(x = d_c, ymin = conf.low, ymax = conf.high, fill = sampling_round), alpha = 0.16, colour = NA) +
  scale_fill_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  scale_x_continuous(limits = c(3,25), breaks = dc_breaks, expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 14), expand = c(0,0)) +
  labs(x = "Distinctiveness*Condition Score",
       y = "Family-level Taxon Richness",
       #      title = "Pitfall Family-level taxon richness vs Distinctiveness*Condition",
       colour = "sampling_round", fill = "sampling_round", shape = "sampling_round", linetype = "sampling_round") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(1.4, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10), #change legend text font size
        axis.title = element_text(size = 15),  # axis labels (x and y titles)
        axis.text  = element_text(size = 12))
PitfallFamilyTaxonRich

ggsave(
  "Figure_PitfallFamilyTaxonRich_vs_BNG.png",
  plot = PitfallFamilyTaxonRich,
  width = 6.5,
  height = 5,
  dpi = 800)


#Order taxon richness
dc_mean <- mean(PitfallOrdersReady$d_c, na.rm = TRUE)
dc_sd   <- sd(PitfallOrdersReady$d_c, na.rm = TRUE)

pred <- ggpredict(m_richness_pitfall_orders, terms = c("d_c_sc [all]", "sampling_round"))
pred_df <- as.data.frame(pred)

#Convert the scaled predictor back to raw d_c
pred_df <- pred_df %>%
  rename(d_c_sc = x, sampling_round = group) %>%
  mutate(d_c = d_c_sc * dc_sd + dc_mean)

#Determine breaks 
max_dc <- max(PitfallOrdersReady$d_c, na.rm = TRUE)
dc_breaks <- seq(0, ceiling(max_dc / 2) * 2, by = 2)

PitfallOrdersTaxonRich <- ggplot() +
  geom_jitter(data = PitfallOrdersReady, aes(x = d_c, y = taxon_richness, color = sampling_round, shape = sampling_round),
              width = jitter_width, alpha = point_alpha, size = point_size) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_line(data = pred_df, aes(x = d_c, y = predicted, colour = sampling_round, linetype = sampling_round), linewidth = 1) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_ribbon(data = pred_df, aes(x = d_c, ymin = conf.low, ymax = conf.high, fill = sampling_round), alpha = 0.16, colour = NA) +
  scale_fill_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  scale_x_continuous(limits = c(3,25), breaks = dc_breaks, expand = c(0,0)) +
  scale_y_continuous(limits = c(3, 14), expand = c(0,0)) +
  labs(x = "Distinctiveness*Condition Score",
       y = "Order-level Taxon Richness",
       #      title = "Pitfall Order-level taxon richness vs Distinctiveness*Condition",
       colour = "sampling_round", fill = "sampling_round", shape = "sampling_round", linetype = "sampling_round") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(1.4, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10), #change legend text font size
        axis.title = element_text(size = 15),  # axis labels (x and y titles)
        axis.text  = element_text(size = 12))
PitfallOrdersTaxonRich

ggsave(
  "Figure_PitfallOrdersTaxonRich_vs_BNG.png",
  plot = PitfallOrdersTaxonRich,
  width = 6.5,
  height = 5,
  dpi = 800)


############ 3.2 Hedge Taxon Richness GLMMs ----
#GLMMs for hedge beating taxon richness with compois distribution
m_richness_hedge_orders <- glmmTMB(taxon_richness ~ d_c_sc + sampling_round + (1 | pitfall),
                             family = compois,
                             data = HedgeOrdersReady)

m_richness_hedge_orders_int <- glmmTMB(taxon_richness ~ d_c_sc * sampling_round + (1 | pitfall),
                                 family = compois,
                                 data = HedgeOrdersReady)

m_richness_hedge_fams <- glmmTMB(taxon_richness ~ d_c_sc + sampling_round + (1 | pitfall),
                           family = compois,
                           data = HedgeFamiliesReady)

m_richness_hedge_fams_int <- glmmTMB(taxon_richness ~ d_c_sc * sampling_round + (1 | pitfall),
                               family = compois,
                               data = HedgeFamiliesReady)

#Additive model AIC score higher for orders but lower for families
AIC(m_richness_hedge_orders)
AIC(m_richness_hedge_orders_int)
AIC(m_richness_hedge_fams)
AIC(m_richness_hedge_fams_int)
anova(m_richness_hedge_fams,m_richness_hedge_fams_int)


#Both models pass diagnostic checks
summary(m_richness_hedge_orders)
extract_irr(m_richness_hedge_fams_int)
summary(m_richness_hedge_fams_int)
summary(m_richness_hedge_fams)

#Interaction model despite having lower AIC violates several assumptions, so additive model retained
check_glmm(m_richness_hedge_fams)

#Sensitivity Analysis - adjacent land use
m_richness_hedge_fams2 <- glmmTMB(taxon_richness ~ d_c_sc + sampling_round + land_use_bin + (1 | pitfall),
                            family = compois,
                            data = HedgeFamiliesReady)

#Incorporating adjacent land use doesn't decrease AIC, no sig dif with ANOVA
summary(m_richness_hedge_fams2)
AIC(m_richness_hedge_fams2,m_richness_hedge_fams)
anova(m_richness_hedge_fams,m_richness_hedge_fams2)



#Family taxon richness
dc_mean <- mean(HedgeFamiliesReady$d_c, na.rm = TRUE)
dc_sd   <- sd(HedgeFamiliesReady$d_c, na.rm = TRUE)

pred <- ggpredict(m_richness_hedge_fams_int, terms = c("d_c_sc [all]", "sampling_round"))
pred_df <- as.data.frame(pred)

#convert the scaled predictor back to raw d_c
pred_df <- pred_df %>%
  rename(d_c_sc = x, sampling_round = group) %>%
  mutate(d_c = d_c_sc * dc_sd + dc_mean)

#determine breaks (multiples of 2) up to a sensible max
max_dc <- max(HedgeFamiliesReady$d_c, na.rm = TRUE)
dc_breaks <- seq(0, ceiling(max_dc / 2) * 2, by = 2)

HedgeFamilyTaxonRich <- ggplot() +
  geom_jitter(data = HedgeFamiliesReady, aes(x = d_c, y = taxon_richness, color = sampling_round, shape = sampling_round),
              width = jitter_width, alpha = point_alpha, size = point_size) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_line(data = pred_df, aes(x = d_c, y = predicted, colour = sampling_round, linetype = sampling_round), linewidth = 1) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_ribbon(data = pred_df, aes(x = d_c, ymin = conf.low, ymax = conf.high, fill = sampling_round), alpha = 0.16, colour = NA) +
  scale_fill_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  scale_x_continuous(limits = c(3,25), breaks = dc_breaks, expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 17), expand = c(0,0)) +
  labs(x = "Distinctiveness*Condition Score",
       y = "Family-level Taxon Richness",
       #      title = "Pitfall Family-level taxon richness vs Distinctiveness*Condition",
       colour = "sampling_round", fill = "sampling_round", shape = "sampling_round", linetype = "sampling_round") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(1.4, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10), #change legend text font size
        axis.title = element_text(size = 15),  # axis labels (x and y titles)
        axis.text  = element_text(size = 12))
HedgeFamilyTaxonRich

ggsave(
  "Figure_HedgeFamilyTaxonRich_vs_BNG.png",
  plot = HedgeFamilyTaxonRich,
  width = 6.5,
  height = 5,
  dpi = 800)




#Order taxon richness
dc_mean <- mean(HedgeOrdersReady$d_c, na.rm = TRUE)
dc_sd   <- sd(HedgeOrdersReady$d_c, na.rm = TRUE)

pred <- ggpredict(m_richness_hedge_orders, terms = c("d_c_sc [all]", "sampling_round"))
pred_df <- as.data.frame(pred)

#convert the scaled predictor back to raw d_c
pred_df <- pred_df %>%
  rename(d_c_sc = x, sampling_round = group) %>%
  mutate(d_c = d_c_sc * dc_sd + dc_mean)

#determine breaks (multiples of 2) up to a sensible max
max_dc <- max(HedgeOrdersReady$d_c, na.rm = TRUE)
dc_breaks <- seq(0, ceiling(max_dc / 2) * 2, by = 2)

HedgeOrdersTaxonRich <- ggplot() +
  geom_jitter(data = HedgeOrdersReady, aes(x = d_c, y = taxon_richness, color = sampling_round, shape = sampling_round),
              width = jitter_width, alpha = point_alpha, size = point_size) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_line(data = pred_df, aes(x = d_c, y = predicted, colour = sampling_round, linetype = sampling_round), linewidth = 1) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_ribbon(data = pred_df, aes(x = d_c, ymin = conf.low, ymax = conf.high, fill = sampling_round), alpha = 0.16, colour = NA) +
  scale_fill_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  scale_x_continuous(limits = c(3,25), breaks = dc_breaks, expand = c(0,0)) +
  scale_y_continuous(limits = c(3, 11), expand = c(0,0)) +
  labs(x = "Distinctiveness*Condition Score",
       y = "Order-level Taxon Richness",
       #      title = "Hedge Order-level taxon richness vs Distinctiveness*Condition",
       colour = "sampling_round", fill = "sampling_round", shape = "sampling_round", linetype = "sampling_round") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(1.4, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10), #change legend text font size
        axis.title = element_text(size = 15),  # axis labels (x and y titles)
        axis.text  = element_text(size = 12))
HedgeOrdersTaxonRich

ggsave(
  "Figure_HedgeOrdersTaxonRich_vs_BNG.png",
  plot = HedgeOrdersTaxonRich,
  width = 6.5,
  height = 5,
  dpi = 800)

############ 3.3 Carabid Species Richness GLMMs ----
#GLMMs for carabid taxon richness with compois distribution
m_richness_carab <- glmmTMB(taxon_richness ~ d_c_sc + sampling_round + (1 | pitfall),
                      family = compois,
                      data = CarabidDataReady)

m_richness_carab2 <- glmmTMB(taxon_richness ~ d_c_sc + sampling_round + land_use_bin + (1 | pitfall),
                       family = compois,
                       data = CarabidDataReady)

m_richness_carab_int <- glmmTMB(taxon_richness ~ d_c_sc * sampling_round + land_use_bin + (1 | pitfall),
                          family = compois,
                          data = CarabidDataReady)

#Models not significantly different, additive retained
AIC(m_richness_carab)
AIC(m_richness_carab2)
AIC(m_richness_carab_int)
anova(m_richness_carab, m_richness_carab2, m_richness_carab_int)

#Both models pass diagnostic checks
summary(m_richness_carab)
summary(m_richness_carab2)
summary(m_richness_carab_int)

extract_irr(m_richness_carab)

check_glmm(m_richness_carab)


#Recover scaling info for d_c (raw) from the same data used to fit the model
dc_mean <- mean(CarabidDataReady$d_c, na.rm = TRUE)
dc_sd   <- sd(CarabidDataReady$d_c, na.rm = TRUE)

#Determine breaks (multiples of 2) up to a sensible max
max_dc <- max(CarabidDataReady$d_c, na.rm = TRUE)
dc_breaks <- seq(0, ceiling(max_dc / 2) * 2, by = 2)

#Get model predictions across the observed range of d_c_sc and by sampling_round
pred <- ggpredict(m_richness_carab, terms = c("d_c_sc [all]", "sampling_round"))
pred_df <- as.data.frame(pred) %>%
  rename(d_c_sc = x, sampling_round = group) %>%
  mutate(d_c = d_c_sc * dc_sd + dc_mean)

#Original points (use the same dataframe used to fit the model)
plot_df <- CarabidDataReady

#Choose sensible x/y limits (optional — here we use data range +/- small margin)
x_min <- min(plot_df$d_c, na.rm = TRUE); x_max <- max(plot_df$d_c, na.rm = TRUE)
y_min <- min(plot_df$taxon_richness, na.rm = TRUE); y_max <- max(plot_df$taxon_richness, na.rm = TRUE)

#Plot
CarabidTaxonRich <- ggplot() +
  geom_jitter(data = plot_df, aes(x = d_c, y = taxon_richness, color = sampling_round, shape = sampling_round),
              width = 0.2, height = 0.2, alpha = point_alpha, size = point_size) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_line(data = pred_df, aes(x = d_c, y = predicted, colour = sampling_round, linetype = sampling_round), linewidth = 1) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_ribbon(data = pred_df, aes(x = d_c, ymin = conf.low, ymax = conf.high, fill = sampling_round), alpha = 0.2, colour = NA) +
  scale_fill_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  scale_shape_manual(values = c("One" = 16, "Two" = 17, "Three" = 15)) +
  scale_x_continuous(limits = c(3,25), breaks = dc_breaks, expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.5, 7.5), expand = c(0,0)) +
  labs(x = "Distinctiveness*Condition Score",
       y = "Carabid Species Richness",
       #title = "Carabid species richness vs BNG score",
       colour = "sampling_round", fill = "sampling_round", shape = "sampling_round", linetype = "sampling_round") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(1.4, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10), #change legend text font size
        axis.title = element_text(size = 15),  # axis labels (x and y titles)
        axis.text  = element_text(size = 12)) 
CarabidTaxonRich

#I had an issue with getting the colours to allocate to the same sampling rounds as my other plots, this seemed to work
lvl <- c("One", "Two", "Three")

CarabidDataReady$sampling_round <- factor(CarabidDataReady$sampling_round, levels = lvl)
pred_df$sampling_round <- factor(pred_df$sampling_round, levels = lvl)

CarabidTaxonRich <- ggplot() +
  geom_jitter(data = CarabidDataReady,
              aes(x = d_c, y = taxon_richness, color = sampling_round, shape = sampling_round),
              width = jitter_width, alpha = point_alpha, size = point_size) +
  geom_line(data = pred_df,
            aes(x = d_c, y = predicted, colour = sampling_round, linetype = sampling_round),
            linewidth = 1) +
  geom_ribbon(data = pred_df,
              aes(x = d_c, ymin = conf.low, ymax = conf.high, fill = sampling_round),
              alpha = 0.2, colour = NA) +
  scale_color_manual(values = c("One" = "#ffb000", "Two" = "#dc267f", "Three" = "#648fff"),
                     breaks = lvl, limits = lvl) +
  scale_fill_manual(values = c("One" = "#ffb000", "Two" = "#dc267f", "Three" = "#648fff"),
                    breaks = lvl, limits = lvl) +
  scale_shape_manual(values = c("One" = 16, "Two" = 17, "Three" = 15),
                     breaks = lvl, limits = lvl) +
  #scale_linetype_manual(values = c("One" = "solid", "Two" = "dashed", "Three" = "dashed"),
  #                      breaks = lvl, limits = lvl) +
  scale_x_continuous(limits = c(3, 25), breaks = dc_breaks, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 7.2), expand = c(0, 0)) +
  labs(x = "Distinctiveness*Condition Score",
       y = "Carabid Species Richness",
       colour = "sampling_round", shape = "sampling_round",
       fill = "sampling_round", linetype = "sampling_round") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(1.4, "cm"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))
CarabidTaxonRich


ggsave(
  "Figure_CarabidTaxonRich_vs_BNG.png",
  plot = CarabidTaxonRich,
  width = 6.5,
  height = 5,
  dpi = 800)


#Shannon models
CarabidDataShannon <- CarabidDataReady %>%
  filter(!is.na(carab_shannon))

m_carab_shannon <- lmerTest::lmer(
  carab_shannon ~ d_c_sc + sampling_round + (1 | pitfall),
  data = CarabidDataShannon)

m_carab_shannon2 <- lmerTest::lmer(
  carab_shannon ~ d_c_sc + sampling_round + land_use_bin + (1 | pitfall),
  data = CarabidDataShannon)

m_carab_shannon_int <- lmerTest::lmer(
  carab_shannon ~ d_c_sc * sampling_round + (1 | pitfall),
  data = CarabidDataShannon)

#Additive model has lower AIC (127 vs 134)
AIC(m_carab_shannon)
AIC(m_carab_shannon2)
AIC(m_carab_shannon_int)

#Both models pass everything apart from zero inflation, but probably still useable 
summary(m_carab_shannon)
summary(m_carab_shannon_int)

check_glmm(m_carab_shannon)


#Carabid Shannon diversity vs D×C (coloured by sampling_round)
# recover scaling used for d_c
dc_mean <- mean(CarabidDataShannon$d_c, na.rm = TRUE)
dc_sd   <- sd(CarabidDataShannon$d_c, na.rm = TRUE)

#Determine breaks (multiples of 2) up to a sensible max
max_dc <- max(CarabidDataShannon$d_c, na.rm = TRUE)
dc_breaks <- seq(0, ceiling(max_dc / 2) * 2, by = 2)

#Get predictions across d_c_sc and by sampling_round
pred <- ggpredict(m_carab_shannon, terms = c("d_c_sc [all]", "sampling_round"))
pred_df <- as.data.frame(pred) %>%
  rename(d_c_sc = x, sampling_round = group) %>%
  mutate(d_c = d_c_sc * dc_sd + dc_mean)

#Original points (only those with shannon present)
plot_df <- CarabidDataShannon %>% filter(!is.na(carab_shannon))

#Sensible axis ranges (small margins)
x_min <- min(plot_df$d_c, na.rm = TRUE); x_max <- max(plot_df$d_c, na.rm = TRUE)
y_min <- min(plot_df$carab_shannon, na.rm = TRUE); y_max <- max(plot_df$carab_shannon, na.rm = TRUE)

#Plot
CarabidShannon <- ggplot() +
  geom_jitter(data = plot_df, aes(x = d_c, y = carab_shannon, color = sampling_round, shape = sampling_round),
              width = jitter_width, alpha = point_alpha, size = point_size) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_line(data = pred_df, aes(x = d_c, y = predicted, colour = sampling_round, linetype = sampling_round), linewidth = 1) +
  scale_color_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  geom_ribbon(data = pred_df, aes(x = d_c, ymin = conf.low, ymax = conf.high, fill = sampling_round), alpha = 0.16, colour = NA) +
  scale_fill_manual(values = c("One"="#ffb000","Two"="#dc267f","Three"="#648fff")) +
  scale_x_continuous(limits = c(3,25), breaks = dc_breaks, expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.1, 1.6), expand = c(0,0)) +
  labs(x = "Distinctiveness*Condition Score",
       y = "Carabid Shannon Diversity",
       #title = "Carabid Shannon diversity vs BNG score",
       colour = "sampling_round", fill = "sampling_round", shape = "sampling_round", linetype = "sampling_round") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(1.4, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10), #change legend text font size
        axis.title = element_text(size = 15),  # axis labels (x and y titles)
        axis.text  = element_text(size = 12))   # tick labels (numbers)
CarabidShannon

lvl <- c("One", "Two", "Three")

CarabidDataReady$sampling_round <- factor(CarabidDataReady$sampling_round, levels = lvl)
pred_df$sampling_round <- factor(pred_df$sampling_round, levels = lvl)

CarabidShannon <- ggplot() +
  geom_jitter(data = CarabidDataReady,
              aes(x = d_c, y = carab_shannon, color = sampling_round, shape = sampling_round),
              width = jitter_width, alpha = point_alpha, size = point_size) +
  geom_line(data = pred_df,
            aes(x = d_c, y = predicted, colour = sampling_round, linetype = sampling_round),
            linewidth = 1) +
  geom_ribbon(data = pred_df,
              aes(x = d_c, ymin = conf.low, ymax = conf.high, fill = sampling_round),
              alpha = 0.2, colour = NA) +
  scale_color_manual(values = c("One" = "#ffb000", "Two" = "#dc267f", "Three" = "#648fff"),
                     breaks = lvl, limits = lvl) +
  scale_fill_manual(values = c("One" = "#ffb000", "Two" = "#dc267f", "Three" = "#648fff"),
                    breaks = lvl, limits = lvl) +
  scale_shape_manual(values = c("One" = 16, "Two" = 17, "Three" = 15),
                     breaks = lvl, limits = lvl) +
  #scale_linetype_manual(values = c("One" = "solid", "Two" = "dashed", "Three" = "dashed"),
  #                      breaks = lvl, limits = lvl) +
  scale_x_continuous(limits = c(3, 25), breaks = dc_breaks, expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.1, 1.6), expand = c(0, 0)) +
  labs(x = "Distinctiveness*Condition Score",
       y = "Carabid Shannon Diversity",
       colour = "sampling_round", shape = "sampling_round",
       fill = "sampling_round", linetype = "sampling_round") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(1.4, "cm"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))
CarabidShannon


ggsave(
  "Figure_CarabidShannon_vs_BNG.png",
  plot = CarabidShannon,
  width = 6.5,
  height = 5,
  dpi = 800)

############ 4.1 Pitfall Orders PERMANOVA and NMDS plots ---- 
#Mantel Test, compare community dissimilarity matrix with spatial matrix 
PitfallMean <- PitfallOrdersReady %>%
  group_by(pitfall) %>%
  summarise(
    across(all_of(taxon_cols_pitfall_orders), mean),
    latitude = first(latitude),
    longitude = first(longitude),
    .groups = "drop")

#Community matrix 
comm_mean <- as.matrix(PitfallMean[, taxon_cols_pitfall_orders])
comm_hel_mean <- decostand(comm_mean, method = "hellinger")
bray_mean <- vegdist(comm_hel_mean, method = "bray")

#Spatial matrix
spat_dist <- geosphere::distm(
  PitfallMean[, c("longitude","latitude")],
  fun = geosphere::distHaversine)
spat_dist <- as.dist(spat_dist)

set.seed(123)
mantel(bray_mean, spat_dist, method = "pearson", permutations = 999)


#PERMANOVA 
comm <- as.matrix(PitfallOrdersReady[, taxon_cols_pitfall_orders])
comm_hel <- decostand(comm, method = "hellinger")
bray_dist <- vegdist(comm_hel, method = "bray")
PitfallOrdersReady$pitfall <- as.factor(PitfallOrdersReady$pitfall)
PitfallOrdersReady$sampling_round <- as.factor(PitfallOrdersReady$sampling_round)
PitfallOrdersReady$condition <- as.factor(PitfallOrdersReady$condition)

set.seed(123)
perm <- how(nperm = 9999)
setBlocks(perm) <- PitfallOrdersReady$pitfall #Set blocks by pitfall due to non-independence

#Passes diagnostic test
disp <- betadisper(bray_dist, PitfallOrdersReady$pitfall)
anova(disp)

adonis_result_add <- adonis2(bray_dist ~ habitat_type + condition + sampling_round + land_use_bin,
                             data = PitfallOrdersReady,
                             permutations = perm,
                             by = "margin")
adonis_result_add

adonis_result_int <- adonis2(bray_dist ~ habitat_type * condition + sampling_round + land_use_bin,
                             data = PitfallOrdersReady,
                             permutations = perm,
                             by = "margin")
adonis_result_int

perm_free <- how(nperm = 9999)  # free permutations
adonis_pitfall <- adonis2(bray_dist ~ pitfall, data = PitfallOrdersReady, permutations = perm_free, by = "margin")
cat("\nPERMANOVA testing pitfall identity (no blocking):\n"); print(adonis_pitfall)

#Simper Test
simper(comm_hel,
       permutations = 999)


#NMDS (3D fit for better stress; will plot axes 1 & 2)
#NMDS plotting workflow refined with AI-assisted troubleshooting
set.seed(111)
nmds <- metaMDS(comm_hel, distance = "bray", k = 3, trymax = 100, trace = FALSE)

#scores -> plot_df (scores + metadata)
site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$sample_id <- rownames(site_scores)
meta <- PitfallOrdersReady %>% mutate(sample_id = rownames(.))
plot_df <- left_join(site_scores, meta, by = "sample_id")

#map habitat->distinctiveness and create ordered factor with labels
m <- c(
  "Native hedgerow"=2,
  "Species-rich native hedgerow"=4,
  "Species-rich native hedgerow with trees"=6,
  "Species-rich native hedgerow with trees - associated with bank or ditch"=8)

plot_df <- plot_df %>%
  mutate(distinctiveness = factor(m[as.character(habitat_type)],
                                  levels=c(2,4,6,8),
                                  labels=c(
                                    "Native hedgerow (2)",
                                    "Species-rich native hedgerow (4)",
                                    "Species-rich native hedgerow with trees (6)",
                                    "SR hedgerow + trees + bank/ditch (8)")))

#colour vector and convex hulls
cols <- c("#ffb000","#fe6100","#dc267f","#785ef0")
names(cols) <- levels(plot_df$distinctiveness)
hulls <- plot_df %>%
  group_by(distinctiveness) %>%
  group_modify(~ if(nrow(.x)>=3) .x[chull(.x$NMDS1,.x$NMDS2),] else tibble()) %>%
  ungroup()

#Dissertation plot, optional line type change around convex hulls
PitfallOrdersNMDS <- ggplot(plot_df, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = distinctiveness, shape = sampling_round), size = 2.5, alpha = 0.9) +
  stat_ellipse(aes(colour = distinctiveness),
               type = "t",
               linewidth = 1.2,
               level = 0.95) + 
  scale_colour_manual(values = cols, guide = "none") +
  scale_linetype_manual(values = c("solid", "dashed", "twodash", "dotdash"), guide = "none") +
  scale_shape_discrete(guide = "none") + 
  theme_minimal() +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9))
PitfallOrdersNMDS

ggsave(
  "Figure_PitfallOrdersNMDS_vs_BNG.png",
  plot = PitfallOrdersNMDS,
  width = 5,
  height = 5,
  dpi = 600)

############ 4.2 Pitfall Families PERMANOVA and NMDS plots ---- 
#Mantel Test
PitfallMeanFamily <- PitfallFamiliesReady %>%
  group_by(pitfall) %>%
  summarise(
    across(all_of(taxon_cols_pitfall_families), mean),
    latitude = first(latitude),
    longitude = first(longitude),
    .groups = "drop")

comm_mean <- as.matrix(PitfallMeanFamily[, taxon_cols_pitfall_families])
comm_hel_mean <- decostand(comm_mean, method = "hellinger")
bray_mean <- vegdist(comm_hel_mean, method = "bray")

spat_dist <- geosphere::distm(
  PitfallMeanFamily[, c("longitude","latitude")],
  fun = geosphere::distHaversine)

spat_dist <- as.dist(spat_dist)

set.seed(123)
mantel(bray_mean, spat_dist, method = "pearson", permutations = 999)

#PERMANOVA 
comm <- as.matrix(PitfallFamiliesReady[, taxon_cols_pitfall_families])
comm_hel <- decostand(comm, method = "hellinger")
bray_dist <- vegdist(comm_hel, method = "bray")
PitfallFamiliesReady$pitfall <- as.factor(PitfallFamiliesReady$pitfall)
PitfallFamiliesReady$sampling_round <- as.factor(PitfallFamiliesReady$sampling_round)
PitfallFamiliesReady$condition <- as.factor(PitfallFamiliesReady$condition)

set.seed(123)
perm <- how(nperm = 9999)
setBlocks(perm) <- PitfallFamiliesReady$pitfall

#Passes diagnostic test
disp <- betadisper(bray_dist, PitfallFamiliesReady$pitfall)
anova(disp)

adonis_result_add <- adonis2(bray_dist ~ habitat_type + condition + sampling_round + land_use_bin,
                             data = PitfallFamiliesReady,
                             permutations = perm,
                             by = "margin")
adonis_result_add

adonis_result_int <- adonis2(bray_dist ~ habitat_type * condition + sampling_round + land_use_bin,
                             data = PitfallFamiliesReady,
                             permutations = perm,
                             by = "margin")
adonis_result_int

perm_free <- how(nperm = 9999)  # free permutations
adonis_pitfall <- adonis2(bray_dist ~ pitfall, data = PitfallFamiliesReady, permutations = perm_free, by = "margin")
cat("\nPERMANOVA testing pitfall identity (no blocking):\n"); print(adonis_pitfall)

#Simper test
simper(comm_hel,
       permutations = 999)

#NMDS (3D fit for better stress; will plot axes 1 & 2)
set.seed(111)
nmds <- metaMDS(comm_hel, distance = "bray", k = 3, trymax = 100, trace = FALSE)

#scores -> plot_df (scores + metadata)
site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$sample_id <- rownames(site_scores)
meta <- PitfallFamiliesReady %>% mutate(sample_id = rownames(.))
plot_df <- left_join(site_scores, meta, by = "sample_id")

#map habitat->distinctiveness and create ordered factor with labels
m <- c(
  "Native hedgerow"=2,
  "Species-rich native hedgerow"=4,
  "Species-rich native hedgerow with trees"=6,
  "Species-rich native hedgerow with trees - associated with bank or ditch"=8)

plot_df <- plot_df %>%
  mutate(distinctiveness = factor(m[as.character(habitat_type)],
                                  levels=c(2,4,6,8),
                                  labels=c(
                                    "Native hedgerow (2)",
                                    "Species-rich native hedgerow (4)",
                                    "Species-rich native hedgerow with trees (6)",
                                    "SR hedgerow + trees + bank/ditch (8)")))

#colour vector and convex hulls
cols <- c("#ffb000","#fe6100","#dc267f","#785ef0")
names(cols) <- levels(plot_df$distinctiveness)
hulls <- plot_df %>%
  group_by(distinctiveness) %>%
  group_modify(~ if(nrow(.x)>=3) .x[chull(.x$NMDS1,.x$NMDS2),] else tibble()) %>%
  ungroup()


#Dissertation plot
PitfallFamiliesNMDS <- ggplot(plot_df, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = distinctiveness, shape = sampling_round), size = 2.5, alpha = 0.9) +
  stat_ellipse(aes(colour = distinctiveness),
               type = "t",
               linewidth = 1.2,
               level = 0.95) + 
  scale_colour_manual(values = cols, guide = "none") +
  scale_linetype_manual(values = c("solid", "dashed", "twodash", "dotdash"), guide = "none") +
  scale_shape_discrete(guide = "none") + 
  theme_minimal() +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9))
PitfallFamiliesNMDS

ggsave(
  "Figure_PitfallFamiliesNMDS_vs_BNG.png",
  plot = PitfallFamiliesNMDS,
  width = 5,
  height = 5,
  dpi = 600)

############ 4.3 Hedge Orders PERMANOVA and NMDS plots ---- 
#Mantel Test
HedgeMeanOrders <- HedgeOrdersReady %>%
  group_by(pitfall) %>%
  summarise(
    across(all_of(taxon_cols_hedge_orders), mean),
    latitude = first(latitude),
    longitude = first(longitude),
    .groups = "drop")

comm_mean <- as.matrix(HedgeMeanOrders[, taxon_cols_hedge_orders])
comm_hel_mean <- decostand(comm_mean, method = "hellinger")
bray_mean <- vegdist(comm_hel_mean, method = "bray")

spat_dist <- geosphere::distm(
  HedgeMeanOrders[, c("longitude","latitude")],
  fun = geosphere::distHaversine)

spat_dist <- as.dist(spat_dist)

set.seed(123)
mantel(bray_mean, spat_dist, method = "pearson", permutations = 999)

#PERMANOVA 
comm <- as.matrix(HedgeOrdersReady[, taxon_cols_hedge_orders])
comm_hel <- decostand(comm, method = "hellinger")
bray_dist <- vegdist(comm_hel, method = "bray")
HedgeOrdersReady$pitfall <- as.factor(HedgeOrdersReady$pitfall)
HedgeOrdersReady$sampling_round <- as.factor(HedgeOrdersReady$sampling_round)
HedgeOrdersReady$condition <- as.factor(HedgeOrdersReady$condition)

set.seed(123)
perm <- how(nperm = 9999)
setBlocks(perm) <- HedgeOrdersReady$pitfall

#Passes diagnostic test
disp <- betadisper(bray_dist, HedgeOrdersReady$pitfall)
anova(disp)

adonis_result_add <- adonis2(bray_dist ~ habitat_type + condition + sampling_round + land_use_bin,
                             data = HedgeOrdersReady,
                             permutations = perm,
                             by = "margin")
adonis_result_add

adonis_result_int <- adonis2(bray_dist ~ habitat_type * condition + sampling_round + land_use_bin,
                             data = HedgeOrdersReady,
                             permutations = perm,
                             by = "margin")
adonis_result_int

perm_free <- how(nperm = 9999)  # free permutations
adonis_pitfall <- adonis2(bray_dist ~ pitfall, data = HedgeOrdersReady, permutations = perm_free, by = "margin")
cat("\nPERMANOVA testing pitfall identity (no blocking):\n"); print(adonis_pitfall)

#Simper test
simper(comm_hel,
       permutations = 999)

#NMDS (3D fit for better stress; will plot axes 1 & 2)
set.seed(111)
nmds <- metaMDS(comm_hel, distance = "bray", k = 3, trymax = 100, trace = FALSE)

#scores -> plot_df (scores + metadata)
site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$sample_id <- rownames(site_scores)
meta <- HedgeOrdersReady %>% mutate(sample_id = rownames(.))
plot_df <- left_join(site_scores, meta, by = "sample_id")

#Map habitat->distinctiveness and create ordered factor with labels
m <- c(
  "Native hedgerow"=2,
  "Species-rich native hedgerow"=4,
  "Species-rich native hedgerow with trees"=6,
  "Species-rich native hedgerow with trees - associated with bank or ditch"=8)

plot_df <- plot_df %>%
  mutate(distinctiveness = factor(m[as.character(habitat_type)],
                                  levels=c(2,4,6,8),
                                  labels=c(
                                    "Native hedgerow (2)",
                                    "Species-rich native hedgerow (4)",
                                    "Species-rich native hedgerow with trees (6)",
                                    "SR hedgerow + trees + bank/ditch (8)")))

#colour vector and convex hulls
cols <- c("#ffb000","#fe6100","#dc267f","#785ef0")
names(cols) <- levels(plot_df$distinctiveness)
hulls <- plot_df %>%
  group_by(distinctiveness) %>%
  group_modify(~ if(nrow(.x)>=3) .x[chull(.x$NMDS1,.x$NMDS2),] else tibble()) %>%
  ungroup()


#Dissertation plot
HedgeOrdersNMDS <- ggplot(plot_df, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = distinctiveness, shape = sampling_round), size = 2.5, alpha = 0.9) +
  stat_ellipse(aes(colour = distinctiveness),
               type = "t",
               linewidth = 1.2,
               level = 0.95) + 
  scale_colour_manual(values = cols, guide = "none") +
  scale_linetype_manual(values = c("solid", "dashed", "twodash", "dotdash"), guide = "none") +
  scale_shape_discrete(guide = "none") + 
  theme_minimal() +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9))
HedgeOrdersNMDS

ggsave(
  "Figure_HedgeOrdersNMDS_vs_BNG.png",
  plot = HedgeOrdersNMDS,
  width = 5,
  height = 5,
  dpi = 600)

############ 4.4 Hedge Families PERMANOVA and NMDS plots ---- 
#Mantel Test
HedgeMeanFamilies <- HedgeFamiliesReady %>%
  group_by(pitfall) %>%
  summarise(
    across(all_of(taxon_cols_hedge_families), mean),
    latitude = first(latitude),
    longitude = first(longitude),
    .groups = "drop")

comm_mean <- as.matrix(HedgeMeanFamilies[, taxon_cols_hedge_families])
comm_hel_mean <- decostand(comm_mean, method = "hellinger")
bray_mean <- vegdist(comm_hel_mean, method = "bray")

spat_dist <- geosphere::distm(
  HedgeMeanFamilies[, c("longitude","latitude")],
  fun = geosphere::distHaversine)

spat_dist <- as.dist(spat_dist)

set.seed(123)
mantel(bray_mean, spat_dist, method = "pearson", permutations = 999)

#PERMANOVA 
comm <- as.matrix(HedgeFamiliesReady[, taxon_cols_hedge_families])
comm_hel <- decostand(comm, method = "hellinger")
bray_dist <- vegdist(comm_hel, method = "bray")
HedgeFamiliesReady$pitfall <- as.factor(HedgeFamiliesReady$pitfall)
HedgeFamiliesReady$sampling_round <- as.factor(HedgeFamiliesReady$sampling_round)
HedgeFamiliesReady$condition <- as.factor(HedgeFamiliesReady$condition)

set.seed(123)
perm <- how(nperm = 9999)
setBlocks(perm) <- HedgeFamiliesReady$pitfall

#Passes diagnostic test
disp <- betadisper(bray_dist, HedgeFamiliesReady$pitfall)
anova(disp)

adonis_result_add <- adonis2(bray_dist ~ habitat_type + condition + sampling_round + land_use_bin,
                             data = HedgeFamiliesReady,
                             permutations = perm,
                             by = "margin")
adonis_result_add

adonis_result_int <- adonis2(bray_dist ~ habitat_type * condition + sampling_round + land_use_bin,
                             data = HedgeFamiliesReady,
                             permutations = perm,
                             by = "margin")
adonis_result_int

perm_free <- how(nperm = 9999)  # free permutations
adonis_pitfall <- adonis2(bray_dist ~ pitfall, data = HedgeFamiliesReady, permutations = perm_free, by = "margin")
cat("\nPERMANOVA testing pitfall identity (no blocking):\n"); print(adonis_pitfall)

#Simper
simper(comm_hel,
       permutations = 999)

#NMDS (3D fit for better stress; will plot axes 1 & 2)
set.seed(111)
nmds <- metaMDS(comm_hel, distance = "bray", k = 3, trymax = 100, trace = FALSE)

#scores -> plot_df (scores + metadata)
site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$sample_id <- rownames(site_scores)
meta <- HedgeFamiliesReady %>% mutate(sample_id = rownames(.))
plot_df <- left_join(site_scores, meta, by = "sample_id")

#map habitat->distinctiveness and create ordered factor with labels
m <- c(
  "Native hedgerow"=2,
  "Species-rich native hedgerow"=4,
  "Species-rich native hedgerow with trees"=6,
  "Species-rich native hedgerow with trees - associated with bank or ditch"=8)

plot_df <- plot_df %>%
  mutate(distinctiveness = factor(m[as.character(habitat_type)],
                                  levels=c(2,4,6,8),
                                  labels=c(
                                    "Native hedgerow (2)",
                                    "Species-rich native hedgerow (4)",
                                    "Species-rich native hedgerow with trees (6)",
                                    "Species-rich native hedgerow with trees + bank/ditch (8)")))

#colour vector and convex hulls
cols <- c("#ffb000","#fe6100","#dc267f","#785ef0")
names(cols) <- levels(plot_df$distinctiveness)
hulls <- plot_df %>%
  group_by(distinctiveness) %>%
  group_modify(~ if(nrow(.x)>=3) .x[chull(.x$NMDS1,.x$NMDS2),] else tibble()) %>%
  ungroup()

#Dissertation plot
HedgeFamiliesNMDS <- ggplot(plot_df, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = distinctiveness, shape = sampling_round), size = 2.5, alpha = 0.9) +
  stat_ellipse(aes(colour = distinctiveness),
               type = "t",
               linewidth = 1.2,
               level = 0.95) + 
  scale_colour_manual(values = cols, guide = "none") +
  scale_linetype_manual(values = c("solid", "dashed", "twodash", "dotdash"), guide = "none") +
  scale_shape_discrete(guide = "none") + 
  theme_minimal() +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9))
HedgeFamiliesNMDS

ggsave(
  "Figure_HedgeFamiliesNMDS_vs_BNG.png",
  plot = HedgeFamiliesNMDS,
  width = 5,
  height = 5,
  dpi = 600)

############ 4.5 Carabids PERMANOVA and NMDS plots ----
#Mantel Test
MeanCarabid <- CarabidDataReady %>%
  group_by(pitfall) %>%
  summarise(
    across(all_of(carabid_cols), mean),
    latitude = first(latitude),
    longitude = first(longitude),
    .groups = "drop")

comm_mean <- as.matrix(MeanCarabid[, carabid_cols])
comm_hel_mean <- decostand(comm_mean, method = "hellinger")
bray_mean <- vegdist(comm_hel_mean, method = "bray")

spat_dist <- geosphere::distm(
  MeanCarabid[, c("longitude","latitude")],
  fun = geosphere::distHaversine)

spat_dist <- as.dist(spat_dist)

set.seed(123)
mantel(bray_mean, spat_dist, method = "pearson", permutations = 999)

#PERMANOVA 
comm <- as.matrix(CarabidDataReady[, carabid_cols])

#Remove empty samples
keep <- rowSums(comm) > 0
comm <- comm[keep, ]
CarabidDataReady <- CarabidDataReady[keep, ]

#Replace NA if present
comm[is.na(comm)] <- 0

#Transform + distance
comm_hel <- decostand(comm, method = "hellinger")
bray_dist <- vegdist(comm_hel, method = "bray")
CarabidDataReady$pitfall <- as.factor(CarabidDataReady$pitfall)
CarabidDataReady$sampling_round <- as.factor(CarabidDataReady$sampling_round)
CarabidDataReady$condition <- as.factor(CarabidDataReady$condition)

set.seed(123)
perm <- how(nperm = 9999)
setBlocks(perm) <- CarabidDataReady$pitfall

#Passes diagnostic test
disp <- betadisper(bray_dist, CarabidDataReady$pitfall)
anova(disp)

adonis_result_add <- adonis2(bray_dist ~ habitat_type + condition + sampling_round + land_use_bin,
                             data = CarabidDataReady,
                             permutations = perm,
                             by = "margin")
adonis_result_add

adonis_result_int <- adonis2(bray_dist ~ habitat_type * condition + sampling_round + land_use_bin,
                             data = CarabidDataReady,
                             permutations = perm,
                             by = "margin")
adonis_result_int

perm_free <- how(nperm = 9999)  # free permutations
adonis_pitfall <- adonis2(bray_dist ~ pitfall, data = CarabidDataReady, permutations = perm_free, by = "margin")
cat("\nPERMANOVA testing pitfall identity (no blocking):\n"); print(adonis_pitfall)

#Simper test
simper(comm_hel,
       permutations = 999)

#NMDS (3D fit for better stress; will plot axes 1 & 2)
set.seed(111)
nmds <- metaMDS(comm_hel, distance = "bray", k = 3, trymax = 100, trace = FALSE)
nmds
#scores -> plot_df (scores + metadata)
site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$sample_id <- rownames(site_scores)
meta <- CarabidDataReady %>% mutate(sample_id = rownames(.))
plot_df <- left_join(site_scores, meta, by = "sample_id")

#map habitat->distinctiveness and create ordered factor with labels
m <- c(
  "Native hedgerow"=2,
  "Species-rich native hedgerow"=4,
  "Species-rich native hedgerow with trees"=6,
  "Species-rich native hedgerow with trees - associated with bank or ditch"=8)

plot_df <- plot_df %>%
  mutate(distinctiveness = factor(m[as.character(habitat_type)],
                                  levels=c(2,4,6,8),
                                  labels=c(
                                    "Native hedgerow (2)",
                                    "Species-rich native hedgerow (4)",
                                    "Species-rich native hedgerow with trees (6)",
                                    "SR hedgerow + trees + bank/ditch (8)")))

#colour vector and convex hulls
cols <- c("#ffb000","#fe6100","#dc267f","#785ef0")
names(cols) <- levels(plot_df$distinctiveness)
hulls <- plot_df %>%
  group_by(distinctiveness) %>%
  group_modify(~ if(nrow(.x)>=3) .x[chull(.x$NMDS1,.x$NMDS2),] else tibble()) %>%
  ungroup()


#Dissertation plot
CarabidNMDS <- ggplot(plot_df, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = distinctiveness, shape = sampling_round), size = 2.5, alpha = 0.9) +
  stat_ellipse(aes(colour = distinctiveness),
               type = "t",
               linewidth = 1.2,
               level = 0.95) + 
  scale_colour_manual(values = cols, guide = "none") +
  scale_linetype_manual(values = c("solid", "dashed", "twodash", "dotdash"), guide = "none") +
  scale_shape_discrete(guide = "none") + 
  theme_minimal() +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9))
CarabidNMDS

ggsave(
  "Figure_CarabidNMDS_vs_BNG.png",
  plot = CarabidNMDS,
  width = 5,
  height = 5,
  dpi = 600)

############ 5.1 Predictive Hedgerow Characteristics Models ----
#Abundance data and models
file <- "MasterDataNew.xlsx"
sheet <- "HedgeCharacteristicsModels"

df_hedge <- read_excel(file, sheet = sheet) %>%
  clean_names()

#Prepare factors and derived land-use variable
df_hedge <- df_hedge %>%
  mutate(
    pitfall = factor(pitfall),
    sampling_round = factor(sampling_round),
    adjacent_land_use = factor(adjacent_land_use),
    direction = factor(direction),
    land_use_bin = case_when(
      adjacent_land_use %in% c("Wheat Field", "Oat Field", "Barley field") ~ "agricultural_field",
      TRUE ~ "non_agricultural_grassland"
    ) %>% factor(levels = c("agricultural_field", "non_agricultural_grassland")))

#z-scale continuous predictors
zscale <- function(x) as.numeric(scale(x))

cont_vars <- c(
  "width", "height", "ground_flora_rich", "average_height", "max_vegetation_height",
  "vegetation_length_from_hedgerow", "aspect", "total_woody",
  "total_shrubs", "total_climbers", "total_trees")

cont_vars <- intersect(cont_vars, names(df_hedge))

df_hedge <- df_hedge %>%
  mutate(across(all_of(cont_vars), zscale))

#This is the sampling unit
site_var <- "pitfall"

hedge_abund_set <- dredge(
  glmmTMB(total_abundance_h ~ height + width + land_use_bin +
            total_woody +
            vegetation_length_from_hedgerow + sampling_round +
            (1|pitfall),
          data = df_hedge,
          family = nbinom2()),
  rank = "AICc")

#taxon_richness_h_fam for family level, taxon_richness_h_ord for order level
hedge_rich_set <- dredge(
  glmmTMB(taxon_richness_h_ord ~ height + width + land_use_bin +
            total_woody +
            vegetation_length_from_hedgerow + sampling_round +
            (1|pitfall),
          data = df_hedge,
          family = compois()),
  rank = "AICc")


#Taxon richness models
#Read and clean
file_rich <- "MasterDataNew.xlsx"
sheet_rich <- "PitfallCharacteristicsModels"

df_pitfall <- read_excel(file_rich, sheet = sheet_rich) %>%
  clean_names()

#Prepare factors and derived land-use variable
df_pitfall <- df_pitfall %>%
  mutate(
    pitfall = factor(pitfall),
    sampling_round = factor(sampling_round),
    adjacent_land_use = factor(adjacent_land_use),
    direction = factor(direction),
    land_use_bin = case_when(
      adjacent_land_use %in% c("Wheat Field", "Oat Field", "Barley field") ~ "agricultural_field",
      TRUE ~ "non_agricultural_grassland"
    ) %>% factor(levels = c("agricultural_field", "non_agricultural_grassland")))

#z-scale continuous predictors
zscale <- function(x) as.numeric(scale(x))

cont_vars <- c(
  "width", "height", "ground_flora_rich", "average_height", "max_vegetation_height",
  "vegetation_length_from_hedgerow", "aspect", "total_woody",
  "total_shrubs", "total_climbers", "total_trees")

cont_vars <- intersect(cont_vars, names(df_pitfall))

df_pitfall <- df_pitfall %>%
  mutate(across(all_of(cont_vars), zscale))

pitfall_abund_set <- dredge(
  glmmTMB(total_abundance_p ~ height + width + land_use_bin +
            ground_flora_rich + average_height +
            vegetation_length_from_hedgerow + sampling_round +
            (1|pitfall),
          data = df_pitfall,
          family = nbinom2()),
  rank = "AICc")

#taxon_richness_p_fam for family level taxon_richness_p_ord for order level
pitfall_rich_set <- dredge(
  glmmTMB(taxon_richness_p_ord ~ height + width + land_use_bin +
            ground_flora_rich + average_height +
            vegetation_length_from_hedgerow + sampling_round +
            (1|pitfall),
          data = df_pitfall,
          family = compois()),
  rank = "AICc")


#Carabid richness models
file_cara <- "MasterDataNew.xlsx"
sheet_cara <- "CarabidsConducive"

df_cara <- read_excel(file_cara, sheet = sheet_cara) %>%
  clean_names()

#Prepare factors and derived land-use variable
df_cara <- df_cara %>%
  mutate(
    pitfall = factor(pitfall),
    sampling_round = factor(sampling_round),
    adjacent_land_use = factor(adjacent_land_use),
    direction = factor(direction),
    land_use_bin = case_when(
      adjacent_land_use %in% c("Wheat Field", "Oat Field", "Barley field") ~ "agricultural_field",
      TRUE ~ "non_agricultural_grassland"
    ) %>% factor(levels = c("agricultural_field", "non_agricultural_grassland")))

#z-scale only truly continuous predictors
zscale <- function(x) as.numeric(scale(x))

cont_vars <- c(
  "width", "height", "ground_flora_rich", "average_height", "max_vegetation_height",
  "vegetation_length_from_hedgerow", "aspect", "total_woody",
  "total_shrubs", "total_climbers", "total_trees")

cont_vars <- intersect(cont_vars, names(df_cara))

df_cara <- df_cara %>%
  mutate(across(all_of(cont_vars), zscale))

carabid_abund_set <- dredge(
  glmmTMB(total_abundance_c ~ height + width + land_use_bin +
            ground_flora_rich + average_height +
            vegetation_length_from_hedgerow + sampling_round +
            (1|pitfall),
          data = df_cara,
          family = nbinom2()),
  rank = "AICc")

carabid_rich_set <- dredge(
  glmmTMB(richness_c ~ height + width + land_use_bin +
            ground_flora_rich + average_height +
            vegetation_length_from_hedgerow + sampling_round +
            (1|pitfall),
          data = df_cara,
          family = compois()),
  rank = "AICc")


#All models average summaries
avg_pitfall_abund <- model.avg(pitfall_abund_set, subset = delta < 2)
summary(avg_pitfall_abund)
sw(avg_pitfall_abund)

avg_hedge_abund <- model.avg(hedge_abund_set, subset = delta < 2)
summary(avg_hedge_abund)
sw(avg_hedge_abund)

avg_pitfall_rich <- model.avg(pitfall_rich_set, subset = delta < 2)
summary(avg_pitfall_rich)
sw(avg_pitfall_rich)

avg_hedge_rich <- model.avg(hedge_rich_set, subset = delta < 2)
summary(avg_hedge_rich)
sw(avg_hedge_rich)

avg_carab_abund <- model.avg(carabid_abund_set, subset = delta < 2)
summary(avg_carab_abund)
sw(avg_carab_abund)

avg_carab_rich <- model.avg(carabid_rich_set, subset = delta < 2)
summary(avg_carab_rich)
sw(avg_carab_rich)

#Helper to extract IRR from model summaries
extract_irr_avg <- function(avg_model, conf.level = 0.95, digits = 3) {
  
  sm <- summary(avg_model)
  
  if (is.null(sm$coefmat.full)) {
    stop("No full-average coefficients found.")
  }
  
  z <- qnorm(1 - (1 - conf.level)/2)
  
  out <- as.data.frame(sm$coefmat.full) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::rename(
      beta = Estimate,
      se = `Std. Error`,
      p.value = `Pr(>|z|)`,
      z.value = `z value`
    ) %>%
    dplyr::mutate(
      IRR = exp(beta),
      conf.low = exp(beta - z * se),
      conf.high = exp(beta + z * se)
    ) %>%
    dplyr::filter(term != "(Intercept)")
  
  #optional formatting (3 sig fig)
  fmt <- function(x) format(signif(x, digits), trim = TRUE, scientific = FALSE)
  
  out %>%
    dplyr::mutate(
      IRR = fmt(IRR),
      conf.low = fmt(conf.low),
      conf.high = fmt(conf.high),
      `95% CI` = paste0(conf.low, "–", conf.high),
      p.value = fmt(p.value),
      z.value = fmt(z.value)
    ) %>%
    dplyr::select(term, IRR, `95% CI`, p.value, z.value)
}

extract_irr_avg(avg_pitfall_abund)
extract_irr_avg(avg_hedge_abund)
extract_irr_avg(avg_pitfall_rich)
extract_irr_avg(avg_hedge_rich)
extract_irr_avg(avg_carab_abund)
extract_irr_avg(avg_carab_rich)


#Top model diagnostics 
top_model_from_dredge <- function(dredge_obj) {
  mod <- get.models(dredge_obj, subset = 1)[[1]]
  return(mod)
}

pitfall_abund_top <- top_model_from_dredge(pitfall_abund_set)
pitfall_rich_top  <- top_model_from_dredge(pitfall_rich_set)
hedge_abund_top   <- top_model_from_dredge(hedge_abund_set)
hedge_rich_top    <- top_model_from_dredge(hedge_rich_set)
carab_abund_top   <- top_model_from_dredge(carabid_abund_set)
carab_rich_top    <- top_model_from_dredge(carabid_rich_set)

top_models <- list(
  pitfall_abund = pitfall_abund_top,
  pitfall_rich = pitfall_rich_top,
  hedge_abund = hedge_abund_top,
  hedge_rich = hedge_rich_top,
  carab_abund = carab_abund_top,
  carab_rich = carab_rich_top)

#Top models all pass diagnostics
diagnose_model_list(top_models, n = 1000)



#Create effects table
effects <- read_excel("MasterDataNew.xlsx", sheet = "EffectsTable")

#Optional: blank out repeated model names for nicer grouping
effects <- effects %>%
  group_by(Model) %>%
  mutate(Model = if_else(row_number() == 1, Model, "")) %>%
  ungroup()

#Create flextable
ft <- flextable(effects) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  align(part = "header", align = "center") %>%
  align(j = c("Effect size (exp β)", "95% CI", "p-value", "Relative weight", "Support"), 
        align = "center", part = "body") %>%
  align(j = c("Model", "Predictor"), 
        align = "left", part = "header") %>%
  align(j = "Predictor", align = "left", part = "body") %>%
  autofit()

model_rows <- which(effects$Model != "")
model_rows <- model_rows[-1]

ft <- ft %>%
  border(i = model_rows,
         border.top = officer::fp_border(color = "grey70", width = 1),
         part = "body")

ft <- ft %>%
  padding(padding.top = 4, padding.bottom = 4, part = "all") %>%  # row spacing
  line_spacing(space = 1.7, part = "all") %>%                     # text spacing
  fontsize(size = 12, part = "all")

ft


#Export to Word if desired
doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "Effects_Table.docx")

############ 5.2 PERMANOVA tables ----
#Load table (made manually from PERMANOVA model outputs)
effects <- read_excel("MasterDataNew.xlsx", sheet = "PERMANOVA") %>%
  mutate(
    `p-value` = case_when(
      is.na(`p-value`) ~ "-",
      `p-value` < 0.001 ~ "-",
      TRUE ~ sprintf("%.3f", as.numeric(`p-value`))),
    `R²` = sprintf("%.3f", as.numeric(`R²`)))

ft <- flextable(effects) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  align(j = c("Sampling method", "Taxonomic level", "Predictor"),
        align = "left", part = "header") %>%
  align(j = c("R²", "p-value"),
        align = "right", part = "header") %>%
  align(j = c("Sampling method", "Taxonomic level", "Predictor"),
        align = "left", part = "body") %>%
  align(j = c("R²", "p-value"),
        align = "right", part = "body") %>%
  autofit()

model_rows <- which(effects$`Sampling method` != "")
model_rows <- model_rows[-1]

#A couple options here based on desired aesthetics
ft <- ft %>%
  border(i = model_rows,
         border.top = officer::fp_border(color = "grey70", width = 1),
         part = "body") %>%
  width(j = c("Sampling method", "Taxonomic level"), width = 1.4) %>%
  width(j = "Predictor", width = 1.2) %>%
  width(j = c("R²", "p-value"), width = 0.9) %>%
  padding(padding.left = 4, padding.right = 4, part = "all")

ft

ft <- ft %>%
  merge_v(j = c("Sampling method", "Taxonomic level")) %>%
  width(j = c("Sampling method", "Taxonomic level"), width = 1.8) %>%
  width(j = "Predictor", width = 2.0) %>%
  width(j = c("R²", "p-value"), width = 1.1) %>%
  padding(padding.top = 2, padding.bottom = 2, part = "all") %>%
  line_spacing(space = 1.1, part = "all") %>%
  fontsize(size = 13, part = "all")
ft

############ 5.3 Mantel tables ----
#Load effects table
effects <- read_excel("MasterDataNew.xlsx", sheet = "Mantel")

ft <- flextable(effects) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  align(j = c("Sampling method", "Taxonomic level"),
        align = "left", part = "header") %>%
  align(j = c("r (Mantel correlation)", "p-value"),
        align = "right", part = "header") %>%
  align(j = c("Sampling method", "Taxonomic level"),
        align = "left", part = "body") %>%
  align(j = c("r (Mantel correlation)", "p-value"),
        align = "right", part = "body") %>%
  autofit()
ft

model_rows <- which(effects$`Sampling method` != "")
model_rows <- model_rows[-1]

#Different format, grouped by sampling method
ft <- ft %>%
  border(i = model_rows,
         border.top = officer::fp_border(color = "grey70", width = 1),
         part = "body") %>%
  width(j = c("Sampling method", "Taxonomic level"), width = 1.4) %>%
  width(j = c("r (Mantel correlation)", "p-value"), width = 1.3) %>%
  padding(padding.left = 4, padding.right = 4, part = "all")

ft <- ft %>%
  merge_v(j = c("Sampling method", "Taxonomic level")) %>%
  width(j = c("Sampling method", "Taxonomic level"), width = 1.8) %>%
  width(j = "r (Mantel correlation)", width = 2.0) %>%
  width(j = c("r (Mantel correlation)"), width = 2) %>%
  width(j = c("p-value"), width = 1.2) %>%
  padding(padding.top = 2, padding.bottom = 2, part = "all") %>%
  line_spacing(space = 1.1, part = "all") %>%
  fontsize(size = 13, part = "all")
ft
