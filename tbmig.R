## Analysis code for TBMigrants by Schwalb et al. 2023
## Distributed under CC BY 4.0

# Packages ====
pacman::p_load(rio, here, tidyverse, janitor, skimr, RColorBrewer, 
               survival, survminer, ggplot2, patchwork)

# Data ====
tbmig <- import(here("data","tbmigrants.xlsx"))

# Data curation ====
tbmig <- clean_names(tbmig) # Janitor function to clean column names
skim(tbmig) # Skim through database to curate

tbmig <- tbmig %>%
  mutate(id = as.numeric(str_remove_all(code, "S")), # Remove "S" character and change name to "ID"
         agegp = replace(age, age == 11, 1), # Clean extreme value and create age group variable
         time_a = as.numeric(str_replace_all(time_a,"[^-.0-9]","")), # Clean Time_A variable 
         time_b = as.numeric(str_replace_all(time_b,"[^-.0-9]","")), # Clean Time_B variable
         travel = replace(travel, id == 208, 1),# Individual with travel time registered
         travel = as.numeric(str_replace_all(travel,"[^-.0-9]","")), # Clean Travel variable
         time_x = as.numeric(str_replace_all(time_x,"[^-.0-9]","")), # Clean Time_X variable 
         time_y = as.numeric(str_replace_all(time_y,"[^-.0-9]","")), # Clean Time_Y variable
         time_z = if_else(time_x < time_a, time_x, time_a, time_a), # Create Time_Z variable, for most recent time
         contact = as.numeric(str_detect(hcw_1st_met, "contact")), # Create contact variable
         actscr = active_1_or_passive_0, # Change name of variable
         etoh = as.numeric(replace(et_oh, et_oh == "poss", 1)), # Clean values and change to numeric
         hiv = as.numeric(str_replace_all(hiv,"[^-.0-9]","")), # Clean HIV variable 
         uz = as.numeric(str_replace_all(uz,"[^-.0-9]","")), # Clean UZ variable
         mlz = as.numeric(str_replace_all(mlz,"[^-.0-9]","")), # Clean MLZ variable
         inh_res = as.numeric(str_detect(epi_contact, "outbreak")), # Create INH-resistant outbreak variable
         inh_res = (replace(inh_res, is.na(inh_res), 0)), # Replace NAs 
         tbdx = 1) %>% # Create dummy variable for TB diagnosis
  select(id,tbdx,sex,agegp,time_a,time_b,travel,time_x,time_y,time_z,contact,actscr,etoh,dm,nfa,hiv,cavitation,uz,mlz,inh_res) %>% # Select variables
  mutate(across(.cols=c(sex,agegp,travel,actscr,contact,etoh,dm,nfa,hiv,cavitation,uz,mlz,inh_res), .fns=as.factor)) %>% # Coerce into factors
  arrange(id) # Order by ID
export(tbmig, here("data","tbmig - clean.xlsx"))

tbmig <- tbmig %>% # (n = 369)
  filter(!(is.na(travel))) %>% # Individuals without info on travel (n=11)
  filter(!(is.na(time_x) & is.na(time_y) & travel == 1)) %>% # Individuals with travel with no specified dates (n=10)
  filter(!(is.na((time_a) & is.na((time_b))))) %>% # Individuals without data on time since migration (n=4)
  filter(!(inh_res == 1)) # Removing those part of INH-resistant outbreak (n=24)

# Descriptive analysis ====
skim(tbmig) # Skim through database after data curation
# 356 individuals, all S+ PTB with travel and follow-up information

# No missing values from following variables
tbmig %>% tabyl(sex)
tbmig %>% tabyl(agegp)
tbmig %>% tabyl(actscr)
tbmig %>% tabyl(contact)
tbmig %>% tabyl(dm)
tbmig %>% tabyl(etoh)
tbmig %>% tabyl(nfa)
tbmig %>% tabyl(hiv)
tbmig %>% tabyl(travel)
tbmig %>% tabyl(cavitation)
tbmig %>% tabyl(inh_res)
tbmig %>% tabyl(travel, actscr)

# Sub-group descriptive analysis ====
tbmigsub <- tbmig %>% # (n = 320)
  filter(time_a >= 2*365) # Cases after 2 years (n = 249)

tbmigsub %>% tabyl(travel)

# Survival analysis ====
sum(tbmig$time_a, na.rm = TRUE)/length(tbmig)/365.25 # Total follow-up (person-years)

# Survival from migration
survmig <- Surv(time = tbmig$time_a, event = tbmig$tbdx)
survmigfit <- survfit(survmig ~ 1)
survmigfitq <- quantile(survmigfit, probs = c(0.5, 0.8, 1.0), conf.int = TRUE)
survmigfitq$quantile/365.25 # Best estimate
survmigfitq$lower/365.25 # Lower bound
survmigfitq$upper/365.25 # Upper bound

ggsurvplot(survmigfit, data.frame(time=tbmig$time_a, tbmig$tbdx), fun = "event", surv.median.line = "hv")

# Survival from travel
survtrav <- Surv(time = tbmig$time_x, event = tbmig$tbdx)
survtravfit <- survfit(survtrav ~ 1)
survtravfitq <- quantile(survtravfit, probs = c(0.5, 0.8, 1.0), conf.int = TRUE)
survtravfitq$quantile/365.25 # Best estimate
survtravfitq$lower/365.25 # Lower bound
survtravfitq$upper/365.25 # Upper bound

ggsurvplot(survtravfit, data.frame(time=tbmig$time_x, tbmig$tbdx), fun = "event", surv.median.line = "hv")

fits <- list(TRAV = survtravfit,  MIG = survmigfit)

km <- ggsurvplot(fits, ggtheme = theme_classic(), palette = c("#FC8D62","#7570B3"), conf.int = TRUE,
           xlab = "Years", ylab = "TB diagnosis", fun = "event", combine = TRUE,
           legend = "bottom", legend.title = "Follow-up",  legend.labs = c("travel","migration"),
           risk.table = "abs_pct", fontsize = 3, ylim = c(0, 1), surv.scale = c("percent"),
           xscale = 365.25, break.x.by = 365.25*2, xlim = c(0,365.25*20), axes.offset = FALSE)

kmplot <- km$plot + 
  geom_segment(aes(x = 2363, y = 0.5, xend = 0, yend = 0.5), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 2363, y = 0.5, xend = 2363, yend = 0), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 231, y = 0.5, xend = 231, yend = 0), linetype = "dashed", color = "grey") + 
  geom_segment(aes(x = 5904, y = 0.8, xend = 0, yend = 0.8), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 5904, y = 0.8, xend = 5904, yend = 0), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 902, y = 0.8, xend = 902, yend = 0), linetype = "dashed", color = "grey")
  
layout <- "A
           A
           A
           A
           A
           A
           B"

tiff("Figure.tiff",  width = 7, height = 5, units = 'in', res = 200) 
kmplot / km$table + plot_layout(design = layout)
dev.off()

# Sensitivity analysis ====
# No active case finding
tbmigsa <- tbmig %>% # (n = 320)
  filter(!(actscr == 1)) # Remove individuals detected through active screening (n=26)

sum(tbmigsa$time_a, na.rm = TRUE)/length(tbmigsa)/365.25 # Total follow-up (person-years)

# Survival from migration
survmig <- Surv(time = tbmigsa$time_a, event = tbmigsa$tbdx)
survmigfit <- survfit(survmig ~ 1)
survmigfitq <- quantile(survmigfit, probs = c(0.5, 0.8, 1.0), conf.int = TRUE)
survmigfitq$quantile/365.25 # Best estimate
survmigfitq$lower/365.25 # Lower bound
survmigfitq$upper/365.25 # Upper bound

ggsurvplot(survmigfit, data.frame(time=tbmigsa$time_a, tbmigsa$tbdx), fun = "event", surv.median.line = "hv")

# Survival from travel
survtrav <- Surv(time = tbmigsa$time_x, event = tbmigsa$tbdx)
survtravfit <- survfit(survtrav ~ 1)
survtravfitq <- quantile(survtravfit, probs = c(0.5, 0.8, 1.0), conf.int = TRUE)
survtravfitq$quantile/365.25 # Best estimate
survtravfitq$lower/365.25 # Lower bound
survtravfitq$upper/365.25 # Upper bound

ggsurvplot(survtravfit, data.frame(time=tbmigsa$time_x, tbmigsa$tbdx), fun = "event", surv.median.line = "hv")

fits <- list(TRAV = survtravfit, MIG = survmigfit)

km <- ggsurvplot(fits, ggtheme = theme_classic(), palette = c("#D95F02","#7570B3"), conf.int = TRUE,
                 xlab = "Years", ylab = "TB diagnosis", fun = "event", combine = TRUE,
                 legend = "bottom", legend.title = "Follow-up",  legend.labs = c("travel", "migration"),
                 risk.table = "abs_pct", fontsize = 3, ylim = c(0, 1), surv.scale = c("percent"),
                 xscale = 365.25, break.x.by = 365.25*2, xlim = c(0,365.25*20), axes.offset = FALSE)

kmplot <- km$plot + 
  geom_segment(aes(x = 2420, y = 0.5, xend = 0, yend = 0.5), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 2420, y = 0.5, xend = 2420, yend = 0), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 277, y = 0.5, xend = 277, yend = 0), linetype = "dashed", color = "grey") + 
  geom_segment(aes(x = 5963, y = 0.8, xend = 0, yend = 0.8), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 5963, y = 0.8, xend = 5963, yend = 0), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 991, y = 0.8, xend = 991, yend = 0), linetype = "dashed", color = "grey")

kmplot / km$table + plot_layout(design = layout)
  
tiff("SFigure 1.tiff",  width = 7, height = 5, units = 'in', res = 200) 
kmplot / km$table + plot_layout(design = layout)
dev.off()

# Follow-up since migration for travellers
tbmigtrav <- tbmig %>% # (n = 320)
  filter(travel == 1) # Individuals who reported travel (n = 134)

sum(tbmigtrav$time_a, na.rm = TRUE)/length(tbmigtrav)/365.25 # Total follow-up (person-years)

# Survival from migration
survmig <- Surv(time = tbmigtrav$time_a, event = tbmigtrav$tbdx)
survmigfit <- survfit(survmig ~ 1)
survmigfitq <- quantile(survmigfit, probs = c(0.5, 0.8, 1.0), conf.int = TRUE)
survmigfitq$quantile/365.25 # Best estimate
survmigfitq$lower/365.25 # Lower bound
survmigfitq$upper/365.25 # Upper bound

ggsurvplot(survmigfit, data.frame(time=tbmigtrav$time_a, tbmigtrav$tbdx), fun = "event", surv.median.line = "hv")

# Follow-up since travel (excluding non-travellers)
survtrav <- Surv(time = tbmigtrav$time_x, event = tbmigtrav$tbdx)
survtravfit <- survfit(survtrav ~ 1)
survtravfitq <- quantile(survtravfit, probs = c(0.5, 0.8, 1.0), conf.int = TRUE)
survtravfitq$quantile/365.25 # Best estimate
survtravfitq$lower/365.25 # Lower bound
survtravfitq$upper/365.25 # Upper bound

ggsurvplot(survtravfit, data.frame(time=tbmigtrav$time_x, tbmigtrav$tbdx), fun = "event", surv.median.line = "hv")

fits <- list(TRAV = survtravfit, MIG = survmigfit)

km <- ggsurvplot(fits, ggtheme = theme_classic(), palette = c("#D95F02","#7570B3"), conf.int = TRUE,
                 xlab = "Years", ylab = "TB diagnosis", fun = "event", combine = TRUE,
                 legend = "bottom", legend.title = "Follow-up",  legend.labs = c("travel", "migration"),
                 risk.table = "abs_pct", fontsize = 3, ylim = c(0, 1), surv.scale = c("percent"),
                 xscale = 365.25, break.x.by = 365.25*2, xlim = c(0,365.25*20), axes.offset = FALSE)

kmplot <- km$plot + 
  geom_segment(aes(x = 4241, y = 0.5, xend = 0, yend = 0.5), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 4241, y = 0.5, xend = 4241, yend = 0), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 231, y = 0.5, xend = 231, yend = 0), linetype = "dashed", color = "grey") + 
  geom_segment(aes(x = 8202, y = 0.8, xend = 0, yend = 0.8), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 8202, y = 0.8, xend = 8202, yend = 0), linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 902, y = 0.8, xend = 902, yend = 0), linetype = "dashed", color = "grey")

kmplot / km$table + plot_layout(design = layout)

tiff("SFigure 2.tiff",  width = 7, height = 5, units = 'in', res = 200) 
kmplot / km$table + plot_layout(design = layout)
dev.off()
