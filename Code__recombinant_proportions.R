###############################################################################
### HIV MOL EPI RECOMB 1990-2024 - Data analysis - Recombinant proportions
### Ayisha Khalid
### R Version 4.4.1

### Created: 17 January 2026
###############################################################################


### Set up ####################################################################

## Libraries
library(tidyverse)
library(readxl)
library(writexl)

options(scipen=999) 



### Data #####################################################################

# Country-specific annual proportions of HIV-1 recombinants
# Sample data for East Africa
recomb_country_data <- read_xlsx("~/Desktop/HIV Analysis/Code/GitHub/Paper 2/sample_data_EA_recomb_country_dist.xlsx")


### Distribution of recombinants ##############################################

#### Country-level ####

# Recombinant roportions
recomb_country <- recomb_country_data %>%
  group_by(year_period, region, country) %>%
  mutate(recomb_total = sum(variant_n, na.rm = TRUE), # total recombinants in country
         plhiv_recomb_total_country = sum(plhiv_variant_country)) %>% # total PLHIV with recombinants in country
  mutate(recomb_pct = (variant_n/recomb_total)*100, na.rm = TRUE) # proportion of recombinant among all recombinant infections

colnames(recomb_country)
recomb_country <- recomb_country %>% 
  select(region, country, year_period, variant, variant_n, recomb_total, 
         recomb_pct, plhiv_variant_country, plhiv_recomb_total_country)

# Variance
recomb_country <- recomb_country %>%
  mutate(p_i = variant_n/recomb_total, 
         var_country = p_i * (1 - p_i) / recomb_total) 



#### Region-level ####

## Recombinant proportions
recomb_region <- recomb_country %>%
  group_by(year_period, region, variant) %>%
  mutate(plhiv_variant_region = sum(plhiv_variant_country, na.rm = TRUE)) 

recomb_region <- recomb_region %>% 
  group_by(year_period, region) %>%
  mutate(plhiv_region_data = sum(plhiv_variant_country, na.rm = TRUE )) %>% ungroup() # total PLHIV in countries with available data

recomb_region <- recomb_region %>%
  group_by(year_period, region, variant) %>%
  summarise(
    plhiv_variant_region = unique(plhiv_variant_region), 
    plhiv_region_data = unique(plhiv_region_data), 
    variant_pct_region = (plhiv_variant_region / plhiv_region_data)*100,
    
    # Confidence intervals
    var_reg = sum( var_country * ( (plhiv_recomb_total_country / unique(plhiv_region_data))^2 ), na.rm = TRUE ),
    se_reg  = sqrt(var_reg),
    ci_low     = pmax(0, (plhiv_variant_region / plhiv_region_data) - 1.96 * se_reg),
    ci_high    = pmin(1, (plhiv_variant_region / plhiv_region_data) + 1.96 * se_reg),
    ci_low_pct_reg  = 100 * ci_low,
    ci_high_pct_reg = 100 * ci_high,
    .groups = "drop")


## Number of PLHIV with recombinants, scaled to entire region including countries without available data
plhiv_recomb_region_adj_data <- recomb_country_data %>%
  distinct(year_period, region, plhiv_recomb_region_unaids)

recomb_region <- left_join(recomb_region, plhiv_recomb_region_adj_data)  
recomb_region <- recomb_region %>%
  mutate(plhiv_variant_region_adj = (variant_pct_region/100)*plhiv_recomb_region_unaids)



# #### Global ####

# Example code below

# recomb_global <- recomb_region %>%
#   group_by(year_period, variant) %>%
#   summarise(
#     ### Number of PLHIV
#     ## Sum region data (ADJUSTED) to get number of PLHIV by variant globally
#     plhiv_variant_global = sum(plhiv_variant_region_adj, na.rm = TRUE),
#     ### Variant proportions
#     ## First, sum to get total number of PLHIV
#     plhiv_global = sum(plhiv_recomb_region_unaids, na.rm = TRUE),
#     ## Second, divide to get proportion (global weighted proportion)
#     variant_pct_global = (plhiv_variant_global / plhiv_global) * 100,
#     
#     # Confidence intervals 
#     var_global = sum( var_reg * ( (plhiv_recomb_region_unaids / unique(plhiv_global))^2 ), na.rm = TRUE ),
#     se_global  = sqrt(var_global),
#     ci_low     = pmax(0, (plhiv_variant_global / plhiv_global) - 1.96 * se_global),
#     ci_high    = pmin(1, (plhiv_variant_global / plhiv_global) + 1.96 * se_global),
#     ci_low_pct_global  = 100 * ci_low,
#     ci_high_pct_global = 100 * ci_high,
#     .groups = "drop"
#   )




### Distribution of other CRFs ################################################

### Country-level ####

# Subset
othcrfs_country <- recomb_country_data %>% 
  filter(variant != "crf01_ae") %>% 
  filter(variant != "crf02_ag") %>%
  filter(variant != "crf07_bc") %>%
  filter(variant != "urfs") %>%
  filter(variant != "unspec_recombs") 

# Other CRF proportions
othcrfs_country <- othcrfs_country %>%
  group_by(year_period, region, country) %>%
  mutate(othcrfs_total = sum(variant_n), 
         plhiv_othcrfs_total_country = sum(plhiv_variant_country)) %>%
  mutate(othcrfs_pct = (variant_n/othcrfs_total)*100) 
colnames(othcrfs_country)
othcrfs_country <- othcrfs_country %>% select(region, country, year_period, variant, variant_n, othcrfs_total, othcrfs_pct, plhiv_variant_country, plhiv_othcrfs_total_country)

# Variance
othcrfs_country <- othcrfs_country %>%
  mutate(p_i = variant_n/othcrfs_total, 
         var_country = p_i * (1 - p_i) / othcrfs_total) 



#### Region-level ####

## Other CRF proportions
othcrfs_region <- othcrfs_country %>%
  group_by(year_period, region, variant) %>%
  mutate(plhiv_variant_region = sum(plhiv_variant_country, na.rm = TRUE)) 

othcrfs_region <- othcrfs_region %>% 
  group_by(year_period, region) %>%
  mutate(plhiv_region_data = sum(plhiv_variant_country, na.rm = TRUE )) %>% ungroup() 

othcrfs_region <- othcrfs_region %>%
  group_by(year_period, region, variant) %>%
  summarise(
    plhiv_variant_region = unique(plhiv_variant_region),  
    plhiv_region_data = unique(plhiv_region_data), 
    variant_pct_region = (plhiv_variant_region / plhiv_region_data)*100,
    
    # Confidence intervals
    var_reg = sum( var_country * ( (plhiv_othcrfs_total_country / unique(plhiv_region_data))^2 ), na.rm = TRUE ),
    se_reg  = sqrt(var_reg),
    ci_low     = pmax(0, (plhiv_variant_region / plhiv_region_data) - 1.96 * se_reg),
    ci_high    = pmin(1, (plhiv_variant_region / plhiv_region_data) + 1.96 * se_reg),
    ci_low_pct_reg  = 100 * ci_low,
    ci_high_pct_reg = 100 * ci_high,
    .groups = "drop")


## Number of PLHIV with other CRFs, scaled to entire region including countries without available data
plhiv_othcrfs_region_adj_data <- recomb_country_data %>%
  distinct(year_period, region, plhiv_othcrfs_region_unaids)

othcrfs_region <- left_join(othcrfs_region, plhiv_othcrfs_region_adj_data)  
othcrfs_region <- othcrfs_region %>%
  mutate(plhiv_variant_region_adj = (variant_pct_region/100)*plhiv_othcrfs_region_unaids)



#### Global-level ####

# # Example code below

# othcrfs_global <- othcrfs_region %>%
#   group_by(year_period, variant) %>%
#   summarise(
#     ### Number of PLHIV
#     ## Sum region data (ADJUSTED) to get number of PLHIV by variant globally
#     plhiv_variant_global = sum(plhiv_variant_region_adj, na.rm = TRUE),
#     ### Variant proportions
#     ## First, sum to get total number of PLHIV
#     plhiv_global = sum(plhiv_othcrfs_region_unaids, na.rm = TRUE),
#     ## Second, divide to get proportion (global weighted proportion)
#     variant_pct_global = (plhiv_variant_global / plhiv_global) * 100,
#     
#     # Confidence intervals 
#     var_global = sum( var_reg * ( (plhiv_othcrfs_region_unaids / unique(plhiv_global))^2 ), na.rm = TRUE ),
#     se_global  = sqrt(var_global),
#     ci_low     = pmax(0, (plhiv_variant_global / plhiv_global) - 1.96 * se_global),
#     ci_high    = pmin(1, (plhiv_variant_global / plhiv_global) + 1.96 * se_global),
#     ci_low_pct_global  = 100 * ci_low,
#     ci_high_pct_global = 100 * ci_high,
#     .groups = "drop"
#   )



