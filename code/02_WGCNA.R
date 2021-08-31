#### WGCNA ANALYSIS ON METABOLOMICS ####

#### 1. LIBRARIES ####
source("wgcna_functions.R")

#### 2. DATA ####
print("-----Reading Data-----")
p180  <- read_metabolite_data()
nmr   <- read_metabolite_data(metabolite_type = "nmr")
qtpad <- read_qtpad_data()

#Split males and females
print("-----Stratifying Data-----")
stratified_p180 <- stratify_by_sex(p180, qtpad)
stratified_nmr <- stratify_by_sex(nmr, qtpad)

#### 3. CHOOSE THE SOFT-THRESHOLDING POWER
print("-----Choosing soft-thresholding power-----")
choose_sf_power(p180,
                plotname = "wgcna_power_p180")
choose_sf_power(nmr,
                plotname = "wgcna_power_nmr",
                suffix = "NMR")
# Based on the analysis we chose a power of 9 for p180,
# and 11 for nmr

#### 4 CHECK SCALE-FREE TOPOLOGY
check_scale_free(p180,
                 soft_power = 9)
check_scale_free(nmr,
                 plotname = "scale_free_check_nmr",
                 soft_power = 11)

#### 4. COMPARE MODULE PRESERVATION
#print("-----Comparing modules by sex-----")
#mp_p180 <- compare_modules(stratified_p180,
#                           soft_power = 7)
#mp_nmr  <- compare_modules(stratified_nmr,
#                           soft_power = 19,
#                           plotname = "sex_module_comparison_nmr")

#### 5. NETWORK CONSTRUCTION AND MODULE DETECTION ####
print("-----Network construction and module detection-----")
wgcna_p180 <- compute_wgcna(p180,
                            soft_power = 9)
wgcna_nmr <- compute_wgcna(nmr,
                           soft_power = 11,
                           plotname = "wgcna_dendrogram_nmr",
                           plot_modules = "wgcna_module_tree_nmr",
                           min_module_size = 10)

#### 6. VISUALIZATION ####
print("-----Generating visualizations-----")
plot_heatmap(wgcna_p180)
plot_heatmap(wgcna_nmr,
             suffix = "nmr")

#### 7. SAVING FILES AND EXPORTING ####
print("-----Saving files-----")
save_wgcna(wgcna_p180)
save_wgcna(wgcna_nmr,
           suffix = "nmr")

#### 8. PRINT SESSION ####
sessionInfo()
