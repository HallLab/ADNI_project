#### WGCNA ANALYSIS ON P180 ####

#### 0. LIBRARIES ####
library(WGCNA)
library(dplyr)

#### 0.5 PARALLEL PROCESSING ####
enableWGCNAThreads()

#### 1. FUNCTIONS ####

read_metabolite_data <- function(metabolite_type='p180'){
     # Read the cleaned metabolite data
     # 
     # Parameters
     # ----------
     # metabolite_type: string of either p180 or nmr 

     # Returns
     # ----------
     # metabolite_data: dataframe of metabolites

     if(metabolite_type=='p180'){
          metabolite_data = read.csv('../results/p180_cleaned.csv')
     } else if(metabolite_type=='nmr'){
          metabolite_data = read.csv('../results/nmr_cleaned.csv')
     } else {
          print('Metabolite type do not recognized')
     }

     return(metabolite_data)     
}

read_qtpad_data <- function(){
     # Read the qt_pad data and apply basic cleaning
     #
     # Returns
     # ----------
     # qtpad_data: dataframe of qt_pad demographic information

     qtpad_data = read.csv('../data/ADNI_adnimerge_20170629_QT-freeze.csv')
     qtpad_data = qtpad_data %>% 
                  subset(VISCODE == 'bl') %>% 
                  select(c('RID','PTGENDER'))
     return(qtpad_data)
}

stratify_by_sex <- function(metabolite_data,
                            qtpad_data){
     # Stratified the metabolite_data by sex based on qtpad_data
     # Qtpad_data must contain PTGENDER column
     #
     # Parameters
     # ----------
     # metabolite_data: dataframe with metabolite concentration values
     # qtpad_data: qtpad dataframe with demographic information
     #
     # Returns
     # ----------
     # stratified_metabolties: list of dataframes with metabolites stratified
     #                         females first
     #

     #Split males and females
     metabolites_males = merge(metabolite_data, qtpad_data) %>% 
                         subset(PTGENDER=='Male') %>% 
                         select(-PTGENDER)
     metabolites_females = merge(metabolite_data, qtpad_data) %>% 
                           subset(PTGENDER=='Female') %>% 
                           select(-PTGENDER)

     stratified_metabolties <- list('females' = metabolites_females,
                                    'males' = metabolites_males)

     return(stratified_metabolties)
}

compare_modules <- function(stratified_metabolties,
                            softPower=6,
                            plotname='sex_module_comparison_p180'){
     # Compare the modules from metabolite data stratified by sex
     # Female data is used as reference
     #
     # Parameters
     # ----------
     # stratified_metabolties: list of dataframe stratified by sex
     # softPower: int with the sof tpower
     # plotname: str with name of plot

     # Returns
     # ----------
     # mp: module preservation information

     # Female modules
     dat1 = stratified_metabolties$females
     wgcna_females = compute_wgcna(dat1,
                                   plotname=NULL,
                                   plot_modules=NULL)

     # Male modules
     dat2 = stratified_metabolties$males
     wgcna_males = compute_wgcna(dat2,
                                 plotname=NULL,
                                 plot_modules=NULL)
     # PLOTS
     filename = paste0('../results/plots/',
                       plotname,
                       '.pdf')
     pdf(file = filename)
     layout(matrix(c(1:4), 4, 1), heights = rep(c(0.8, 0.2), 2))
     # Plot the female dendrogram
     plotDendroAndColors(wgcna_females[[2]], 
                         wgcna_females[[3]],
                         "Female modules",
                         main = "Female gene dendrogram and module colors",
                         dendroLabels = FALSE,
                         setLayout = FALSE,
                         marAll = c(2,10,3,0),
                         cex.colorLabels = 1.4,
                         cex.main = 2,
                         cex.lab = 1.4,
                         cex.axis = 1.4,
                         addGuide = TRUE)
     # Plot the male dendrogram with female module colors
     plotDendroAndColors(wgcna_males[[2]],
                         wgcna_females[[3]],
                         "Female modules",
                         main = "Male gene dendrogram and female module colors",
                         dendroLabels = FALSE,
                         setLayout = FALSE,
                         marAll = c(2,10,3,0),
                         cex.colorLabels = 1.4,
                         cex.main = 2,
                         cex.lab = 1.4,
                         cex.axis = 1.4,
                         addGuide = TRUE)
     dev.off()
     
     dat1 = dat1 %>%
            select(-RID)
     dat2 = dat2 %>%
            select(-RID)
     multiExpr = list(Female = list(data = dat1), 
                      Male = list(data = dat2))
     multiColor = list(Female = wgcna_females[[3]],
                       Male = wgcna_males[[3]])
     mp = modulePreservation(multiExpr, 
                             multiColor, 
                             referenceNetworks = c(1,2),
                             nPermutations = 1000,
                             randomSeed = 1,
                             verbose = 3)
     ref = 1
     test = 2
     statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1],
                      mp$preservation$observed[[ref]][[test]][, -1])
     statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1],
                    mp$preservation$Z[[ref]][[test]][, -1])
     # Compare preservation to quality:
     print(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
           signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
     return(mp)
}

choose_sf_power <- function(metabolite_data,
                            plotname='wgcna_power_plot_p180'){
     # Go through the process of selecting a set of powers, and 
     # plot the results
     #
     # Parameters
     # ----------
     # metabolite_data: dataframe with the metabolites
     # plotname: str with the name of the plots
     #

     data = metabolite_data %>%
            select(-RID)
     powers = c(c(1:10), 
                seq(from = 12, 
                    to=20, 
                    by=2))
     sft = pickSoftThreshold(data, 
                             powerVector = powers, 
                             verbose = 5)
     # Plot the results
     cex1 = 0.9;
     # Scale-free topology fit index as a function of the soft-thresholding power
     filename = paste0('../results/plots/',
                       plotname,
                       '.pdf')
     pdf(file = filename)
     plot(sft$fitIndices[,1], 
          -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
          xlab="Soft Threshold (power)",
          ylab="Scale Free Topology Model Fit,signed R^2",
          type="n",
          main = paste("Scale independence"))
     text(sft$fitIndices[,1],
          -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
          labels=powers,
          cex=cex1,
          col="red");# this line corresponds to using an R^2 cut-off of h
     abline(h=0.90,
            col="red")
     # Mean connectivity as a function of the soft-thresholding power
     plot(sft$fitIndices[,1],
          sft$fitIndices[,5],
          xlab="Soft Threshold (power)",
          ylab="Mean Connectivity",
          type="n",
          main = paste("Mean connectivity"))
     text(sft$fitIndices[,1],
          sft$fitIndices[,5],
          labels=powers,
          cex=cex1,
          col="red")
     dev.off()
 
}

compute_wgcna <- function(metabolite_data,
                          plotname='wgcna_clustering_p180',
                          plot_modules='wgcna_clustering_modules_p180',
                          minModuleSize=5,
                          softPower=6){
     # Compute WGCNA
     # We set minimum module size to 5, based on previous studies
     # Returns list of TOM, tree, colors, eigenmetabolites, and RID
     #
     # Parameters
     # ----------
     # metabolite_data: dataframe with metabolite concentration values
     # plotname: str with name of the clutering plot
     # plot_modules: str with the name of clustering of modules
     # minModuleSize: int with minimun size of module
     # softPower: int with sof power
     # 
     # Returns
     # ----------
     # wgcna: list of TOM, tree, colors, eigenmetabolites, and RID

     # Compute the TOM 
     data = metabolite_data %>%
            select(-RID)

     adjacency = adjacency(data, 
                           power = softPower)
     TOM = TOMsimilarity(adjacency)
     dissTOM = 1-TOM

     # Call the hierarchical clustering function
     prelim_tree = hclust(as.dist(dissTOM),
                          method = "average")
     # Module identification using dynamic tree cut:
     prelim_modules = cutreeDynamic(dendro = prelim_tree,
                                    distM = dissTOM,
                                    deepSplit = 2,
                                    pamRespectsDendro = FALSE,
                                    minClusterSize = minModuleSize)
     # Convert numeric lables into colors
     prelim_colors = labels2colors(prelim_modules)
         
     ## Merging modules that are very similar
     # Calculate eigengenes
     MEList = moduleEigengenes(data,
                               colors = prelim_colors)
     #MEs = MEList$eigengenes
     # Calculate dissimilarity of module eigengenes
     MEDiss = 1-cor(MEList$eigengene);
     # Cluster module eigengenes
     METree = hclust(as.dist(MEDiss),
                     method = "average");
     # The height cut of 0.25 means modules with correlations greater than 0.75
     MEDissThres = 0.25
     
     # Call an automatic merging function
     merge = mergeCloseModules(data,
                               prelim_colors,
                               cutHeight = MEDissThres,
                               verbose = 3)
     # The merged module colors
     final_colors = merge$colors;
     # Eigengenes of the new merged modules:
     final_MEs = merge$newMEs;
     
     print(table(final_colors))
     # Construct numerical labels corresponding to the colors
     #colorOrder = c("grey",
     #               standardColors(50));
     #moduleLabels = match(moduleColors, colorOrder)-1;
     #MEs = mergedMEs;

     #### PLOTS
     if(!is.null(plotname)){
          # Plot dendrogram with initial and final modules
          filename = paste0('../results/plots/',
                            plotname,
                            '.pdf')
          pdf(file = filename)
          plotDendroAndColors(prelim_tree,
                              cbind(prelim_colors, final_colors),
                              c("Dynamic Tree Cut","Merged dynamic"),
                              dendroLabels = FALSE,
                              hang = 0.03,
                              addGuide = TRUE,
                              guideHang = 0.05)
          dev.off()
     }
     
     if(!is.null(plot_modules)){
          filename = paste0('../results/plots/',
                            plot_modules,
                            '.pdf')
          pdf(file = filename)
          plot(METree,
               main = "Clustering of module eigenmetabolites",
               xlab = "",
               sub = "")
          # Plot the cut line into the dendrogram
          abline(h=MEDissThres,
                 col = "red")
          dev.off()
     }

     wgcna = list('TOM' = TOM,
                  'tree' = prelim_tree,
                  'colors' = final_colors,
                  'eigenmetabolites' = final_MEs,
                  'RID' = metabolite_data$RID)

     return(wgcna)
}

plot_heatmap <- function(TOM,
                         tree,
                         colors,
                         plotname='heatmap_p180'){
     # Plot the TOM as a heatmap

     plotDiss = (TOM)
     diag(plotDiss) = NA
     filename = paste0('../results/plots/',
                       plotname,
                       '.pdf')
     pdf(file = filename)
     TOMplot(plotDiss, 
             tree,
             colors,
             main = "Network heatmap plot")
     dev.off()
}

save_wgcna <- function(wgcna,
                       suffix='p180'){
     # Save eigenmetabolites with RID and module colors to 
     # external file
     #
     # Parameters
     # ----------
     # wgcna: list from compute_wgcna
     # suffix: str to add to filenames
     #
     #

     ME_filename = paste0('../results/eigenmetabolites',
                          suffix,
                          '.csv')
     data = paste(wgcna[[4]],
                  wgcna[[5]])
     write.csv(data,
               file=ME_filename,
               row.names=FALSE)

     colors_name = paste0('../results/module_colors',
                          suffix,
                          '.csv')
     write.table(wgcna[[3]],
                 file=colors_name,
                 row.names=FALSE,
                 col.names=FALSE,
                 quote=FALSE)
}

#### 2. DATA ####
print('-----Reading Data-----')
p180  = read_metabolite_data()
nmr   = read_metabolite_data(metabolite_type='nmr')
qtpad = read_qtpad_data()
#Split males and females
print('-----Stratifying Data-----')
stratified_p180 = stratify_by_sex(p180, qtpad)
stratified_nmr = stratify_by_sex(nmr, qtpad)

#### 3. COMPARE MODULE PRESERVATION
print('-----Comparing modules by sex-----')
mp_p180 = compare_modules(stratified_p180)
mp_nmr  = compare_modules(stratified_nmr,
                          plotname='sex_module_comparison_nmr')

#### 4. CHOOSE THE SOFT-THRESHOLDING POWER
# Can't 
print('-----Choosing soft-thresholding power-----')
choose_sf_power(p180)
choose_sf_power(nmr,
                plotname='wgcna_power_plot_nmr')
# Based on the analysis we chose a power of 6

#### 5. NETWORK CONSTRUCTION AND MODULE DETECTION ####
print('-----Network construction and module detection-----')
wgcna_p180 <- compute_wgcna(p180)
wgcna_nmr <- compute_wgcna(nmr,
                           plotname='wgcna_clustering_nmr',
                           plot_modules='wgcna_clustering_modules_p180')

#### 6. HEATMAP VISUALIZATION ####
print('-----Generating visualizations-----')
plot_heatmap(wgcna_p180[[1]],
             wgcna_p180[[2]],
             wgcna_p180[[3]])

plot_heatmap(wgcna_nmr[[1]],
             wgcna_nmr[[2]],
             wgcna_nmr[[3]],
             plotname='heatmap_nmr')

#### 7. SAVING FILES AND EXPORTING ####
print('-----Saving files-----')
save_wgcna(wgcna_p180)
save_wgcna(wgcna_nmr,
           suffix='nmr')

# Export to cytoscape
#modules = 'turquoise'
#in_module = wgcna[[3]] == modules
#
#cyt = exportNetworkToCytoscape(wgcna[[1]][in_module, in_module],
#                         edgeFile = paste("../results/CytoscapeInput-edges-",
#                                          paste(modules,
#                                                collapse="-"),
#                                          ".txt",
#                                          sep=""),
#                         nodeFile = paste("../results/CytoscapeInput-nodes-",
#                                          paste(modules,
#                                                collapse="-"),
#                                          ".txt",
#                                          sep=""),
#                         weighted = TRUE,
#                         threshold = 0.05,
#                         nodeNames = names(p180[in_module]),
#                         nodeAttr = wgcna[[3]][in_module])
#