#### WGCNA ANALYSIS ON P180 ####

#### 0. LIBRARIES ####
library(WGCNA)
library(dplyr)

#### 1. FUNCTIONS ####

choose_sf_power <- function(data,
                            plotname='wgcna_power_plot'){
     # Go through the process of selecting a set of powers, and 
     # plot the results
     #
     # Parameters
     # ----------
     # data: dataframe with the metabolites

     # Call the network topology analysis function
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

compute_TOM <- function(data,
                        plotname='TOM',
                        softPower=6){

     # Compute the TOM and plot the heatmap
     adjacency = adjacency(data, 
                           power = softPower)
     TOM = TOMsimilarity(adjacency)
     if(!is.null(plotname)){
          filename = paste0('../results/plots/',
                       plotname,
                       '.pdf')
          pdf(file=filename)
          heatmap(TOM, Rowv = NA, Colv = NA)
          dev.off()
     }
     
     return(TOM)
}

compute_wgcna <- function(data,
                          plotname='wgcna_clustering',
                          plot_modules='wgcna_clustering_modules',
                          minModuleSize=5,
                          softPower=6){
     # Compute WGCNA
     # We set minimum module size to 5, based on previous studies
     # Returns list of TOM, tree, colors, and eigenmetabolites

     # Compute the TOM 
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
               main = "Clustering of module eigengenes",
               xlab = "",
               sub = "")
          # Plot the cut line into the dendrogram
          abline(h=MEDissThres,
                 col = "red")
          dev.off()
     }

     return(list(TOM, prelim_tree, final_colors, final_MEs))

}

compare_modules <- function(dat1,
                            dat2,
                            softPower=6,
                            plotname='female_male_module_comparison'){
     # Compare the modules from two datasets
     # dat1 is used as reference

     # Female modules
     wgcna_females = compute_wgcna(dat1,
                                   plotname=NULL,
                                   plot_modules=NULL)

     # Male modules
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

#### 2. DATA ####
p180  = read.csv('../results/p180_cleaned.csv')
qtpad = read.csv('../data/ADNI_adnimerge_20170629_QT-freeze.csv')
qtpad = qtpad %>% 
        subset(VISCODE == 'bl') %>% 
        select(c('RID','PTGENDER'))

#Split males and females
p180_males   = merge(p180, qtpad) %>% 
               subset(PTGENDER=='Male') %>% 
               select(-c(RID, PTGENDER))
p180_females = merge(p180, qtpad) %>% 
               subset(PTGENDER=='Female') %>% 
               select(-c(RID, PTGENDER))
p180 = p180 %>% select(-RID)

#### 3. COMPARE MODULE PRESERVATION
mp = compare_modules(p180_females, 
                     p180_males)

#### 4. CHOOSE THE SOFT-THRESHOLDING POWER
choose_sf_power(p180)
# Based on the analysis we chose a power of 6

#### 5. NETWORK CONSTRUCTION AND MODULE DETECTION ####
wgcna <- compute_wgcna(p180)

# Save the eigenmetabolites
write.csv(wgcna[[4]],
          file="../results/eigenmetabolites.csv",
          row.names=FALSE)
          
write.table(wgcna[[3]],
           file="../results/module_colors.csv",
           row.names=FALSE,
           col.names=FALSE,
           quote=FALSE)
# Save module colors and labels for use in subsequent parts
#save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")
