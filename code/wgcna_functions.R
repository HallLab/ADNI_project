#### 0. SETTINGS AND LIBRARIES ####
library(WGCNA)
library(dplyr)
library(CoExpNets)
enableWGCNAThreads()

#### 1. FUNCTIONS ####

read_metabolite_data <- function(metabolite_type="p180") {
     # Read the cleaned metabolite data
     #
     # Parameters
     # ----------
     # metabolite_type: string of either p180 or nmr
     # Returns
     # ----------
     # metabolite_data: dataframe of metabolites
     #

     if (metabolite_type == "p180") {
          metabolite_data <- read.csv("../results/p180_cleaned.csv")
     } else if (metabolite_type == "nmr") {
          metabolite_data <- read.csv("../results/nmr_cleaned.csv")
     } else {
          print("Metabolite type not recognized")
     }
     return(metabolite_data)
}

read_qtpad_data <- function() {
     # Read the qt_pad data and apply basic cleaning
     #
     # Returns
     # ----------
     # qtpad_data: dataframe of qt_pad demographic information
     #

     qtpad_data <- read.csv("../data/ADNI_adnimerge_20170629_QT-freeze.csv")
     qtpad_data <- qtpad_data %>%
                   subset(VISCODE == "bl") %>%
                   select(c("RID", "PTGENDER"))
     return(qtpad_data)
}

stratify_by_sex <- function(metabolite_data,
                            qtpad_data) {
     # Stratifiy the metabolite_data by sex based on qtpad_data.
     # Normalizes each column after the stratification.
     # Qtpad_data must contain PTGENDER column
     #
     # Parameters
     # ----------
     # metabolite_data: dataframe with metabolite concentration values
     # qtpad_data: qtpad dataframe with demographic information
     #
     # Returns
     # ----------
     # stratified_metabolites: list of dataframes with metabolites stratified
     #                         females first
     #

     #Split males and females
     n_cols <- ncol(metabolite_data)
     metabolites_males <- merge(metabolite_data, qtpad_data) %>%
                          subset(PTGENDER == "Male") %>%
                          select(-PTGENDER)
     metabolites_males[, 2:n_cols] <- metabolites_males[, 2:n_cols] %>%
                                      scale()

     metabolites_females <- merge(metabolite_data, qtpad_data) %>%
                            subset(PTGENDER == "Female") %>%
                            select(-PTGENDER)
     metabolites_females[, 2:n_cols] <- metabolites_females[, 2:n_cols] %>%
                                        scale()

     stratified_metabolites <- list("females" = metabolites_females,
                                    "males" = metabolites_males)

     return(stratified_metabolites)
}

compare_modules <- function(stratified_metabolites,
                            soft_power=6,
                            plotname="sex_module_comparison_p180") {
     # Compare the modules from metabolite data stratified by sex
     # Female data is used as reference
     #
     # Parameters
     # ----------
     # stratified_metabolites: list of dataframe stratified by sex
     # soft_power: int with the soft power
     # plotname: str with name of plot
     #
     # Returns
     # ----------
     # mp: module preservation information
     #

     # Female modules
     dat1 <- stratified_metabolites$females
     wgcna_females <- compute_wgcna(dat1,
                                    plotname = NULL,
                                    plot_modules = NULL,
                                    soft_power = soft_power)

     # Male modules
     dat2 <- stratified_metabolites$males
     wgcna_males <- compute_wgcna(dat2,
                                  plotname = NULL,
                                  plot_modules = NULL,
                                  soft_power = soft_power)
     # PLOTS
     filename <- paste0("../results/plots/",
                       plotname,
                       ".pdf")
     pdf(file = filename)
     layout(matrix(c(1:4), 4, 1), heights = rep(c(0.8, 0.2), 2))
     # Plot the female dendrogram
     plotDendroAndColors(wgcna_females[[2]],
                         wgcna_females[[3]],
                         "Female modules",
                         main = "Female gene dendrogram and module colors",
                         dendroLabels = FALSE,
                         setLayout = FALSE,
                         marAll = c(2, 10, 3, 0),
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
                         marAll = c(2, 10, 3, 0),
                         cex.colorLabels = 1.4,
                         cex.main = 2,
                         cex.lab = 1.4,
                         cex.axis = 1.4,
                         addGuide = TRUE)
     dev.off()

     dat1 <- dat1 %>%
             select(-RID)
     dat2 <- dat2 %>%
             select(-RID)
     multi_expr <- list(Female = list(data = dat1),
                        Male = list(data = dat2))
     multi_color <- list(Female = wgcna_females[[3]],
                         Male = wgcna_males[[3]])
     mp <- modulePreservation(multi_expr,
                              multi_color,
                              referenceNetworks = c(1, 2),
                              nPermutations = 500,
                              randomSeed = 1,
                              verbose = 3)
     ref <- 1
     test <- 2
     stats_obs <- cbind(mp$quality$observed[[ref]][[test]][, -1],
                        mp$preservation$observed[[ref]][[test]][, -1])
     stats_z <- cbind(mp$quality$Z[[ref]][[test]][, -1],
                      mp$preservation$Z[[ref]][[test]][, -1])
     # Compare preservation to quality:
     print(cbind(stats_obs[, c("medianRank.pres", "medianRank.qual")],
           signif(stats_z[, c("Zsummary.pres", "Zsummary.qual")], 2)))
     return(mp)
}

choose_sf_power <- function(metabolites,
                            plotname="wgcna_power_p180",
                            stratified=NULL) {
     # Go through the process of selecting a set of powers, and
     # plot the results
     #
     # Parameters
     # ----------
     # metabolites: dataframe with metabolite concentration
     #              values stratified by sex or not
     # plotname: str with the name of the plots
     # stratified: whether the metabolite input is stratified or not
     #             (the function will try to guess, use it only if guessing
     #             is not correct)
     #

     if (class(metabolites) == "list") {
          stratified <- TRUE
     } else if (class(metabolites) == "data.frame") {
          stratified <- FALSE
     } else {
          print("The class of metabolites is not correct")
          stop()
     }

     powers <- c(2:30)

     get_power_info <- function(data,
                                powers) {
          data <- data %>%
                  select(-RID)
          sft <- pickSoftThreshold(data,
                                   RsquaredCut = 0.70,
                                   powerVector = powers)
          return(sft)
     }

     plot_sf_power <- function(power_table,
                               filename,
                               powers) {
          cex1 <- 0.9
          labels <- c("females",
                      "males")
          colors <- c("red",
                      "black")

          if (class(power_table) == "list") {
               x <- power_table[[2]][, 1]
               y <- -sign(power_table[[1]][, 3]) *
                          power_table[[1]][, 2]
               y_males <- -sign(power_table[[2]][, 3]) *
                                power_table[[2]][, 2]

          } else if (class(power_table) == "data.frame") {
               x <- power_table[, 1]
               y <- -sign(power_table[, 3]) *
                          power_table[, 2]
          }

          y_min <- min(y) - 0.2

          pdf(file = filename)
          plot(x,
               y,
               xlab = "Soft Threshold (power)",
               ylab = "Scale Free Topology Model Fit,signed R^2",
               type = "n",
               main = paste("Scale independence"),
               ylim = c(y_min, 1))
          abline(h = 0.8,
                 col = "red")
          text(x,
               y,
               labels = powers,
               cex = cex1,
               col = colors[1])
          if (class(power_table) == "list") {
               text(x,
                    y_males,
                    labels = powers,
                    cex = cex1,
                    col = colors[2])
               legend("bottomright",
                      legend = labels,
                      col = colors,
                      pch = 20)
          }

          if (class(power_table) == "list") {
               y <- power_table[[1]][, 5]
               y_males <- power_table[[2]][, 5]

          } else if (class(power_table) == "data.frame") {
               y <- power_table[, 5]
          }

          plot(x,
               y,
               xlab = "Soft Threshold (power)",
               ylab = "Mean Connectivity",
               type = "n",
               main = paste("Mean connectivity"))
          text(x,
               y,
               labels = powers,
               cex = cex1,
               col = colors[1])
          if (class(power_table) == "list") {
               text(x,
                    y_males,
                    labels = powers,
                    cex = cex1,
                    col = colors[2])
               legend("topright",
                      legend = labels,
                      col = colors,
                      pch = 20)
          }
          dev.off()
     }

     if (stratified == TRUE) {
          n_sets <- length(metabolites)
          power_tables <- vector(mode = "list",
                                 length = n_sets)
          for (i in 1:n_sets) {
               power_tables[[i]] <- get_power_info(metabolites[[i]],
                                                   powers)
          }

          # Plot the results
          filename <- paste0("../results/plots/",
                             plotname,
                             ".pdf")
          plot_sf_power(power_tables[[2]],
                        filename,
                        powers)
     } else if (stratified == FALSE) {
          power_tables <- get_power_info(metabolites,
                                         powers)
          # Plot the results
          filename <- paste0("../results/plots/",
                             plotname,
                             ".pdf")
          plot_sf_power(power_tables[[2]],
                        filename,
                        powers)
     }
     statement <- paste0("The best power is ",
                         as.character(power_tables[[1]]))
     print(statement)
}

check_scale_free <- function(metabolites,
                             plotname="scale_free_check_p180",
                             soft_power=6) {
     # Check the scale free topology criteria by plotting
     # the network connectivities
     #
     # Parameters
     # ----------
     # metabolites: dataframe with metabolite concentration
     #              values stratified by sex or not
     # plotname: str with the name of the plots
     # soft_power: selected soft-threshold to use
     #

     filename <- paste0("../results/plots/",
                        plotname,
                        ".pdf")
     data <- metabolites %>%
             select(-RID)
     ADJ1 <- abs(cor(data,
                     use = "p"))^soft_power
     k <- as.vector(apply(ADJ1,
                          2,
                          sum,
                          na.rm = T))
     pdf(file = filename)
     hist(k)
     scaleFreePlot(k,
                   main = "Check scale free topology\n")
     dev.off()
}

compute_wgcna <- function(metabolites,
                          plotname="wgcna_dendrogram_p180",
                          plot_modules="wgcna_module_tree_p180",
                          min_module_size=5,
                          soft_power=6) {
     # Compute WGCNA or consensus WGCNA depending on input
     # We set minimum module size to 5, based on previous studies
     # Returns TOM, tree, colors, eigenmetabolites, and RID
     #
     # Parameters
     # ----------
     # metabolites: dataframe with metabolite concentration
     #              values stratified by sex or not
     # plotname: str with name of the dendrogram plot
     # plot_modules: str with the name of the module tree
     # minModuleSize: int with minimun size of module
     # softPower: int with sof power
     #
     # Returns
     # ----------
     # wgcna: TOM, tree, colors, eigenmetabolites, and RID
     #
     if (class(metabolites) == "list") {
          # Generate multiExpr dataset for WGCNA
          n_sets <- length(metabolites)
          multi_expr <- vector(mode = "list",
                               length = n_sets)
          for (i in 1:n_sets) {
               data <- metabolites[[i]] %>%
                       select(-RID)
               multi_expr[[i]] <- list(data = data)
          }

          # Compute the TOM
          n_metabolites <- length(multi_expr[[1]]$data)
          tom_list   <- vector(mode = "list",
                               length = n_sets)
          for (i in 1:n_sets) {
             adjacency <- adjacency(multi_expr[[i]]$data,
                                    power = soft_power)
             tom_list[[i]] <- TOMsimilarity(adjacency)
          }

          # Scale the TOM for compatibility
          scale_percentile <- 0.95
          set.seed(12345)
          n_samples <- as.integer(1 / (1 - scale_percentile) * 400)
          scale_sample <- sample(n_metabolites * (n_metabolites - 1) / 2,
                                 size = n_samples)
          tom_scaling_samples <- list()
          scale_quant <- rep(1, n_sets)
          scale_powers <- rep(1, n_sets)
          for (i in 1:n_sets) {
               tom_scaling_samples[[i]] <- as.dist(tom_list[[i]])[scale_sample]
               # Calculate 95th percentile
               scale_quant[i] <- quantile(tom_scaling_samples[[i]],
                                          probs = scale_percentile,
                                          type = 8)
               # Scale the male TOM
               if (i > 1) {
                  scale_powers[i] <- log(scale_quant[1]) / log(scale_quant[i])
                  tom_list[[i]] <- tom_list[[i]] ^ scale_powers[i]
               }
          }
          # Perhaps generate a quantile plot to check changes in TOM

          # Calculate consensus TOM
          tom <- pmin(tom_list[[1]],
                      tom_list[[2]])
     } else {
          data <- metabolites %>%
                  select(-RID)
          adjacency <- adjacency(data,
                                 power = soft_power)
          tom <- TOMsimilarity(adjacency)
     }
     diss_tom <- 1 - tom
     # Call the hierarchical clustering function
     prelim_tree <- hclust(as.dist(diss_tom),
                           method = "average")
     # Module identification using dynamic tree cut:
     prelim_modules <- cutreeDynamic(dendro = prelim_tree,
                                     distM = diss_tom,
                                     minClusterSize = min_module_size)
     # Convert numeric lables into colors
     prelim_colors <- labels2colors(prelim_modules)

     ## Merging modules that are very similar
     # The height cut of 0.10 means modules with correlations greater than 0.90
     me_diss_thres <- 0.10
     if (class(metabolites) == "list") {
          # Calculate eigengenes
          module_mes <- multiSetMEs(multi_expr,
                                    colors = NULL,
                                    universalColors = prelim_colors)

          # Consensus dissimilarity
          module_diss <- consensusMEDissimilarity(module_mes)
          # Call an automatic merging function
          merge <- mergeCloseModules(multi_expr,
                                     prelim_colors,
                                     cutHeight = me_diss_thres,
                                     verbose = 3)
     } else {
          mes <- moduleEigengenes(data,
                                  colors = prelim_colors)
          # Calculate dissimilarity of module eigengenes
          module_diss <- 1 - cor(mes$eigengenes)
          # Call an automatic merging function
          merge <- mergeCloseModules(data,
                                     prelim_colors,
                                     cutHeight = me_diss_thres,
                                     verbose = 3)
     }
     # Cluster modules tree
     me_tree <- hclust(as.dist(module_diss),
                       method = "average")

     # Keep merged colors for plot
     merged_colors <- merge$colors
     
     #Change names for kmeans
     print('-----Running K-means-----')
     names(merge)[1] <- "moduleColors"
     names(merge)[6] <- "MEs"
     #Apply KMeans
     k_means <- applyKMeans('metabolites',
                            merge,
                            data,
                            n.iterations = 100,
                            debug = F,
                            n.debug = 500,
                            net.type = "unsigned",
                            min.exchanged.genes = 20,
                            excludeGrey = F,
                            silent = F)
     #Compute new MEs
     n_partitions <- length(k_means$partitions)
     final_colors <- k_means$partitions[[n_partitions]]
     final_mes <- moduleEigengenes(data,
                                   colors = final_colors)
     final_mes <- final_mes$eigengenes

     #Compute module membership
     module_membership <- as.data.frame(cor(data,
                                            final_mes,
                                            use = "p"))

     #### PLOTS
     if (!is.null(plotname)) {
          # Plot dendrogram with initial and final modules
          filename <- paste0("../results/plots/",
                             plotname,
                             ".pdf")
          pdf(file = filename)
          plotDendroAndColors(prelim_tree,
                              cbind(prelim_colors,
                                    merged_colors,
                                    final_colors),
                              c("Dynamic Tree Cut",
                                "Merged dynamic",
                                "K-means"),
                              dendroLabels = FALSE,
                              hang = 0.03,
                              addGuide = TRUE,
                              guideHang = 0.05)
          dev.off()
     }

     if (!is.null(plot_modules)) {
          filename <- paste0("../results/plots/",
                             plot_modules,
                             ".pdf")
          pdf(file = filename)
          plot(me_tree,
               main = "Clustering of module eigenmetabolites",
               xlab = "",
               sub = "")
          # Plot the cut line into the dendrogram
          abline(h = me_diss_thres,
                 col = "red")
          dev.off()
     }

     if (class(metabolites) == "list") {
        rid <- metabolites[[1]]$RID
     } else {
        rid <- metabolites$RID
     }

     wgcna <- list("TOM" = tom,
                   "tree" = prelim_tree,
                   "colors" = final_colors,
                   "eigenmetabolites" = final_mes,
                   "RID" = rid,
                   "MM" = module_membership)
     return(wgcna)
}

plot_heatmap <- function(tom,
                         tree,
                         colors,
                         plotname="heatmap_p180") {
     # Plot the TOM as a heatmap

     plot_diss <- (tom)
     diag(plot_diss) <- NA
     filename <- paste0("../results/plots/",
                        plotname,
                        ".pdf")
     pdf(file = filename)
     TOMplot(plot_diss,
             tree,
             colors,
             main = "Network heatmap plot")
     dev.off()
}

save_wgcna <- function(wgcna,
                       suffix="p180") {
     # Save eigenmetabolites with RID and module colors to
     # external file
     #
     # Parameters
     # ----------
     # wgcna: list from compute_wgcna
     # suffix: str to add to filenames
     #

     me_filename <- paste0("../results/eigenmetabolites_",
                           suffix,
                           ".csv")

     if (class(wgcna$eigenmetabolites) == "data.frame") {
          eigenmetabolites <- wgcna$eigenmetabolites
     } else if (class(wgcna$eigenmetabolites) == "list") {
          eigenmetabolites <- wgcna$eigenmetabolites[[1]]$data
     }

     data <- cbind(wgcna$RID,
                   eigenmetabolites)
     colnames(data)[1] <- "RID"
     write.csv(data,
               file = me_filename,
               row.names = FALSE)

     colors_name <- paste0("../results/module_colors_",
                           suffix,
                           ".csv")
     color_table <- data.frame(cbind(wgcna$colors,
                                     rownames(wgcna$MM)))
     write.table(color_table,
                 file = colors_name,
                 row.names = FALSE,
                 col.names = FALSE,
                 quote = FALSE,
                 sep = ",")

     mm_filename <- paste0("../results/MM_",
                           suffix,
                           ".csv")
     write.csv(wgcna$MM,
               file = mm_filename)
}

export_to_cytoscape <- function(wgcna,
                                module,
                                node_names) {
     # Export data for cytoscape visualization
     #
     # Parameters
     # ----------
     # wgcna: wgcna result from compute_wgcna
     # module: name of module to export
     # node_names: name of nodes
     #

     in_module <- wgcna[[3]] == module
     edge_filename <- paste0("../results/CytoscapeInput-edges-",
                              module,
                              ".txt")
     node_filename <- paste0("../results/CytoscapeInput-nodes-",
                             module,
                             ".txt")

     cyt <- exportNetworkToCytoscape(wgcna[[1]][in_module, in_module],
                                     edgeFile = edge_filename,
                                     nodeFile = node_filename,
                                     weighted = TRUE,
                                     threshold = 0.05,
                                     nodeNames = node_names,
                                     nodeAttr = wgcna[[3]][in_module])
}
