library(shiny)
library(ROCR)
library(ggvis)
library(d3heatmap)
library(reshape2)
library(DESeq2)
library(edgeR)
library(phyloseq)
library(stats)
library(PathoStat)
library(plotly)
library(webshot)
library(vegan)
library(dplyr)
library(ape)






# Converts decimal percentage to string with specified digits
pct2str <- function(v, digits=2) {sprintf(paste0('%.',digits,'f'), v*100)}

# get RA from counts
getRelativeAbundance <- function(df){
  ra.out <- apply(df, 2, function(x) round(x/sum(x), digits = 4))
  return(ra.out)
}

getLogCPM <- function(df){
  logCPM.out <- apply(df, 2, function(y) log10(y*1e6/sum(y) + 1))
  return(logCPM.out)
}

vals <- reactiveValues(
    shiny.input = getShinyOption("pathostat.shinyInput"),
    shiny.input.backup = getShinyOption("pathostat.shinyInput"),
    taxdata = NULL,
    taxcountdata = NULL
)

updateTaxLevel <- function(){
  shinyInput <- vals$shiny.input
  pstat <- shinyInput$pstat

  updateSelectInput(session, "taxl",
                    choices = colnames(pstat@tax_table@.Data))
  updateSelectInput(session, "taxl.alpha",
                    choices = colnames(pstat@tax_table@.Data))
  updateSelectInput(session, "taxl.beta",
                    choices = colnames(pstat@tax_table@.Data))
  updateSelectInput(session, "taxl.pca",
                    choices = colnames(pstat@tax_table@.Data))
  updateSelectInput(session, "taxl.da",
                    choices = colnames(pstat@tax_table@.Data))
  updateSelectInput(session, "taxl.edger",
                    choices = colnames(pstat@tax_table@.Data))
  updateSelectInput(session, "taxl.pa",
                    choices = colnames(pstat@tax_table@.Data))
  updateSelectInput(session, "sra_taxlev",
                    choices = colnames(pstat@tax_table@.Data))
  updateSelectInput(session, "hmra_taxlev",
                    choices = colnames(pstat@tax_table@.Data))
  updateSelectInput(session, "taxl_single_species",
                    choices = colnames(pstat@tax_table@.Data))
}

# update samples
updateSample <- function(){
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    updateSelectInput(session, "filterSample",
                      choices = colnames(pstat@otu_table@.Data))
}



#Update covariate names
updateCovariate <- function(){
    shinyInput <- vals$shiny.input
    pstat <- shinyInput$pstat
    covariates <- colnames(sample_data(pstat))
    # choose the covariates that has less than 8 levels
    covariates.colorbar <- c()
    for (i in 1:length(covariates)){
        num.levels <- length(unique(sample_data(pstat)[[covariates[i]]]))
        if (num.levels < 8){
            covariates.colorbar <- c(covariates.colorbar, covariates[i])
        }
    }
    # choose the covariates that has 2 levels
    covariates.two.levels <- c()
    for (i in 1:length(covariates)){
        num.levels <- length(unique(sample_data(pstat)[[covariates[i]]]))
        if (num.levels == 2){
            covariates.two.levels <- c(covariates.two.levels, covariates[i])
        }
    }

    updateSelectInput(session, "select_covariate_condition_biomarker",
                      choices = covariates)
    updateSelectInput(session, "select_single_species_condition",
                      choices = covariates.colorbar)
    updateSelectInput(session, "select_target_condition_biomarker",
                      choices = covariates.colorbar)
    updateSelectInput(session, "select_condition_sample_filter",
                      choices = c("Read Number", covariates))
    updateSelectInput(session, "select_condition_sample_filter_micro",
                      choices = c("Taxon elements number", covariates))
    updateSelectInput(session, "select_condition_sample_filter_sidebar",
                      choices = c("Read Number", covariates))
    updateSelectInput(session, "select_condition_sample_distribution",
                      choices = covariates)
    updateSelectInput(session, "select_condition",
                      choices = covariates)
    updateSelectInput(session, "select_heatmap_condition_1",
                      choices = covariates.colorbar)
    updateSelectInput(session, "select_heatmap_condition_2",
                      choices = covariates.colorbar)
    updateSelectInput(session, "select_alpha_div_condition",
                      choices = covariates.colorbar)
    updateSelectInput(session, "select_beta_condition",
                      choices = covariates.two.levels)
    updateSelectInput(session, "select_beta_heatmap_condition_1",
                      choices = covariates.colorbar)
    updateSelectInput(session, "select_beta_heatmap_condition_2",
                      choices = covariates.colorbar)
    updateSelectInput(session, "select_pca_color",
                      choices = covariates)
    updateSelectInput(session, "select_pca_shape",
                      choices = covariates.colorbar)
    updateSelectInput(session, "da.condition",
                      choices = covariates)
    updateSelectInput(session, "edger.condition",
                      choices = covariates.colorbar)
    updateSelectInput(session, "da.condition.covariate",
                      choices = covariates)
    updateSelectInput(session, "pa.condition",
                      choices = covariates.colorbar)
    updateSelectInput(session, "sra_select_conditions",
                      choices = covariates)
    updateSelectInput(session, "gra_select_conditions",
                      choices = covariates)
    updateSelectInput(session, "hmra_select_conditions",
                      choices = covariates)
}
