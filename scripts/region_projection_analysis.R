library(tidyverse)
library(Seurat)
library(qs)
library(argparse)
library(rstan)
library(fitdistrplus)


utils <- new.env()
source("scripts/utils/utils.R", local = utils)

plot_utils <- new.env()
source("scripts/utils/plot_utils.R", local = plot_utils)



calculate_fisher_region_cell_property_enrichment <- function(regions,cell_properties){

    region_unique <- unique(na.omit(regions))
    cell_property_unique <- unique(na.omit(cell_properties))

    enrichment_matrix <- matrix(0, nrow = length(region_unique), ncol = length(cell_property_unique))

    rownames(enrichment_matrix) <- sort(region_unique)
    colnames(enrichment_matrix) <- sort(cell_property_unique)

    for (region in region_unique){
        for (cell_property in cell_property_unique){
            regions_test <- case_when(regions == region ~ region, TRUE ~ 'Others')
            cell_properties_test <- case_when(cell_properties == cell_property ~ cell_property, TRUE ~ 'Others')
            
            enrichment_matrix[region,cell_property] <- table(regions_test,cell_properties_test) %>% 
                as.data.frame() %>% 
                arrange(desc(regions_test != "Others"),desc(cell_properties_test != "Others")) %>% 
                pivot_wider(names_from = cell_properties_test, values_from = Freq) %>% 
                dplyr::select(-regions_test) %>% 
                fisher.test(alternative="greater") %>% 
                .[["p.value"]]
            
        }
    }

    return(enrichment_matrix)
}


plot_region_cell_property_enrichment <- function(enrichment_matrix, plot_folder){
    
    max_value <- -log(min(enrichment_matrix[enrichment_matrix != 0]))*2
    
    #Take inverse and log of enrichment matrix
    enrichment_matrix_inv <- -log(enrichment_matrix) %>% ifelse(is.infinite(.),max_value,.) %>% as.data.frame()

    #Pivot longe r to alow plotting
    enrichment_matrix_inv$regions <- rownames(enrichment_matrix_inv)
    enrichment_matrix_inv_long <- enrichment_matrix_inv %>% pivot_longer(-regions,names_to = 'projection', values_to = 'inv_p_value')

    enrichment_matrix_inv_long$inv_p_value_plot <- ""
    to_add <- enrichment_matrix_inv_long$inv_p_value > -log(0.05/length(enrichment_matrix_inv_long$inv_p_value))
    enrichment_matrix_inv_long$inv_p_value_plot[to_add] <- round(enrichment_matrix_inv_long$inv_p_value[to_add])

    enrichment_matrix_inv_long$projection <- factor(
        enrichment_matrix_inv_long$projection, 
        levels = c("BST","PAG", "BST+PAG", "LHA", "PBN", "LHA+PBN", "PVT", "CEA", "PVT+CEA", "Other")
    )

    
    # TODO: Køre med faerre region labels naar jeg faar dem af bernd. Eller lave et plot af hver maaske. 
    #Kan bare køre functioner to gange. en for hver region opløsning. Klart det paeneste og bedste
    region_cell_property_enrichment_plot <- enrichment_matrix_inv_long %>%
        filter(!projection == "Other") %>%
        ggplot() + 
        #geom_col(aes(x = projection, y = inv_p_value, fill = regions), position = "dodge") +
        geom_col(aes(x = regions, y = inv_p_value)) + facet_grid(~projection) + coord_flip() + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        geom_text(aes(x = regions, y = inv_p_value + 200, label = inv_p_value_plot))


    ggsave(paste0(plot_folder,"/enrichment_cell_by_region_corrected.png"),region_cell_property_enrichment_plot, width = 4000, height = 3000, units = "px")

}


#Following results are the main results for the region-projection analysis. Should be checked using know literature
plot_region_projection_absolute <- function(meta_data, plot_folder){

    plotting_colors <- meta_data$projection %>% unique() %>% plot_utils$generate_label_colors(.)

    region_projection_absolute_plot <- meta_data %>% 
        filter(projection != "Other") %>%
        ggplot() + geom_bar(aes(x = region_corrected, fill = projection), position = 'dodge') + 
        scale_fill_manual(values = plotting_colors) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    ggsave(paste0(plot_folder,"/region_projection_absolute_plot.png"),region_projection_absolute_plot, width = 4000, height = 3000, units = "px")
}

plot_region_projection_relative <- function(meta_data, plot_folder){

    plotting_colors <- meta_data$projection %>% unique() %>% plot_utils$generate_label_colors(.)

    region_projection_relative_plot <- meta_data %>% 
        group_by(run_index,region_corrected, projection) %>% 
        summarize(count = n()) %>% 
        mutate(total_count = sum(count),relative_count = count / total_count) %>% 
        filter(projection != "Other") %>%
        ggplot + geom_col(aes(x =  region_corrected, y = relative_count, fill = projection), position = "dodge") + 
        scale_fill_manual(values = plotting_colors) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    ggsave(paste0(plot_folder,"/region_projection_relative_plot.png"),region_projection_relative_plot, width = 4000, height = 3000, units = "px")
}


calculate_relative_counts <- function(meta_data){

    relative_count_df <- meta_data %>% 
        group_by(run_index,section_index, region_reduced, projection) %>% 
        summarize(count = n()) %>% 
        mutate(total_count = sum(count), relative_count = count / total_count) %>% 
        filter(projection != "Other", region_reduced != "Other")

    return(relative_count_df)
}

prepare_data_for_bayesian_model <- function(relative_count_df, region, projection_){

    relative_count_df_filtered <- relative_count_df %>% 
        filter(region_reduced == region, projection == projection_) 

    M <- relative_count_df_filtered %>% nrow()
    theta <- relative_count_df_filtered$relative_count
    n <- relative_count_df_filtered$total_count
    y <- relative_count_df_filtered$count

    return(list(
        M = M, 
        theta = theta,
        n = n,
        y = y
    ))
}

initialise_stan_model <- function(){

    # Stan model code (as a string)
    stan_model_code <- "
    data {
        int<lower=0> M;  // number of samples
        vector<lower=0,upper=1>[M] theta;  // parameters from Beta distribution (samples)
        int<lower=0> n[M];  // number of trials for each sample
        int<lower=0> y[M];  // observed successes for each sample
    }
    parameters {
        real<lower=0> alpha;  // shape parameter of prior for Beta distribution
        real<lower=0> beta;   // rate parameter of prior for Beta distribution
    }
    model {
        // Priors for alpha and beta
        alpha ~ gamma(0.001, 0.001);
        beta ~ gamma(0.001, 0.001);
    
        // Likelihood
        for (i in 1:M) {
            // Prior for each theta_i from Beta distribution
            theta[i] ~ beta(alpha, beta);
            
            // Likelihood for each sample
            y[i] ~ binomial(n[i], theta[i]);
        }
    }
    "

    # Compile the model
    stan_model_ <- rstan::stan_model(model_code = stan_model_code)

    return(stan_model_)
}

fit_bayesian_model <- function(spatial_object, data_to_fit){

    # Sample from the posterior
    fit <- sampling(stan_model, data = data_to_fit, iter = 2000, warmup = 500, chains = 4)

    return(fit)
}

marginalise_beta_dist <- function(fit,num_draws_per_sample){

    # Assuming `fit` is your Stan fit object
    posterior_samples <- rstan::extract(fit)
    alpha_samples <- posterior_samples$alpha
    beta_samples <- posterior_samples$beta

    # Number of posterior samples
    num_posterior_samples <- length(alpha_samples)

    # Initialize an empty vector to store the samples
    beta_draws <- numeric(num_posterior_samples * num_draws_per_sample)

    # Draw samples
    for (i in 1:num_posterior_samples) {
    beta_draws[((i - 1) * num_draws_per_sample + 1):(i * num_draws_per_sample)] <-
        rbeta(num_draws_per_sample, alpha_samples[i], beta_samples[i])
    }

    return(beta_draws)
}


fit_beta_MLE <- function(data_to_fit){

    samples <- data_to_fit$theta
    # Fit a beta distribution to the sample data using MLE
    fit <- fitdist(samples, "beta")

}

plot_region_projection_dist <- function(beta_dists, relative_count_df, plot_path){
    #beta_dists_sub <- beta_dists %>% 
    #    slice_sample(n = 1000000)
    
    # And point out the n highest
    quants <- beta_dists %>% 
        group_by(region, projection) %>% 
        summarize(q2.5 = quantile(samples, c(0.025)))
    
    beta_dist_region_projection_plot <- ggplot() + 
        geom_violin(data = beta_dists, mapping = aes(x = region, y = samples, fill = projection)) + facet_grid(~projection) + 
        geom_point(data = quants, mapping = aes(x = region, y = q2.5), size = 10, col = "red", shape = "-") + 
        geom_point(data = relative_count_df, mapping = aes(x = region_reduced, y = relative_count)) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        ylim(0,0.2)
    
    ggsave(plot_path, beta_dist_region_projection_plot, width = 4000, height = 3000, units = "px")
}

create_parser <- function(){
    # Input: NA
    #
    # Creates and outputs parser with arguments that allows interaction via CLI
    #
    # Output: Parser

    # Define parser
    parser <- ArgumentParser(description = "Script to generate spatial object from ROIs and transcripts data file")

    # Add argument to parser
    parser$add_argument("-o","--output-folder-path", type = "character", help = "Relative path to folder to save output")
    parser$add_argument("-p","--plot-folder-path", type = "character", help = "Relative path to folder to save plots")
    parser$add_argument("-sop","--spatial-object-path", type = "character", help="Relative path to spatial object in seurat format")

    return(parser)
}


main <- function(){
    # Input: NA
    #
    # Main function that defines parser argument, loads transcripts and ROIs, creates count matrix, creates spatial seurat object, saves the spatial seurat object
    #
    # Output: NA

    # Create parser and parse arguments
    parser <- create_parser()
    args <- parser$parse_args()

    # Retrieve arguments into corresponding variables
    path_to_output_folder <- args$output_folder_path
    path_to_plot_folder <- args$plot_folder_path
    spatial_object_path <- args$spatial_object_path
    
    dir.create(paste0(path_to_output_folder),showWarnings = FALSE)
    dir.create(paste0(path_to_plot_folder),showWarnings = FALSE)
    
    cat(paste0(utils$current_time()," Saving all output files to ",path_to_output_folder,"\n"))
    cat(paste0(utils$current_time()," Saving all plots to ",path_to_plot_folder,"\n"))

    cat(paste0(utils$current_time()," Loading spatial object from ",spatial_object_path,"\n"))
    spatial_object <- qread(spatial_object_path)

    # Summary plots

    cat(paste0(utils$current_time()," Generating region-projection plots with absolute values\n"))
    plot_region_projection_absolute(spatial_object@meta.data, path_to_plot_folder)

    cat(paste0(utils$current_time()," Generating region-projection plots with relative values\n"))
    plot_region_projection_relative(spatial_object@meta.data, path_to_plot_folder)

    # Fisher test

    cat(paste0(utils$current_time()," Performing enrichment test for region-projection pairs\n"))
    enrichment_matrix <- calculate_fisher_region_cell_property_enrichment(spatial_object$region_corrected, spatial_object$projection)

    cat(paste0(utils$current_time()," Generating plot from enrichment test for region-projection pairs\n"))
    plot_region_cell_property_enrichment(enrichment_matrix, path_to_plot_folder)

    # Bayesian model

    cat(paste0(utils$current_time()," Initialising stan model for bayesian model\n"))
    stan_model <- initialise_stan_model()
    
    relative_count_df <- calculate_relative_counts(spatial_object@meta.data)
    
    beta_dists_bayesian <- data.frame()
    beta_dists_MLE <- data.frame()
    
    regions <- relative_count_df$region_reduced %>% unique()
    projections <- relative_count_df$projection %>% unique()
    region_projection_pairs <- expand.grid(regions, projections) %>% 
        filter(Var1 != "Other")

    cat(paste0(utils$current_time()," Fitting data using bayesian model and MLE for all region_projection pairs:\n"))

    for (row_num in 1:nrow(region_projection_pairs)){

        region <- region_projection_pairs[row_num, "Var1"]
        projection <- region_projection_pairs[row_num, "Var2"]

        data_to_fit <- prepare_data_for_bayesian_model(relative_count_df,region,projection)

        if (data_to_fit$M > 1){
            fit_bayesian <- sampling(stan_model, data = data_to_fit, iter = 10000, warmup = 500, chains = 4, verbose = FALSE, refresh = 0)
            beta_dist_samples_bayesian <- marginalise_beta_dist(fit_bayesian, 100)
            
            beta_dist_bayesian <- data.frame(samples = beta_dist_samples_bayesian, region = region, projection = projection)
            beta_dists_bayesian <- rbind(beta_dists_bayesian, beta_dist_bayesian)

            fit_MLE <- fitdist(data_to_fit$theta, "beta")
            beta_dist_samples_MLE <- rbeta(1000, fit_MLE$estimate[1], fit_MLE$estimate[2])

            beta_dist_MLE <- data.frame(samples = beta_dist_samples_MLE, region = region, projection = projection)

            beta_dists_MLE <- rbind(beta_dists_MLE, beta_dist_MLE)
        }

        utils$loop_progress(row_num, nrow(region_projection_pairs))

    }

    cat(paste0(utils$current_time()," Generating region-projection plot on bayesian model fit\n"))
    plot_region_projection_dist(beta_dists_bayesian, relative_count_df, paste0(path_to_plot_folder,"/projection_region_plot_bayesian"))

    cat(paste0(utils$current_time()," Generating region-projection plot on MLE fit\n"))
    plot_region_projection_dist(beta_dists_MLE, relative_count_df, paste0(path_to_plot_folder,"/region_projection_plot_MLE.png"))

}

main()
