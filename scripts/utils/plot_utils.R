library(tidyverse)
library(scales)

generate_label_colors <- function(labels, grey_label = NULL){
    labels_unique <- labels %>% unique()
    
    # Generate colors
    plotting_colors <- labels_unique %>% length() %>% hue_pal()(.)
    names(plotting_colors) <- labels_unique

    if (!is.null(grey_label)){
        plotting_colors[grey_label] <- "#7F7F7F"
    }

    plotting_colors["Other"] <- "#7F7F7F"
    plotting_colors["Others"] <- "#7F7F7F"

    return(plotting_colors)
}

