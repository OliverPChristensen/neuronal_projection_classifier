current_time <- function(){
    # Input: NA
    #
    # Extract and format current time
    #
    # Output: Current time

    # Extract and format current time
    return(paste0("[",format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),"]"))
}

loop_progress <- function(current_iteration, total_number_iteration){
    # Calculate the iterations that correspond to each 10% mark
    percent_marks <- ceiling(seq(0.1, 1.0, by = 0.1) * total_number_iteration)

    # Check progress and print every 10%
    if (current_iteration %in% percent_marks) {
        percent_completed <- sum(current_iteration >= percent_marks)*10
        return(cat(paste0(utils$current_time(),"  - Completed ",percent_completed, "%\n")))
    }
}