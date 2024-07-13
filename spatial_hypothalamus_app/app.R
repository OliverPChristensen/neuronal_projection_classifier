library(shiny)
library(tidyverse)
library(scales)
library(plotly)
library(sf)
library(qs)
#library(reactable)
library(cowplot)

shiny_data <- qread("./data/shiny_data.qs")

transcript_menu_items <- shiny_data$transcript_object$transcript %>% unique() %>% sort()


run_menu_items <- shiny_data$menu_converter$projection %>% unique()
bregma_levels <- shiny_data$menu_converter %>% 
    arrange(run_index) %>% 
    .[["bregma_level"]]



ui <- fluidPage(

    titlePanel("Spatial Hypothalamus Data"),
    sidebarLayout(
        sidebarPanel(

            radioButtons("run_menu_item", 
                label = "Projection:",
                choices = run_menu_items,
                selected = run_menu_items[1],
                inline = TRUE
            ),

            selectizeInput("section_menu_item", 
                label = "Bregma level:",
                choices = NULL, 
                selected = "",
                options = list(placeholder = 'Select a bregma level', onInitialize = I('function() { this.setValue(""); }'))
            ),

            selectizeInput("transcript_menu_item", 
                label = "Transcripts:",
                choices = transcript_menu_items,
                multiple = TRUE,
                options = list(maxItems = 9)
            ),

            checkboxInput("show_regions", 
                label = "Show anatomical regions", 
                value = TRUE
            ),
            
            checkboxInput("interactive", 
                label = "Use interactive plotting", 
                value = FALSE
            ),

            checkboxInput("downsample", 
                label = "Downsample highly abundant transcripts when using interactive plotting", 
                value = TRUE
            ),

            actionButton("render_plot", 
                label = "Plot", 
                disabled = TRUE
            ),

            HTML("<br><br>"),

            uiOutput("downsample_warning")


        ),

        mainPanel(
            tabsetPanel(
                tabPanel("Spatial plot",
                    splitLayout(
                        cellWidths = c("800px", "300px"),
                        cellArgs = list(style = "padding: 0px; margin: 0px 10px;"),
                        uiOutput("plot_ui"),
                        plotOutput("legend")
                    )

                ),

                tabPanel("Summary barplots",

                    tags$div(style = "margin-top: 50px;"),
                    
                    plotOutput("cell_region_barplot", width = "1100px", height = "800px"),

                    tags$div(style = "margin-top: 50px;"),
                    
                    #plotOutput("transcript_region_barplot", width = "1100px", height = "1500px")
                    uiOutput("dynamic_plot")
                )
            )
        )
    )
)

server <- function(input, output, session) {

    observeEvent(input$run_menu_item, {
        print(paste(Sys.time(),"update_section_menu"))
        updateSelectizeInput(session, "section_menu_item", choices = NULL, selected = "", server = TRUE)

        section_menu_items <- switch(input$run_menu_item,
            "BST and PAG" = sort(bregma_levels[1:8]),
            "LHA and PBN" = sort(bregma_levels[9:16]),
            "PVT and CEA" = sort(bregma_levels[17:24])
        )

        updateSelectizeInput(session, "section_menu_item", choices = section_menu_items, selected = "", server = TRUE)

    })

    observe({
        if (input$section_menu_item == "") {
            updateActionButton(session, "render_plot", label = "Plot", disabled = TRUE)
        } else {
            updateActionButton(session, "render_plot", label = "Plot", disabled = FALSE)
        }
    })

    zoom_config <- reactiveValues(
        previous_zoom_state = NULL,
        previous_section = list(
            run_menu_item = "NULL",
            section_menu_item = "NULL",
            interactive = "NULL"
        )
    )

    # This needs to be separat from plotly_plot because prev zoom state needs to be changed to NULL. not enough to just check when creating plotly plot
    observeEvent(input$render_plot, {
        print(paste(Sys.time(),"reset_zoom"))
        if (zoom_config$previous_section$run_menu_item != input$run_menu_item | zoom_config$previous_section$section_menu_item != input$section_menu_item | zoom_config$previous_section$interactive != input$interactive){
            print("zoom state null")
            zoom_config$previous_zoom_state <- NULL
        }
        
        zoom_config$previous_section <- list(
            run_menu_item = input$run_menu_item,
            section_menu_item = input$section_menu_item,
            interactive = input$interactive 
        )

    })

    #This has to be observeEvent because needs to countiuously moniter zoom, because zoom state is lost later
    observeEvent(event_data("plotly_relayout"),{
        print(paste(Sys.time(),"update_plot"))

        zoom_state <- event_data("plotly_relayout")

        if (all(c("xaxis.range[0]", "xaxis.range[1]", "yaxis.range[0]", "yaxis.range[1]") %in% names(zoom_state))){
            
            zoom_config$previous_zoom_state <- list("x1" = zoom_state$`xaxis.range[0]`, "x2" = zoom_state$`xaxis.range[1]`, "y1" = zoom_state$`yaxis.range[0]`, "y2" = zoom_state$`yaxis.range[1]`)

        }

    })

    filter_data <- eventReactive(input$render_plot, {
        print(paste(Sys.time(),"filter_data start"))

        selected_run_index <- shiny_data$menu_converter %>%
            filter(projection == input$run_menu_item, bregma_level == input$section_menu_item) %>% 
            .[["run_index"]]

        selected_section_index <- shiny_data$menu_converter %>%
            filter(projection == input$run_menu_item, bregma_level == input$section_menu_item) %>% 
            .[["section_index"]]
    


        cell_object <- shiny_data$cell_object %>% filter(run_index == selected_run_index, section_index == selected_section_index)

        if (input$interactive){
            concant <- c("Virus1","Virus2","Cell area",input$transcript_menu_item)
            #concant <- c("section_index",input$transcript_menu_item)

            cell_object_centroid <- shiny_data$cell_object_centroid %>% 
                filter(run_index == selected_run_index, section_index == selected_section_index)
            
            print(paste(Sys.time(),"tooltip start"))
            cell_object_centroid_tooltip <- cell_object_centroid %>% 
                st_drop_geometry() %>% 
                mutate(tooltip = pmap_chr(select(., all_of(concant)), ~ paste(concant, c(...), sep = ": ", collapse = "<br>")))

            cell_object_centroid$tooltip <- cell_object_centroid_tooltip$tooltip
            print(paste(Sys.time(),"tooltip start"))
        } else {
            cell_object_centroid = NULL
        }


        if (length(input$transcript_menu_item) > 0){
            if (input$interactive & input$downsample){
                transcript_object <- shiny_data$transcript_object_downsampled %>% filter(run_index == selected_run_index, section_index == selected_section_index, transcript %in% input$transcript_menu_item)
            } else {
                transcript_object <- shiny_data$transcript_object %>% filter(run_index == selected_run_index, section_index == selected_section_index, transcript %in% input$transcript_menu_item)
            }
            transcript_region_table <- shiny_data$transcript_region_table %>% filter(run_index == selected_run_index, section_index == selected_section_index, transcript %in% input$transcript_menu_item)
        } else {
            transcript_object <- NULL
            transcript_region_table <- NULL
        }

        if (length(input$show_regions) > 0){
            region_object <- shiny_data$region_object %>% filter(run_index == selected_run_index, section_index == selected_section_index)
        } else {
            region_object <- NULL
        }

        print(paste(Sys.time(),"filter_data end"))

        return(list(
            cell_object = cell_object, 
            cell_object_centroid = cell_object_centroid,
            transcript_object = transcript_object, 
            region_object = region_object, 
            transcript_region_table = transcript_region_table
        ))
    })



    create_ggplot <- eventReactive(input$render_plot, {
        filtered_data <- filter_data()
        print(paste(Sys.time(),"create_ggplot start"))

        projection_unique <- filtered_data$cell_object$projection %>% unique()
        
        # Generate colors
        plotting_colors <- projection_unique %>% length() %>% hue_pal()(.)
        names(plotting_colors) <- projection_unique
        
        plotting_colors["Other"] <- "#7F7F7F"
        
        
        plot <- ggplot() + 
            geom_sf(data = filtered_data$cell_object, aes(fill = projection), lwd = 0.1) + 
            scale_fill_manual(values = plotting_colors)
        
        if (!is.null(filtered_data$transcript_object)){
            plot <- plot + 
                geom_sf(data = filtered_data$transcript_object, aes(color = transcript), size = 0.5, shape = 15) + 
                scale_colour_brewer(palette = "Set1") + 
                guides(color = guide_legend(override.aes = list(size = 6)))
        }
        
        if (input$show_regions){
            plot <- plot + 
                geom_sf(data = filtered_data$region_object,fill = NA, color = "white", lwd = 0.1) + 
                geom_sf_text(data = filtered_data$region_object, aes(label = region), color = "white", size = 6, fontface = "bold")
        }
        
        if (!is.null(filtered_data$cell_object_centroid)){
            plot <- plot + 
                geom_sf(data = filtered_data$cell_object_centroid, aes(text = tooltip), size = 5, alpha = 0)
        }
        
        legend_plot <- plot + 
            labs(fill = "Projection", color = "Transcript") + 
            theme(text = element_text(color = "black"),
                legend.margin = margin(0, 0, 0, 0),
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 25,face = "bold"))
        
        legend <- get_legend(legend_plot)

        plot_plot <- plot + 
            theme_void() + 
            theme(text = element_text(color = "black"),
                legend.position = "none",
                panel.background = element_rect(fill = "black", color = NA),
                panel.border = element_rect(color = "#595959", fill = NA, linewidth = 2)
            )

        print(paste(Sys.time(),"create_ggplot end"))

        return(list(
            plot = plot_plot,
            legend = legend
        ))
    })

    create_plotly_plot <- eventReactive(input$render_plot, {
        plot_and_legend <- create_ggplot()

        print(paste(Sys.time(),"create_plotly_plot start"))

        if (!is.null(zoom_config$previous_zoom_state)){
            print("zoom_plot")

            plotly_plot <- ggplotly(plot_and_legend$plot, tooltip = "text") %>%
                layout(
                    xaxis = list(range = c(zoom_config$previous_zoom_state$x1, zoom_config$previous_zoom_state$x2)),
                    yaxis = list(range = c(zoom_config$previous_zoom_state$y1, zoom_config$previous_zoom_state$y2))
                )

        } else {
            print("non_zoom_plot")
            plotly_plot <- ggplotly(plot_and_legend$plot, tooltip = "text")
        
        }
        
        plotly_plot <- event_register(plotly_plot, 'plotly_relayout')
        print(paste(Sys.time(),"create_plotly_plot end"))

        return(plotly_plot)
    })

    plot_interactive <- eventReactive(input$render_plot,{
        print(paste(Sys.time(),"update_interactive"))
        return(input$interactive)
    })

    output$plot_ui <- renderUI({
        print(paste(Sys.time(),"render_plot start"))
        if (plot_interactive()) {
            output$spatial_plot <- NULL
            output$spatial_plotly <- renderPlotly({
                create_plotly_plot()
            })
            print(paste(Sys.time(),"render_plot end"))
            plotlyOutput("spatial_plotly", width = "800px", height = "800px")
        } else {
            output$spatial_plotly <- NULL
            output$spatial_plot <- renderPlot({
                create_ggplot()
            })
            print(paste(Sys.time(),"render_plot end"))
            plotOutput("spatial_plot", width = "800px", height = "800px")
        }
    })

    output$legend <- renderPlot({
        plot_and_legend <- create_ggplot()
        legend_plot <- ggplot() +
            theme_void() +
            annotation_custom(plot_and_legend$legend) + 
            theme(legend.position = "bottom")
        legend_plot
    })

    output$cell_region_barplot <- renderPlot({
        filtered_data <- filter_data()
        cell_object <- filtered_data$cell_object

        cell_region_barplot <- cell_object %>%
            group_by(region,projection) %>% 
            summarize(counts = n()) %>% 
            filter(!is.na(region)) %>% 
                ggplot(aes(x = projection, y = counts, fill = region)) +
                geom_col(position = position_dodge()) +
                geom_text(
                    aes(label = counts), 
                    position = position_dodge(width = 0.9), 
                    hjust = -0.5, 
                    size = 4, 
                    color = "black"
                ) + 
            labs(title = "Projection cell count across regions in chosen section") + 
            theme_classic() + 
            theme(
                text = element_text(size = 15, color = "black"),
                plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 17, face = "bold"),
                axis.title.y = element_text(size = 17, face = "bold"),
            ) + coord_flip()

        cell_region_barplot
    })

    plot_height <- eventReactive(input$render_plot,{
        return(paste0(200 + 200*length(input$transcript_menu_item),"px"))
    })

    output$dynamic_plot <- renderUI({
        plotOutput("transcript_region_barplot", width = "1100px", height = plot_height())
    })

    output$transcript_region_barplot <- renderPlot({
        filtered_data <- filter_data()
        
        if (!is.null(filtered_data$transcript_region_table)){
            transcript_region_barplot <- filtered_data$transcript_region_table %>%
                ggplot(aes(x = transcript, y = counts, fill = region)) +
                geom_col(position = position_dodge()) +
                geom_text(
                    aes(label = counts), 
                    position = position_dodge(width = 0.9), 
                    hjust = -0.5, 
                    size = 4, 
                    color = "black"
                ) + 
                labs(title = "Transcript count across regions in chosen section") + 
                theme_classic() + 
                theme(
                    text = element_text(size = 15, color = "black"),
                    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
                    axis.title.x = element_text(size = 17, face = "bold"),
                    axis.title.y = element_text(size = 17, face = "bold"),
                ) + coord_flip()

            transcript_region_barplot
        }
    })

    output$downsample_warning <- renderUI({
        filtered_data <- filter_data()
        if (plot_interactive() & input$downsample & !is.null(filtered_data$transcript_region_table)){
            downsampled_transcripts <- filtered_data$transcript_region_table %>% 
                group_by(transcript) %>%
                summarize(counts1 = sum(counts)) %>% 
                filter(counts1 > 10000) %>%
                .[["transcript"]]
            
            if (length(downsampled_transcripts) > 0){
                HTML(paste0('<p style="color: red; font-size: 15px; font-weight: bold;">',
                    paste0("Downsampled transcripts: ",paste(downsampled_transcripts, collapse = ", ")),
                    '</p>'
                ))
                
            }
        }
        
    })
}

shinyApp(ui = ui, server = server)