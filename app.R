library(shiny)
library(SimpleMating)
library(dplyr)
library(purrr)

source("R/ranking_functions.R")

required_objects <- c(
  "pheno", "haplo.mat", "K", "crossPlan",
  "map", "marker_eff", "Markers"
)

ui <- fluidPage(
  
  titlePanel("SimpleMating App"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      fileInput("rds_file", "Upload SimpleMating RDS file", accept = ".rds"),
      
      uiOutput("analysis_ui"),
      uiOutput("trait_ui"),
      uiOutput("weights_ui"),
      uiOutput("propsel_ui"),
      uiOutput("method_ui"),
      
      actionButton("run_analysis", "Run Analysis"),
      
      downloadButton("download_results", "Download Analysis Results"),
      
      hr(),
      
      h4("Parent Ranking"),
      
      numericInput(
        "number_of_parents",
        "Number of Parents",
        value = 5,
        min = 3,
        max = 5,
        step = 1
      ),
      
      numericInput(
        "rel_cutoff",
        "Relationship Cutoff",
        value = 0,
        min = 0,
        max = 1,
        step = 0.01
      ),
      
      actionButton("run_ranking", "Rank Parent Sets"),
      
      downloadButton("download_ranking", "Download Ranked Parent Sets")
      
    ),
    
    mainPanel(
      
      h3("Data Check"),
      verbatimTextOutput("data_check"),
      
      h3("Selected Analysis"),
      textOutput("selected_analysis"),
      
      h3("Selected Traits"),
      textOutput("selected_traits"),
      
      h3("Weights Check"),
      textOutput("weights_check"),
      
      h3("Results"),
      tableOutput("results_table"),
      
      h3("Ranked Parent Sets"),
      tableOutput("ranking_table"),
      
      h3("Ranking Status"),
      textOutput("ranking_status"),
      
      h3("Analysis Status"),
      textOutput("analysis_status"),
      
    )
  )
)

server <- function(input, output, session) {
  
  uploaded_data <- reactive({
    req(input$rds_file)
    readRDS(input$rds_file$datapath)
  })
  
  data_valid <- reactive({
    data <- uploaded_data()
    all(required_objects %in% names(data))
  })
  
  trait_names <- reactive({
    data <- uploaded_data()
    setdiff(colnames(data$pheno), "Name")
  })
  
  output$data_check <- renderPrint({
    data <- uploaded_data()
    found <- names(data)
    missing <- setdiff(required_objects, found)
    
    cat("Objects found:\n")
    print(found)
    
    cat("\nMissing required objects:\n")
    if (length(missing) == 0) {
      cat("None - data looks valid.\n")
    } else {
      print(missing)
    }
    
    cat("\nAvailable traits:\n")
    print(trait_names())
  })
  
  output$analysis_ui <- renderUI({
    req(data_valid())
    
    tagList(
      selectInput(
        "analysis_type",
        "Analysis Type",
        choices = c(
          "MPV",
          "TGV",
          "Usefulness Additive",
          "Usefulness Additive + Dominance"
        )
      ),
      
      radioButtons(
        "trait_mode",
        "Trait Mode",
        choices = c("Single Trait", "Multi Trait")
      )
    )
  })
  
  output$trait_ui <- renderUI({
    req(data_valid())
    req(input$trait_mode)
    
    traits <- trait_names()
    
    if (input$trait_mode == "Single Trait") {
      selectInput("selected_traits", "Trait", choices = traits, selected = "SY")
    } else {
      selectInput(
        "selected_traits",
        "Traits",
        choices = traits,
        selected = c("DMD", "SY"),
        multiple = TRUE
      )
    }
  })
  
  output$weights_ui <- renderUI({
    req(input$trait_mode)
    
    if (input$trait_mode != "Multi Trait") {
      return(NULL)
    }
    
    req(input$selected_traits)
    
    selected_traits <- input$selected_traits
    default_weight <- round(1 / length(selected_traits), 3)
    
    tagList(
      h4("Trait Weights"),
      lapply(selected_traits, function(trait) {
        numericInput(
          inputId = paste0("weight_", trait),
          label = paste(trait, "weight"),
          value = default_weight,
          min = 0,
          max = 1,
          step = 0.01
        )
      })
    )
  })
  
  output$propsel_ui <- renderUI({
    
    req(input$analysis_type)
    
    if (input$analysis_type %in% c(
      "Usefulness Additive",
      "Usefulness Additive + Dominance"
    )) {
      
      numericInput(
        "propSel",
        "Proportion Selected",
        value = 0.05,
        min = 0.001,
        max = 0.999,
        step = 0.01
      )
      
    }
    
  })
  
  output$method_ui <- renderUI({
    
    req(input$analysis_type)
    
    if (input$analysis_type == "Usefulness Additive + Dominance") {
      
      selectInput(
        "method",
        "Method",
        choices = c("Phased", "NonPhased"),
        selected = "Phased"
      )
      
    }
    
  })
  
  selected_weights <- reactive({
    req(input$trait_mode)
    
    if (input$trait_mode == "Single Trait") {
      return(NULL)
    }
    
    req(input$selected_traits)
    
    weights <- sapply(input$selected_traits, function(trait) {
      input[[paste0("weight_", trait)]]
    })
    
    as.numeric(weights)
  })
  
  weights_valid <- reactive({
    if (input$trait_mode == "Single Trait") {
      return(TRUE)
    }
    
    weights <- selected_weights()
    
    if (any(is.na(weights))) {
      return(FALSE)
    }
    
    abs(sum(weights) - 1) < 0.0001
  })
  
  output$weights_check <- renderText({
    req(input$trait_mode)
    
    if (input$trait_mode == "Single Trait") {
      return("No weights needed for single-trait analysis.")
    }
    
    weights <- selected_weights()
    
    paste0(
      "Weight total: ",
      round(sum(weights), 4),
      ifelse(weights_valid(), " - OK", " - must equal 1")
    )
  })
  
  analysis_results <- eventReactive(input$run_analysis, {
    req(data_valid())
    req(input$analysis_type)
    req(input$trait_mode)
    req(input$selected_traits)
    
    validate(
      need(weights_valid(), "For multi-trait analysis, weights must add up to 1.")
    )
    
    data <- uploaded_data()
    
    if (input$analysis_type == "MPV") {
      
      if (input$trait_mode == "Single Trait") {
        
        Crit <- data.frame(
          Id = data$pheno[, "Name"],
          Criterion = data$pheno[, input$selected_traits]
        )
        
        result <- getMPV(
          MatePlan = data$crossPlan,
          Criterion = Crit,
          K = data$K
        )
        
      } else {
        
        Crit <- data.frame(
          Id = data$pheno[, "Name"],
          data$pheno[, input$selected_traits, drop = FALSE]
        )
        
        result <- getMPV(
          MatePlan = data$crossPlan,
          Criterion = Crit,
          K = data$K,
          Weights = selected_weights()
        )
      }
      
      return(result)
    }
    
    if (input$analysis_type == "TGV") {
      
      if (input$trait_mode == "Single Trait") {
        
        trait <- input$selected_traits
        
        result <- getTGV(
          MatePlan = data$crossPlan,
          Markers = data$Markers,
          addEff = data$marker_eff[, paste0(trait, "_add")],
          domEff = data$marker_eff[, paste0(trait, "_dom")],
          K = data$K
        )
        
      } else {
        
        traits <- input$selected_traits
        
        result <- getTGV(
          MatePlan = data$crossPlan,
          Markers = data$Markers,
          addEff = data$marker_eff[, paste0(traits, "_add"), drop = FALSE],
          domEff = data$marker_eff[, paste0(traits, "_dom"), drop = FALSE],
          K = data$K,
          Weights = selected_weights()
        )
      }
      
      return(result)
    }
    
    if (input$analysis_type == "Usefulness Additive") {
      
      Markers_02 <- data$Markers
      Markers_02[Markers_02 == 1] <- NA
      
      if (input$trait_mode == "Single Trait") {
        
        trait <- input$selected_traits
        
        result <- getUsefA(
          MatePlan = data$crossPlan,
          Markers = Markers_02,
          addEff = data$marker_eff[, paste0(trait, "_add")],
          Map.In = data$map,
          K = data$K,
          propSel = input$propSel
        )
        
      } else {
        
        traits <- input$selected_traits
        
        result <- getUsefA_mt(
          MatePlan = data$crossPlan,
          Markers = Markers_02,
          addEff = data$marker_eff[, paste0(traits, "_add"), drop = FALSE],
          Map.In = data$map,
          K = data$K,
          propSel = input$propSel,
          Weights = selected_weights()
        )
      }
      
      return(result[[2]])
    }
    
    if (input$analysis_type == "Usefulness Additive + Dominance") {
      
      req(input$method)
      
      if (input$method == "Phased") {
        marker_source <- data$haplo.mat
      } else {
        marker_source <- data$Markers
      }
      
      if (input$trait_mode == "Single Trait") {
        
        trait <- input$selected_traits
        
        result <- getUsefAD(
          MatePlan = data$crossPlan,
          Markers = marker_source,
          addEff = data$marker_eff[, paste0(trait, "_add")],
          domEff = data$marker_eff[, paste0(trait, "_dom")],
          Map.In = data$map,
          K = data$K,
          propSel = input$propSel,
          Method = input$method
        )
        
      } else {
        
        traits <- input$selected_traits
        
        result <- getUsefAD_mt(
          MatePlan = data$crossPlan,
          Markers = marker_source,
          addEff = data$marker_eff[, paste0(traits, "_add"), drop = FALSE],
          domEff = data$marker_eff[, paste0(traits, "_dom"), drop = FALSE],
          Map.In = data$map,
          K = data$K,
          propSel = input$propSel,
          Weights = selected_weights(),
          Method = input$method
        )
      }
      
      return(result[[2]])
    }
    
    data.frame(Message = "This analysis type is not connected yet.")
  })
  
  output$selected_analysis <- renderText({
    req(input$analysis_type)
    req(input$trait_mode)
    
    paste(
      "Selected analysis:",
      input$analysis_type,
      "| Trait mode:",
      input$trait_mode
    )
  })
  
  output$selected_traits <- renderText({
    req(input$selected_traits)
    paste("Selected traits:", paste(input$selected_traits, collapse = ", "))
  })
  
  output$analysis_status <- renderText({
    
    if (input$run_analysis == 0) {
      return("No analysis has been run yet.")
    }
    
    paste(input$analysis_type, "analysis completed successfully.")
  })
  
  ranking_results <- eventReactive(input$run_ranking, {
    
    validate(
      need(!is.null(analysis_results()), "Please run an analysis before ranking parent sets.")
    )
    
    results <- analysis_results()
    
    validate(
      need(all(c("Parent1", "Parent2", "Y", "K") %in% colnames(results)),
           "The analysis results are not suitable for parent ranking.")
    )
    
    rank_parent_sets(
      df = results,
      number_of_parents = input$number_of_parents,
      rel_cutoff = input$rel_cutoff
    )
    
  })
  
  output$download_results <- downloadHandler(
    
    filename = function() {
      paste0("analysis_results_", Sys.Date(), ".csv")
    },
    
    content = function(file) {
      write.csv(
        analysis_results(),
        file,
        row.names = FALSE
      )
    }
    
  )
  
  output$download_ranking <- downloadHandler(
    
    filename = function() {
      paste0("ranked_parent_sets_", Sys.Date(), ".csv")
    },
    
    content = function(file) {
      write.csv(
        ranking_results(),
        file,
        row.names = FALSE
      )
    }
    
  )
  
  output$ranking_status <- renderText({
    
    if (input$run_analysis == 0) {
      return("Please run an analysis before ranking parent sets.")
    }
    
    if (input$run_ranking == 0) {
      return("Ranking has not been run yet.")
    }
    
    "Ranking complete."
  })
  
  output$results_table <- renderTable({
    head(analysis_results(), 20)
  })
  
  output$ranking_table <- renderTable({
    head(ranking_results(), 20)
  })
}

shinyApp(ui, server)