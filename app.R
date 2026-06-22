library(shiny)
library(SimpleMating)
library(dplyr)
library(purrr)
library(DT)
library(openxlsx)

options(shiny.maxRequestSize = 100 * 1024^2)

source("R/ranking_functions.R")

required_objects <- c(
  "pheno", "haplo.mat", "K", "crossPlan",
  "map", "marker_eff", "Markers"
)

ui <- fluidPage(
  
  titlePanel(
    title = tags$div(
      "SimpleMating",
      tags$div(
        style = "font-size: 65%; font-weight: normal; color: #555555; margin-top: 3px; line-height: 1.2;",
        tags$b("Package:"), " M.A. Peixoto, R.R. Amadeu, L.L. Bhering, L.F.V. Ferrão, P.R. Muñoz, M.F.R. Resende Jr. | ",
        tags$b("Additional Parent Ranking Code:"), " Leif Skøt | ",
        tags$b("GUI:"), " Lasse Skøt"
      )
    )
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      
      h3("Step 1: Load Data"),
      
      radioButtons(
        "input_mode",
        "Input Mode",
        choices = c(
          "Upload prepared RDS file",
          "Upload individual data files"
        ),
        selected = "Upload prepared RDS file"
      ),
      
      uiOutput("data_upload_ui"),
      
      hr(),
      
      h3("Step 2: Configure Analysis"),
      uiOutput("analysis_ui"),
      uiOutput("trait_ui"),
      uiOutput("weights_ui"),
      uiOutput("propsel_ui"),
      uiOutput("method_ui"),
      
      hr(),
      
      uiOutput("run_analysis_ui"),
      
      uiOutput("parent_ranking_ui"),
      
      hr(),
      
      h3("Step 5: Export Results"),
      
      uiOutput("download_results_ui"),
      uiOutput("download_ranking_ui"),
      uiOutput("download_excel_ui")
      
    ),
    
    mainPanel(
      
      tabsetPanel(
        
        tabPanel(
          "Results",
          h3("Analysis Results"),
          DTOutput("results_table"),
          h3("Analysis Status"),
          textOutput("analysis_status")
        ),
        
        tabPanel(
          "Parent Ranking",
          h3("Ranked Parent Sets"),
          DTOutput("ranking_table"),
          h3("Ranking Status"),
          textOutput("ranking_status")
        ),
        
        tabPanel(
          "Data Check",
          h3("Data Check"),
          verbatimTextOutput("data_check")
        ),
        
        tabPanel(
          "Analysis Summary",
          h3("Analysis Summary"),
          verbatimTextOutput("analysis_summary")
        ),
        
        tabPanel(
          "Selections",
          h3("Selected Analysis"),
          textOutput("selected_analysis"),
          h3("Selected Traits"),
          textOutput("selected_traits"),
          h3("Weights Check"),
          textOutput("weights_check")
        )
        
    )
      
)
      
    )
  )


server <- function(input, output, session) {
  
  analysis_status_text <- reactiveVal("No analysis has been run yet.")
  
  ranking_status_text <- reactiveVal(
    "Please run an analysis before ranking parent sets."
  )
  
  uploaded_data <- reactive({
    
    req(input$input_mode)
    
    if (input$input_mode == "Upload prepared RDS file") {
      
      req(input$rds_file)
      readRDS(input$rds_file$datapath)
      
    } else {
      
      req(input$pheno_file)
      req(input$markers_file)
      req(input$map_file)
      req(input$k_file)
      req(input$marker_eff_file)
      req(input$haplo_file)
      req(input$crossplan_file)
      
      pheno <- read.table(
        input$pheno_file$datapath,
        header = TRUE
      )
      
      Markers <- read.table(
        input$markers_file$datapath,
        header = TRUE,
        row.names = 1
      )
      Markers <- as.matrix(Markers)
      
      map <- read.table(
        input$map_file$datapath,
        header = TRUE
      )
      
      K <- read.table(
        input$k_file$datapath,
        header = TRUE,
        row.names = 1
      )
      K <- as.matrix(K)
      
      marker_eff <- read.table(
        input$marker_eff_file$datapath,
        header = TRUE
      )
      
      haplo.mat <- read.table(
        input$haplo_file$datapath,
        header = TRUE,
        row.names = 1
      )
      haplo.mat <- as.matrix(haplo.mat)
      
      crossPlan <- read.table(
        input$crossplan_file$datapath,
        header = TRUE
      )
      
      list(
        pheno = pheno,
        haplo.mat = haplo.mat,
        K = K,
        crossPlan = crossPlan,
        map = map,
        marker_eff = marker_eff,
        Markers = Markers
      )
    }
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
    
    analysis_status_text("Running analysis...")
    ranking_status_text("Run parent ranking to view ranked parent sets.")
    
    withProgress(
      message = "Running analysis...",
      value = 0,
      {
    
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
        
        analysis_status_text(
          paste(input$analysis_type, "analysis completed successfully.")
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
        
        analysis_status_text(
          paste(input$analysis_type, "analysis completed successfully.")
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
        
        analysis_status_text(
          paste(input$analysis_type, "analysis completed successfully.")
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
        
        analysis_status_text(
          paste(input$analysis_type, "analysis completed successfully.")
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
        
        analysis_status_text(
          paste(input$analysis_type, "analysis completed successfully.")
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
        
        analysis_status_text(
          paste(input$analysis_type, "analysis completed successfully.")
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
        
        analysis_status_text(
          paste(input$analysis_type, "analysis completed successfully.")
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
        
        analysis_status_text(
          paste(input$analysis_type, "analysis completed successfully.")
        )
      }
      
      return(result[[2]])
    }
    
    data.frame(Message = "This analysis type is not connected yet.")
      }
    )
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
  
  output$analysis_summary <- renderPrint({
    
    req(input$run_analysis > 0)
    
    results <- analysis_results()
    
    best_row <- results[which.max(results$Y), ]
    worst_row <- results[which.min(results$Y), ]
    
    cat("ANALYSIS SETTINGS\n")
    cat("-----------------\n")
    cat("Analysis Type:", input$analysis_type, "\n")
    cat("Trait Mode:", input$trait_mode, "\n")
    cat("Selected Traits:", paste(input$selected_traits, collapse = ", "), "\n")
    
    if (input$trait_mode == "Multi Trait") {
      cat("Weights:", paste(selected_weights(), collapse = ", "), "\n")
    }
    
    if (input$analysis_type %in% c(
      "Usefulness Additive",
      "Usefulness Additive + Dominance"
    )) {
      cat("Proportion Selected:", input$propSel, "\n")
    }
    
    if (input$analysis_type == "Usefulness Additive + Dominance") {
      cat("Method:", input$method, "\n")
    }
    
    cat("\nRESULTS OVERVIEW\n")
    cat("----------------\n")
    cat("Number of crosses evaluated:", nrow(results), "\n")
    cat("Best Score:", round(max(results$Y, na.rm = TRUE), 5), "\n")
    cat("Average Score:", round(mean(results$Y, na.rm = TRUE), 5), "\n")
    cat("Worst Score:", round(min(results$Y, na.rm = TRUE), 5), "\n")
    
    cat("\nBEST CROSS\n")
    cat("----------\n")
    cat("Parent 1:", best_row$Parent1, "\n")
    cat("Parent 2:", best_row$Parent2, "\n")
    cat("Score:", round(best_row$Y, 5), "\n")
    cat("Relationship K:", round(best_row$K, 5), "\n")
    
    cat("\nWORST CROSS\n")
    cat("-----------\n")
    cat("Parent 1:", worst_row$Parent1, "\n")
    cat("Parent 2:", worst_row$Parent2, "\n")
    cat("Score:", round(worst_row$Y, 5), "\n")
    cat("Relationship K:", round(worst_row$K, 5), "\n")
    
    cat("\nRELATIONSHIP SUMMARY\n")
    cat("--------------------\n")
    cat("Average K:", round(mean(results$K, na.rm = TRUE), 5), "\n")
    cat("Minimum K:", round(min(results$K, na.rm = TRUE), 5), "\n")
    cat("Maximum K:", round(max(results$K, na.rm = TRUE), 5), "\n")
    
  })
  
  output$download_results_ui <- renderUI({
    
    if (is.null(input$run_analysis) || input$run_analysis == 0) {
      return(NULL)
    }
    
    downloadButton(
      "download_results",
      "Download Analysis Results"
    )
    
  })
  
  output$download_ranking_ui <- renderUI({
    
    if (is.null(input$run_ranking) || input$run_ranking == 0) {
      return(NULL)
    }
    
    downloadButton(
      "download_ranking",
      "Download Ranked Parent Sets"
    )
    
  })
  
  output$download_excel_ui <- renderUI({
    
    if (is.null(input$run_analysis) || input$run_analysis == 0) {
      return(NULL)
    }
    
    downloadButton(
      "download_excel",
      "Download Excel Workbook"
    )
    
  })
  
  output$run_analysis_ui <- renderUI({
    
    req(data_valid())
    
    tagList(
      h3("Step 3: Run Analysis"),
      actionButton(
        "run_analysis",
        "Run Analysis"
      )
    )
    
  })
  
  output$parent_ranking_ui <- renderUI({
    
    if (is.null(input$run_analysis) || input$run_analysis == 0) {
      return(NULL)
    }
    
    tagList(
      hr(),
      
      h3("Step 4: Rank Parent Sets"),
      
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
      
      actionButton(
        "run_ranking",
        "Rank Parent Sets"
      )
    )
  })
  
  output$data_upload_ui <- renderUI({
    
    req(input$input_mode)
    
    if (input$input_mode == "Upload prepared RDS file") {
      
      fileInput(
        "rds_file",
        "Upload SimpleMating RDS file",
        accept = ".rds"
      )
      
    } else {
      
      tagList(
        fileInput("pheno_file", "Phenotype file", accept = c(".txt", ".csv")),
        fileInput("markers_file", "Marker matrix file", accept = c(".txt", ".csv")),
        fileInput("map_file", "Genetic map file", accept = c(".txt", ".csv")),
        fileInput("k_file", "Relationship matrix file", accept = c(".txt", ".csv")),
        fileInput("marker_eff_file", "Marker effects file", accept = c(".txt", ".csv")),
        fileInput("haplo_file", "Haplotype matrix file", accept = c(".txt", ".csv")),
        fileInput("crossplan_file", "Crossplan file", accept = c(".txt", ".csv"))
      )
      
    }
  })
  
    output$analysis_status <- renderText({
      analysis_status_text()
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
      
      ranking_status_text("Ranking parent sets...")
      
      ranked <- withProgress(
        message = "Ranking parent sets...",
        value = 0,
        {
          incProgress(0.2, detail = "Checking parent combinations")
          
          ranked_result <- rank_parent_sets(
            df = results,
            number_of_parents = input$number_of_parents,
            rel_cutoff = input$rel_cutoff
          )
          
          incProgress(0.8, detail = "Preparing ranked table")
          
          ranked_result
        }
      )
      
      ranking_status_text("Ranking complete.")
      
      ranked
    })
  
  output$download_results <- downloadHandler(
    
    filename = function() {
      paste0("analysis_results_", Sys.Date(), ".csv")
    },
    
    content = function(file) {
      
      validate(
        need(input$run_analysis > 0,
             "Please run an analysis before downloading results.")
      )
      
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
      
      validate(
        need(input$run_ranking > 0,
             "Please run parent ranking before downloading ranking results.")
      )
      
      write.csv(
        ranking_results(),
        file,
        row.names = FALSE
      )
    }
    
  )
  
  output$download_excel <- downloadHandler(
    
    filename = function() {
      paste0("SimpleMating_results_", Sys.Date(), ".xlsx")
    },
    
    content = function(file) {
      
      wb <- createWorkbook()
      
      addWorksheet(wb, "Analysis_Results")
      writeData(wb, "Analysis_Results", analysis_results())
      
      if (input$run_ranking > 0) {
        addWorksheet(wb, "Parent_Rankings")
        writeData(wb, "Parent_Rankings", ranking_results())
      }
      
      results <- analysis_results()
      best_row <- results[which.max(results$Y), ]
      
      summary_df <- data.frame(
        Item = c(
          "Analysis Type",
          "Trait Mode",
          "Selected Traits",
          "Crosses Evaluated",
          "Best Parent 1",
          "Best Parent 2",
          "Best Score",
          "Average Score",
          "Worst Score"
        ),
        Value = c(
          input$analysis_type,
          input$trait_mode,
          paste(input$selected_traits, collapse = ", "),
          nrow(results),
          best_row$Parent1,
          best_row$Parent2,
          round(max(results$Y, na.rm = TRUE), 5),
          round(mean(results$Y, na.rm = TRUE), 5),
          round(min(results$Y, na.rm = TRUE), 5)
        )
      )
      
      addWorksheet(wb, "Analysis_Summary")
      writeData(wb, "Analysis_Summary", summary_df)
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  output$ranking_status <- renderText({
    
    if (is.null(input$run_analysis) || input$run_analysis == 0) {
      return("Please run an analysis before ranking parent sets.")
    }
    
    if (is.null(input$run_ranking) || input$run_ranking == 0) {
      return("Run parent ranking to view ranked parent sets.")
    }
    
    ranking_status_text()
  })
  
  output$results_table <- renderDT({
    datatable(
      analysis_results(),
      options = list(pageLength = 20, scrollX = TRUE)
    )
  })
  
  output$ranking_table <- renderDT({
    datatable(
      ranking_results(),
      options = list(pageLength = 20, scrollX = TRUE)
    )
  })
}

shinyApp(ui, server)