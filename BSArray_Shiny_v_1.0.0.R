options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB

library(DT)
library(shiny)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(purrr)
library(IRanges)
library(openxlsx)
library(tibble)

# Load BSAarray functions
source("BSArray_v1.0.0.R")

# Load SNP database once when app starts
snp_db <- read_csv("SNP_master_DB.csv", show_col_types = FALSE) %>%
  mutate(
    Chr_v6 = paste0("Gm", str_pad(as.integer(Chr), 2, pad = "0")),
    Position_v6 = as.numeric(Position_v6)
  )

gene_db <- read_csv("gene_summary_clean.csv", show_col_types = FALSE) %>%
  separate(Interval, into = c("Start", "End"), sep = "-", convert = TRUE) %>%
  mutate(
    Start = as.numeric(Start),
    End = as.numeric(End)
  )

# ---- Helper: summarize merged W-based putative regions ----
summarize_putative_regions <- function(results_obj) {
  regions_df <- results_obj$Shaded_Regions
  windows_df <- results_obj$Results
  markers_df <- results_obj$Significant_Markers
  
  empty_tbl <- tibble(
    Chromosome = character(),
    Interval = character(),
    Start = character(),
    End = character(),
    `Length (Mb)` = character(),
    `Number of Windows` = integer(),
    `Average SNPs/Mb` = numeric(),
    `Peak Marker` = character(),
    `Peak Marker Position` = character(),
    `Peak Marker Score` = character()
  )
  
  if (is.null(regions_df) || nrow(regions_df) == 0) {
    return(empty_tbl)
  }
  
  purrr::pmap_dfr(
    list(regions_df$Chr, regions_df$Region_Start, regions_df$Region_End),
    function(Chr, Region_Start, Region_End) {
      
      region_windows <- windows_df %>%
        dplyr::filter(
          .data$Chr == !!Chr,
          .data$Window_Start <= !!Region_End,
          .data$Window_End >= !!Region_Start
        )
      
      region_markers <- markers_df %>%
        dplyr::filter(
          .data$Chr == !!Chr,
          .data$Region_Start == !!Region_Start,
          .data$Region_End == !!Region_End
        ) %>%
        dplyr::mutate(
          Combined_Score = .data$Abs_Delta_AF * .data$W_SNP
        )
      
      length_mb <- (Region_End - Region_Start + 1) / 1e6
      n_windows <- nrow(region_windows)
      
      avg_snps_mb <- if (n_windows > 0) {
        round(
          mean(
            region_windows$SNP_Count /
              ((region_windows$Window_End - region_windows$Window_Start) / 1e6),
            na.rm = TRUE
          ),
          1
        )
      } else {
        NA_real_
      }
      
      if (nrow(region_markers) > 0) {
        peak_marker_row <- region_markers %>%
          dplyr::arrange(
            dplyr::desc(.data$Combined_Score),
            dplyr::desc(.data$Abs_Delta_AF),
            dplyr::desc(.data$W_SNP),
            .data$Position
          ) %>%
          dplyr::slice(1)
        
        peak_marker <- peak_marker_row$SNP_Name[[1]]
        peak_marker_pos <- format(peak_marker_row$Position[[1]], big.mark = ",")
        peak_combined <- sprintf("%.3f", peak_marker_row$Combined_Score[[1]])
      } else {
        peak_marker <- NA_character_
        peak_marker_pos <- NA_character_
        peak_combined <- NA_character_
      }
      
      start_mb <- round(Region_Start / 1e6, 2)
      end_mb   <- round(Region_End / 1e6, 2)
      interval <- paste0(sprintf("%.2f", start_mb), "-", sprintf("%.2f", end_mb))
      
      tibble(
        Chromosome = as.character(Chr),
        Interval = interval,
        Start = format(Region_Start, big.mark = ","),
        End = format(Region_End, big.mark = ","),
        `Length (Mb)` = sprintf("%.2f", length_mb),
        `Number of Windows` = n_windows,
        `Average SNPs/Mb` = avg_snps_mb,
        `Peak Marker` = peak_marker,
        `Peak Marker Position` = peak_marker_pos,
        `Peak Marker Score` = peak_combined
      )
    }
  )
}

# ---- Helper: clean candidate gene display text ----
clean_gene_display_text <- function(x) {
  x %>%
    stringr::str_replace_all("‚Äî", " — ") %>%
    stringr::str_replace_all("<br\\s*/?>", "; ") %>%
    stringr::str_replace_all("\\s*;\\s*", "; ") %>%
    stringr::str_replace_all(";{2,}", "; ") %>%
    stringr::str_squish()
}

# ---- Helper: strip HTML for export ----
strip_gene_html <- function(x) {
  x %>%
    clean_gene_display_text() %>%
    stringr::str_replace_all("<.*?>", "")
}

rename_classic_output_samples <- function(df, classic_results) {
  df %>%
    dplyr::rename(
      !!classic_results$Bulk_1_Name := R_Bulk,
      !!classic_results$Bulk_2_Name := S_Bulk,
      !!classic_results$Parent_1_Name := R_Parent,
      !!classic_results$Parent_2_Name := S_Parent
    )
}

# ---- Helper: summarize informative intervals (Classic mode) ----
summarize_informative_intervals <- function(classic_results, min_snp = 20, gap_bp = 2e6) {
  df <- classic_results$Informative_SNPs
  
  if (is.null(df) || nrow(df) == 0) {
    return(tibble(
      Chromosome = character(),
      Interval = character(),
      `Interval Size (Mb)` = numeric(),
      `Number of SNPs` = integer(),
      `First Marker` = character(),
      `Last Marker` = character()
    ))
  }
  
  df <- df %>%
    arrange(Chr, Position) %>%
    group_by(Chr) %>%
    mutate(
      gap = Position - lag(Position),
      new_group = ifelse(is.na(gap) | gap > gap_bp, 1, 0),
      group_id = cumsum(new_group)
    ) %>%
    ungroup()
  
  regions <- df %>%
    group_by(Chr, group_id) %>%
    summarise(
      Start = min(Position, na.rm = TRUE),
      End = max(Position, na.rm = TRUE),
      `Number of SNPs` = n(),
      `First Marker` = SNP_Name[1],
      `Last Marker` = SNP_Name[dplyr::n()],
      .groups = "drop"
    ) %>%
    filter(`Number of SNPs` >= min_snp) %>%
    mutate(
      Start_Mb = round(Start / 1e6, 2),
      End_Mb   = round(End / 1e6, 2),
      Interval = paste0(sprintf("%.2f", Start_Mb), "-", sprintf("%.2f", End_Mb)),
      `Interval Size (Mb)` = round((End - Start) / 1e6, 2)
    ) %>%
    select(
      Chromosome = Chr,
      Interval,
      `Interval Size (Mb)`,
      `Number of SNPs`,
      `First Marker`,
      `Last Marker`
    ) %>%
    arrange(Chromosome, Interval)
  
  regions
}

# ---- Helper: keep only informative SNPs within retained classic intervals ----
filter_classic_informative_to_intervals <- function(classic_results, min_snp = 20) {
  intervals <- summarize_informative_intervals(
    classic_results,
    min_snp = min_snp
  )
  
  if (is.null(intervals) || nrow(intervals) == 0) {
    return(tibble())
  }
  
  
  classic_results$Informative_SNPs %>%
    dplyr::inner_join(
      intervals %>%
        dplyr::mutate(
          Start_bp = as.numeric(sub("-.*", "", Interval)) * 1e6,
          End_bp   = as.numeric(sub(".*-", "", Interval)) * 1e6
        ) %>%
        dplyr::select(Chromosome, Start_bp, End_bp),
      by = c("Chr" = "Chromosome"),
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(Position >= Start_bp, Position <= End_bp)
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .sidebar-section {
        border: 1px solid #d9d9d9;
        border-radius: 8px;
        padding: 12px;
        margin-bottom: 15px;
        background-color: #fafafa;
      }
      .sidebar-section-strong {
        border: 2px solid #bfc7d5;
        border-radius: 8px;
        padding: 12px;
        margin-bottom: 15px;
        background-color: #f4f7fb;
      }
      .sidebar-section-title {
        font-weight: 700;
        margin-bottom: 10px;
        font-size: 15px;
      }
      .compact-controls .form-group {
        margin-bottom: 10px;
      }
      .compact-controls .shiny-input-container {
        margin-bottom: 8px;
      }
      .run-button-wrap {
        margin-top: 12px;
      }
      .small-note {
        font-size: 12px;
        color: #777;
        margin-top: 6px;
      }
      table.dataTable td {
        white-space: normal !important;
        vertical-align: top;
      }
    "))
  ),
  
  titlePanel("BSAarray: Automated SNP Array Pipeline for QTL Discovery"),
  p("A Shiny interface for BSAarray and classic parent–bulk matching workflows"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      div(
        class = "sidebar-section-strong",
        div(class = "sidebar-section-title", "Analysis Mode"),
        selectInput(
          "analysis_mode",
          "Choose mode",
          choices = c(
            "BSAarray" = "BSArray",
            "Classic" = "classic"
          ),
          selected = "BSArray"
        )
      ),
      
      conditionalPanel(
        condition = "input.analysis_mode == 'BSArray'",
        div(
          class = "sidebar-section",
          div(class = "sidebar-section-title", "Input Files"),
          fileInput("bulkfile", "Upload BSA file"),
          fileInput("parentfile", "Upload Parents file")
        ),
        div(
          class = "sidebar-section compact-controls",
          div(class = "sidebar-section-title", "Sample Mapping"),
          uiOutput("bulk_select_ui"),
          uiOutput("parent_map_ui"),
          div(
            class = "small-note",
            "Bulk and parent menus are populated automatically from the uploaded GenomeStudio files."
          )
        ),
        div(
          class = "sidebar-section-strong",
          div(class = "sidebar-section-title", "Analysis Parameters"),
          numericInput("window_size", "Window size", value = 1000000),
          numericInput("step_size", "Step size", value = 100000),
          numericInput("min_snp", "Minimum SNP per window", value = 20),
          numericInput(
            "flank_bp",
            "Gene annotation window (bp)",
            value = 10000,
            min = 0,
            step = 5000
          ),
          div(
            class = "run-button-wrap",
            actionButton("run_analysis", "Run Analysis", class = "btn-primary")
          )
        )
      ),
      conditionalPanel(
        condition = "input.analysis_mode == 'classic'",
        div(
          class = "sidebar-section",
          div(class = "sidebar-section-title", "Input Files"),
          fileInput("bulk_matrix_file", "Upload Bulk Genotype File"),
          fileInput("parent_matrix_file", "Upload Parent Genotype File")
        ),
        div(
          class = "sidebar-section compact-controls",
          div(class = "sidebar-section-title", "Sample Mapping"),
          uiOutput("matrix_bulk_select_ui"),
          uiOutput("matrix_parent_select_ui"),
          div(
            class = "small-note",
            "Sample menus are populated automatically from the uploaded GenomeStudio genotype files."
          )
        ),
        div(
          class = "sidebar-section-strong",
          div(class = "sidebar-section-title", "Run"),
          selectInput(
            "classic_matching_model",
            "Classic matching model",
            choices = c(
              "Strict parent-bulk matching" = "strict",
              "Dominant favorable-parent model" = "dominant_favorable",
              "Recessive favorable-parent model" = "recessive_favorable"
            ),
            selected = "strict"
          ),
          div(
            class = "small-note",
            "Use the dominant model when the favorable bulk may contain heterozygotes, such as F2 or F2:3 populations with dominant QTL."
          ),
          numericInput(
            "classic_min_snp",
            "Minimum SNPs per informative interval",
            value = 20,
            min = 1,
            step = 1
          ),
          div(
            class = "run-button-wrap",
            actionButton("run_classic", "Run Classic Mode", class = "btn-primary")
          )
        )
      )
    ),
    
    mainPanel(
      width = 9,
      
      conditionalPanel(
        condition = "input.analysis_mode == 'BSArray'",
        tabsetPanel(
          tabPanel(
            "Genome Overview",
            wellPanel(
              h4("Status"),
              textOutput("status")
            ),
            br(),
            h3("All Chromosomes"),
            downloadButton("download_all_chr_plot", "Download All Chromosomes Plot"),
            downloadButton("download_full_results", "Download Full Results (Excel)"),
            br(), br(),
            plotOutput("sig_plot", height = "800px")
          ),
          
          tabPanel(
            "Chromosome Detail",
            selectInput(
              "chr_scope",
              "Show chromosomes:",
              choices = c(
                "All chromosomes" = "all",
                "≥ Significant (99%)" = "high"
              ),
              selected = "high"
            ),
            uiOutput("chr_summary"),
            selectInput("chr_select", "Select Chromosome", choices = NULL),
            downloadButton("download_chr_plot", "Download Plot"),
            br(), br(),
            plotOutput("chr_plot", height = "600px")
          ),
          
          tabPanel(
            "Putative Regions",
            h4("Putative Regions"),
            p(
              "Putative genomic regions identified from the analysis. The peak marker is selected as the marker with the strongest combined signal, based on both allele-frequency difference and marker-level support.",
              style = "font-size: 14px; color: gray;"
            ),
            downloadButton("download_region_table", "Download Putative Regions"),
            br(), br(),
            DTOutput("region_table")
          ),
          
          tabPanel(
            "Significant Markers",
            h4("Significant Markers (within Putative Regions)"),
            p(
              "Markers located within putative genomic regions. These markers can be used for follow-up validation or marker design.",
              style = "font-size: 14px; color: gray;"
            ),
            downloadButton("download_marker_table", "Download Significant Markers"),
            br(), br(),
            DTOutput("marker_table")
          ),
          
          tabPanel(
            "Top 5 per Region",
            h4("Top 5 Candidate Markers per Putative Region"),
            p(
              "Top-ranked markers are prioritized by allele-frequency difference within each putative region, with marker-level support shown as additional evidence.",
              style = "font-size: 14px; color: gray;"
            ),
            downloadButton("download_top5", "Download Top 5 per Region"),
            br(), br(),
            DTOutput("top5_table")
          ),
          
          tabPanel(
            "Candidate Gene Annotations",
            h4("Candidate Gene Annotations"),
            p(
              "Nearby gene annotations for top-ranked candidate markers within putative genomic regions.",
              style = "font-size: 14px; color: gray;"
            ),
            downloadButton("download_gene_table", "Download Candidate Gene Table"),
            br(), br(),
            DTOutput("gene_table")
          )
        )
      ),
      
      conditionalPanel(
        condition = "input.analysis_mode == 'classic'",
        tabsetPanel(
          tabPanel(
            "Classic Summary",
            wellPanel(
              h4("Status"),
              textOutput("classic_status")
            ),
            br(),
            h4("Putative Informative Intervals"),
            p(
              "Putative genomic regions supported by informative SNPs under strict parent–bulk matching.",
              style = "font-size: 14px; color: gray;"
            ),
            DTOutput("classic_region_table"),
            br(),
            h4("Marker Classification Summary"),
            DTOutput("classic_summary_table")
          ),
          
          tabPanel(
            "Informative SNPs",
            h4("Informative SNPs"),
            p(
              "SNPs where each bulk matches its corresponding parent within retained putative region(s).",
              style = "font-size: 14px; color: gray;"
            ),
            downloadButton("download_classic_informative", "Download Informative SNPs"),
            br(), br(),
            DTOutput("classic_informative_table")
          ),
          
          tabPanel(
            "All Classified SNPs",
            h4("All Classified SNPs"),
            p(
              "All SNPs classified by strict parent–bulk matching logic.",
              style = "font-size: 14px; color: gray;"
            ),
            downloadButton("download_classic_classified", "Download Classified SNPs"),
            br(), br(),
            DTOutput("classic_classified_table")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # ----------------------------
  # BSAarray mode helpers
  # ----------------------------
  parent_choices <- reactive({
    req(input$parentfile)
    get_parent_sample_names(input$parentfile$datapath)
  })
  
  bulk_choices <- reactive({
    req(input$bulkfile)
    get_bulk_sample_names(input$bulkfile$datapath)
  })
  
  output$bulk_select_ui <- renderUI({
    req(input$analysis_mode == "BSArray", bulk_choices())
    
    bulk1_default <- if (length(bulk_choices()) >= 1) bulk_choices()[1] else NULL
    bulk2_pool <- setdiff(bulk_choices(), bulk1_default)
    bulk2_default <- if (length(bulk2_pool) >= 1) bulk2_pool[1] else NULL
    
    tagList(
      selectInput(
        "bulk1_select",
        "Select Bulk A",
        choices = bulk_choices(),
        selected = bulk1_default
      ),
      selectInput(
        "bulk2_select",
        "Select Bulk B",
        choices = bulk2_pool,
        selected = bulk2_default
      )
    )
  })
  
  output$parent_map_ui <- renderUI({
    req(input$analysis_mode == "BSArray", input$parentfile, input$bulk1_select, input$bulk2_select, parent_choices())
    
    p_choices <- parent_choices()
    
    parent1_default <- if (!is.null(input$parent_bulk1) && input$parent_bulk1 %in% p_choices) {
      input$parent_bulk1
    } else if (length(p_choices) >= 1) {
      p_choices[1]
    } else {
      NULL
    }
    
    parent2_pool <- setdiff(p_choices, parent1_default)
    
    parent2_default <- if (!is.null(input$parent_bulk2) && input$parent_bulk2 %in% parent2_pool) {
      input$parent_bulk2
    } else if (length(parent2_pool) >= 1) {
      parent2_pool[1]
    } else {
      NULL
    }
    
    tagList(
      selectInput(
        "parent_bulk1",
        paste("Parent corresponding to", input$bulk1_select),
        choices = p_choices,
        selected = parent1_default
      ),
      selectInput(
        "parent_bulk2",
        paste("Parent corresponding to", input$bulk2_select),
        choices = parent2_pool,
        selected = parent2_default
      )
    )
  })
  
  observeEvent(input$parent_bulk1, {
    req(input$analysis_mode == "BSArray", parent_choices())
    
    new_parent2_choices <- setdiff(parent_choices(), input$parent_bulk1)
    selected_parent2 <- input$parent_bulk2
    
    if (is.null(selected_parent2) || !(selected_parent2 %in% new_parent2_choices)) {
      selected_parent2 <- if (length(new_parent2_choices) >= 1) new_parent2_choices[1] else NULL
    }
    
    updateSelectInput(
      session,
      "parent_bulk2",
      choices = new_parent2_choices,
      selected = selected_parent2
    )
  })
  
  observeEvent(input$bulk1_select, {
    req(input$analysis_mode == "BSArray", bulk_choices())
    
    new_bulk2_choices <- setdiff(bulk_choices(), input$bulk1_select)
    selected_bulk2 <- input$bulk2_select
    
    if (is.null(selected_bulk2) || !(selected_bulk2 %in% new_bulk2_choices)) {
      selected_bulk2 <- if (length(new_bulk2_choices) >= 1) new_bulk2_choices[1] else NULL
    }
    
    updateSelectInput(
      session,
      "bulk2_select",
      choices = new_bulk2_choices,
      selected = selected_bulk2
    )
  })
  
  # ----------------------------
  # Classic mode helpers
  # ----------------------------
  matrix_bulk_choices <- reactive({
    req(input$bulk_matrix_file)
    get_genomestudio_matrix_sample_names(input$bulk_matrix_file$datapath)
  })
  
  matrix_parent_choices <- reactive({
    req(input$parent_matrix_file)
    get_genomestudio_matrix_sample_names(input$parent_matrix_file$datapath)
  })
  
  output$matrix_bulk_select_ui <- renderUI({
    req(input$analysis_mode == "classic", matrix_bulk_choices())
    
    bulk1_default <- if (length(matrix_bulk_choices()) >= 1) matrix_bulk_choices()[1] else NULL
    bulk2_pool <- setdiff(matrix_bulk_choices(), bulk1_default)
    bulk2_default <- if (length(bulk2_pool) >= 1) bulk2_pool[1] else NULL
    
    tagList(
      selectInput(
        "matrix_bulk1",
        "Select Bulk A",
        choices = matrix_bulk_choices(),
        selected = bulk1_default
      ),
      selectInput(
        "matrix_bulk2",
        "Select Bulk B",
        choices = bulk2_pool,
        selected = bulk2_default
      )
    )
  })
  
  output$matrix_parent_select_ui <- renderUI({
    req(input$analysis_mode == "classic", matrix_parent_choices())
    
    parent1_default <- if (length(matrix_parent_choices()) >= 1) matrix_parent_choices()[1] else NULL
    parent2_pool <- setdiff(matrix_parent_choices(), parent1_default)
    parent2_default <- if (length(parent2_pool) >= 1) parent2_pool[1] else NULL
    
    tagList(
      selectInput(
        "matrix_parent1",
        "Select parent corresponding to Bulk A",
        choices = matrix_parent_choices(),
        selected = parent1_default
      ),
      selectInput(
        "matrix_parent2",
        "Select parent corresponding to Bulk B",
        choices = parent2_pool,
        selected = parent2_default
      )
    )
  })
  
  observeEvent(input$matrix_bulk1, {
    req(input$analysis_mode == "classic", matrix_bulk_choices())
    
    new_bulk2_choices <- setdiff(matrix_bulk_choices(), input$matrix_bulk1)
    selected_bulk2 <- input$matrix_bulk2
    
    if (is.null(selected_bulk2) || !(selected_bulk2 %in% new_bulk2_choices)) {
      selected_bulk2 <- if (length(new_bulk2_choices) >= 1) new_bulk2_choices[1] else NULL
    }
    
    updateSelectInput(
      session,
      "matrix_bulk2",
      choices = new_bulk2_choices,
      selected = selected_bulk2
    )
  })
  
  observeEvent(input$matrix_parent1, {
    req(input$analysis_mode == "classic", matrix_parent_choices())
    
    new_parent2_choices <- setdiff(matrix_parent_choices(), input$matrix_parent1)
    selected_parent2 <- input$matrix_parent2
    
    if (is.null(selected_parent2) || !(selected_parent2 %in% new_parent2_choices)) {
      selected_parent2 <- if (length(new_parent2_choices) >= 1) new_parent2_choices[1] else NULL
    }
    
    updateSelectInput(
      session,
      "matrix_parent2",
      choices = new_parent2_choices,
      selected = selected_parent2
    )
  })
  
  # ----------------------------
  # BSAarray mode analysis
  # ----------------------------
  analysis_results <- eventReactive(input$run_analysis, {
    req(
      input$analysis_mode == "BSArray",
      input$bulkfile,
      input$parentfile,
      input$bulk1_select,
      input$bulk2_select,
      input$parent_bulk1,
      input$parent_bulk2
    )
    
    validate(
      need(input$bulk1_select != input$bulk2_select, "Bulk A and Bulk B must be different samples."),
      need(input$parent_bulk1 != input$parent_bulk2, "Parent assignments for Bulk A and Bulk B must be different.")
    )
    
    tryCatch({
      result <- BSArray_run(
        mode = "BSArray",
        bulkfile = input$bulkfile$datapath,
        parentfile = input$parentfile$datapath,
        snp_db = snp_db,
        bulk1_name = input$bulk1_select,
        bulk2_name = input$bulk2_select,
        parent1_name = input$parent_bulk1,
        parent2_name = input$parent_bulk2,
        window_size = input$window_size,
        step_size = input$step_size,
        min_snp = input$min_snp
      )
      
      showNotification(
        "BSAarray analysis completed successfully!",
        type = "message",
        duration = 8
      )
      
      result
    }, error = function(e) {
      showNotification(
        paste("Error during BSAarray analysis:", e$message),
        type = "error",
        duration = 8
      )
      NULL
    })
  })
  
  # ----------------------------
  # Classic mode analysis
  # ----------------------------
  classic_results <- eventReactive(input$run_classic, {
    req(
      input$analysis_mode == "classic",
      input$bulk_matrix_file,
      input$parent_matrix_file,
      input$matrix_bulk1,
      input$matrix_bulk2,
      input$matrix_parent1,
      input$matrix_parent2
    )
    
    validate(
      need(input$matrix_bulk1 != input$matrix_bulk2, "Bulk A and Bulk B must be different samples."),
      need(input$matrix_parent1 != input$matrix_parent2, "Selected parents must be different.")
    )
    
    tryCatch({
      result <- BSArray_run(
        mode = "classic",
        bulk_matrix_file = input$bulk_matrix_file$datapath,
        parent_matrix_file = input$parent_matrix_file$datapath,
        snp_db = snp_db,
        bulk1_name = input$matrix_bulk1,
        bulk2_name = input$matrix_bulk2,
        parent1_name = input$matrix_parent1,
        parent2_name = input$matrix_parent2,
        classic_matching_model = input$classic_matching_model
      )
      
      showNotification(
        "Classic mode analysis completed successfully!",
        type = "message",
        duration = 4
      )
      
      result
    }, error = function(e) {
      showNotification(
        paste("Error during Classic mode analysis:", e$message),
        type = "error",
        duration = 8
      )
      NULL
    })
  })
  
  # ----------------------------
  # BSAarray mode outputs
  # ----------------------------
  observe({
    req(input$analysis_mode == "BSArray", analysis_results(), input$chr_scope)
    
    results_df <- analysis_results()$Results
    
    all_chrs <- results_df %>%
      dplyr::filter(!is.na(Chr)) %>%
      dplyr::distinct(Chr) %>%
      dplyr::pull(Chr) %>%
      as.character()
    
    high_chrs <- results_df %>%
      dplyr::mutate(is_sig = W_stat >= analysis_results()$Threshold$W) %>%
      dplyr::group_by(Chr) %>%
      dplyr::arrange(Window_Start, .by_group = TRUE) %>%
      dplyr::mutate(run = cumsum(c(1, diff(is_sig) != 0))) %>%
      dplyr::filter(is_sig) %>%
      dplyr::count(Chr, run, name = "n_windows") %>%
      dplyr::group_by(Chr) %>%
      dplyr::summarise(max_consecutive = max(n_windows), .groups = "drop") %>%
      dplyr::filter(max_consecutive >= 2) %>%
      dplyr::pull(Chr) %>%
      as.character()
    
    chr_choices <- switch(
      input$chr_scope,
      "all" = all_chrs,
      "high" = high_chrs
    )
    
    if (length(chr_choices) == 0) {
      updateSelectInput(session, "chr_select", choices = character(0), selected = character(0))
      return()
    }
    
    selected_chr <- input$chr_select
    if (is.null(selected_chr) || !(selected_chr %in% chr_choices)) {
      selected_chr <- chr_choices[1]
    }
    
    updateSelectInput(
      session,
      "chr_select",
      choices = chr_choices,
      selected = selected_chr
    )
  })
  
  output$status <- renderText({
    req(input$analysis_mode == "BSArray")
    
    if (input$run_analysis == 0) {
      "Waiting for input files and Run Analysis."
    } else if (is.null(analysis_results())) {
      "Analysis failed. Please check input files."
    } else {
      paste(
        "Analysis completed successfully.",
        "\nBSA file:", input$bulkfile$name,
        "\nParents file:", input$parentfile$name,
        "\nBulk A:", input$bulk1_select,
        "\nBulk B:", input$bulk2_select,
        "\nParent for Bulk A:", input$parent_bulk1,
        "\nParent for Bulk B:", input$parent_bulk2,
        "\nPrimary statistic: W"
      )
    }
  })
  
  output$chr_summary <- renderUI({
    req(input$analysis_mode == "BSArray", analysis_results())
    
    sig_chrs <- analysis_results()$Results %>%
      dplyr::mutate(is_sig = W_stat >= analysis_results()$Threshold$W) %>%
      dplyr::group_by(Chr) %>%
      dplyr::arrange(Window_Start, .by_group = TRUE) %>%
      dplyr::mutate(run = cumsum(c(1, diff(is_sig) != 0))) %>%
      dplyr::filter(is_sig) %>%
      dplyr::count(Chr, run, name = "n_windows") %>%
      dplyr::group_by(Chr) %>%
      dplyr::summarise(max_consecutive = max(n_windows), .groups = "drop") %>%
      dplyr::filter(max_consecutive >= 2) %>%
      dplyr::pull(Chr) %>%
      as.character()
    
    sig_text <- if (length(sig_chrs) > 0) paste(sig_chrs, collapse = ", ") else "None"
    
    tags$div(
      style = "margin-top: 4px; margin-bottom: 12px;",
      tags$b("Significant chromosomes (≥99%): "),
      tags$span(style = "color:#d62728;", sig_text)
    )
  })
  
  region_df_reactive <- reactive({
    req(input$analysis_mode == "BSArray", analysis_results())
    summarize_putative_regions(analysis_results())
  })
  
  marker_df_reactive <- reactive({
    req(input$analysis_mode == "BSArray", analysis_results())
    
    analysis_results()$Significant_Markers %>%
      dplyr::mutate(
        Region_Start_Mb = round(Region_Start / 1e6, 2),
        Region_End_Mb   = round(Region_End / 1e6, 2),
        Abs_Delta_AF = round(Abs_Delta_AF, 3),
        W = round(W_SNP, 3),
        Window_ID = paste0(sprintf("%.2f", Region_Start_Mb), "-", sprintf("%.2f", Region_End_Mb)),
        Position = format(Position, big.mark = ",")
      ) %>%
      dplyr::select(Chr, Window_ID, SNP_Name, Position, Abs_Delta_AF, W) %>%
      dplyr::rename(
        Chromosome = Chr,
        `Window ID` = Window_ID,
        `Marker ID` = SNP_Name,
        `Position (bp)` = Position,
        `Absolute Delta AF` = Abs_Delta_AF
      )
  })
  
  annotated_top5_reactive <- reactive({
    req(input$analysis_mode == "BSArray", analysis_results())
    
    analysis_results()$Top5_Per_Region %>%
      annotate_markers_with_genes(
        gene_db = gene_db,
        flank_bp = input$flank_bp
      ) %>%
      dplyr::mutate(
        Abs_Delta_AF = round(Abs_Delta_AF, 3),
        W = round(W_SNP, 3),
        Position = format(Position, big.mark = ","),
        Gene_Annotation_Display = clean_gene_display_text(Gene_Annotation_Display)
      ) %>%
      dplyr::group_by(Chr, Region_Start, Region_End) %>%
      dplyr::arrange(desc(Abs_Delta_AF), desc(W), .by_group = TRUE) %>%
      dplyr::slice_head(n = 5) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(Chr, Region_Start, Region_End, SNP_Name, .keep_all = TRUE) %>%
      dplyr::mutate(
        Region_Start_Mb = round(Region_Start / 1e6, 2),
        Region_End_Mb   = round(Region_End / 1e6, 2),
        Window_ID = paste0(sprintf("%.2f", Region_Start_Mb), "-", sprintf("%.2f", Region_End_Mb))
      )
  })
  
  top5_df_reactive <- reactive({
    req(input$analysis_mode == "BSArray", annotated_top5_reactive())
    
    annotated_top5_reactive() %>%
      dplyr::group_by(Chr, Region_Start, Region_End) %>%
      dplyr::arrange(desc(Abs_Delta_AF), desc(W), .by_group = TRUE) %>%
      dplyr::mutate(`Rank in Region` = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::select(`Rank in Region`, Chr, Window_ID, SNP_Name, Position, Abs_Delta_AF, W) %>%
      dplyr::rename(
        Chromosome = Chr,
        `Window ID` = Window_ID,
        `Marker ID` = SNP_Name,
        `Position (bp)` = Position,
        `Absolute Delta AF` = Abs_Delta_AF
      )
  })
  
  gene_df_reactive <- reactive({
    req(input$analysis_mode == "BSArray", annotated_top5_reactive())
    
    annotated_top5_reactive() %>%
      dplyr::mutate(
        Closest_Annotation = stringr::str_replace(Closest_Annotation, ";$", ""),
        Closest_Distance_kb = dplyr::if_else(
          is.na(Closest_Distance_kb),
          NA_character_,
          ifelse(Closest_Distance_kb == 0, "<0.01", sprintf("%.2f", Closest_Distance_kb))
        ),
        Closest_Gene = ifelse(
          is.na(Closest_Gene),
          NA,
          paste0(Closest_Gene)
        )
      ) %>%
      dplyr::select(
        Chr, Window_ID, SNP_Name, Position, Abs_Delta_AF, W,
        Closest_Gene, Closest_Distance_kb, Closest_Annotation
      ) %>%
      dplyr::rename(
        Chromosome = Chr,
        `Window ID` = Window_ID,
        `Marker ID` = SNP_Name,
        `Position (bp)` = Position,
        `Absolute Delta AF` = Abs_Delta_AF,
        `Closest Gene` = Closest_Gene,
        `Distance (kb)` = Closest_Distance_kb,
        `Closest Gene Annotation` = Closest_Annotation
      )
  })
  
  gene_export_reactive <- reactive({
    req(input$analysis_mode == "BSArray", annotated_top5_reactive())
    
    annotated_top5_reactive() %>%
      dplyr::mutate(
        Closest_Annotation = stringr::str_replace(Closest_Annotation, ";$", ""),
        Closest_Distance_kb = dplyr::if_else(
          is.na(Closest_Distance_kb),
          NA_character_,
          ifelse(Closest_Distance_kb == 0, "<0.01", sprintf("%.2f", Closest_Distance_kb))
        )
      ) %>%
      dplyr::select(
        Chr, Window_ID, SNP_Name, Position, Abs_Delta_AF, W,
        Closest_Gene, Closest_Distance_kb, Closest_Annotation
      ) %>%
      dplyr::rename(
        Chromosome = Chr,
        `Window ID` = Window_ID,
        `Marker ID` = SNP_Name,
        `Position (bp)` = Position,
        `Absolute Delta AF` = Abs_Delta_AF,
        `Closest Gene` = Closest_Gene,
        `Distance (kb)` = Closest_Distance_kb,
        `Closest Gene Annotation` = Closest_Annotation
      )
  })
  
  output$chr_plot <- renderPlot({
    req(input$analysis_mode == "BSArray", analysis_results(), input$chr_select)
    
    plot_bsarray_signature_regions(
      results = analysis_results(),
      chromosome = input$chr_select,
      start_mb = NULL,
      end_mb = NULL,
      show_thresholds = TRUE,
      show_regions = TRUE,
      title = paste0(
        "Chromosome ",
        stringr::str_remove(input$chr_select, "^Gm"),
        " Genomic Landscape Plot"
      )
    )
  })
  
  output$sig_plot <- renderPlot({
    req(input$analysis_mode == "BSArray", analysis_results())
    validate(need(!is.null(analysis_results()), "No analysis results available."))
    plot_all_chromosomes_shaded(analysis_results())
  })
  
  output$region_table <- renderDT({
    req(input$analysis_mode == "BSArray")
    region_df <- region_df_reactive()
    
    validate(need(nrow(region_df) > 0, "No putative genomic regions detected."))
    
    DT::datatable(
      region_df,
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip"),
      rownames = FALSE
    )
  })
  
  output$marker_table <- renderDT({
    req(input$analysis_mode == "BSArray")
    marker_df <- marker_df_reactive()
    
    validate(need(nrow(marker_df) > 0, "No significant markers found within putative regions."))
    
    DT::datatable(
      marker_df,
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip"),
      rownames = FALSE
    )
  })
  
  output$top5_table <- renderDT({
    req(input$analysis_mode == "BSArray")
    top5_df <- top5_df_reactive()
    
    validate(need(nrow(top5_df) > 0, "No top markers found in putative regions."))
    
    DT::datatable(
      as.data.frame(top5_df),
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
    )
  })
  
  output$gene_table <- renderDT({
    req(input$analysis_mode == "BSArray")
    gene_table_df <- gene_df_reactive()
    
    validate(need(nrow(gene_table_df) > 0, "No annotated candidate genes found."))
    
    DT::datatable(
      as.data.frame(gene_table_df),
      rownames = FALSE,
      escape = FALSE,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = "tip",
        columnDefs = list(
          list(width = "160px", targets = 2),
          list(width = "90px", targets = 7),
          list(width = "300px", targets = 8)
        )
      )
    )
  })
  
  # ----------------------------
  # Classic mode outputs
  # ----------------------------
  output$classic_status <- renderText({
    req(input$analysis_mode == "classic")
    
    if (input$run_classic == 0) {
      "Waiting for genotype files and Run Classic Mode."
    } else if (is.null(classic_results())) {
      "Classic mode analysis failed. Please check input files."
    } else {
      paste(
        "Classic mode analysis completed successfully.",
        "\nBulk genotype file:", input$bulk_matrix_file$name,
        "\nParent genotype file:", input$parent_matrix_file$name,
        "\nBulk A:", input$matrix_bulk1,
        "\nBulk B:", input$matrix_bulk2,
        "\nParent for Bulk A:", input$matrix_parent1,
        "\nParent for Bulk B:", input$matrix_parent2,
        "\nMatching model:", classic_results()$Matching_Model
      )
    }
  })
  
  output$classic_region_table <- renderDT({
    req(input$analysis_mode == "classic", classic_results())
    
    region_df <- summarize_informative_intervals(
      classic_results(),
      min_snp = input$classic_min_snp
    )
    
    
    validate(
      need(nrow(region_df) > 0, "No informative regions detected after applying the minimum SNP threshold.")
    )
    
    DT::datatable(
      region_df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
    )
  })
  
  output$classic_summary_table <- renderDT({
    req(input$analysis_mode == "classic", classic_results())
    
    res <- classic_results()
    
    df <- res$Classified_Data %>%
      dplyr::count(Marker_Class, name = "Count") %>%
      dplyr::arrange(dplyr::desc(Count))
    
    DT::datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
    )
  })
  
  output$classic_informative_table <- renderDT({
    req(input$analysis_mode == "classic", classic_results())
    
    res <- classic_results()
    
    df <- filter_classic_informative_to_intervals(
      res,
      min_snp = input$classic_min_snp
    ) %>%
      dplyr::select(
        Chr, SNP_Name, Position, Allele1, Allele2,
        R_Bulk, S_Bulk, R_Parent, S_Parent, Marker_Class
      ) %>%
      rename_classic_output_samples(res) %>%
      dplyr::mutate(
        Position = format(Position, big.mark = ",")
      ) %>%
      dplyr::rename(
        Chromosome = Chr,
        `Marker ID` = SNP_Name,
        `Position (bp)` = Position,
        `Marker Class` = Marker_Class
      )
    
    validate(need(nrow(df) > 0, "No informative SNPs found in retained informative intervals."))
    
    DT::datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
    )
  })
  
  output$classic_classified_table <- renderDT({
    req(input$analysis_mode == "classic", classic_results())
    
    res <- classic_results()
    
    df <- res$Classified_Data %>%
      dplyr::select(
        Chr, SNP_Name, Position, Allele1, Allele2,
        R_Bulk, S_Bulk, R_Parent, S_Parent, Marker_Class
      ) %>%
      rename_classic_output_samples(res) %>%
      dplyr::mutate(
        Position = format(Position, big.mark = ",")
      ) %>%
      dplyr::rename(
        Chromosome = Chr,
        `Marker ID` = SNP_Name,
        `Position (bp)` = Position,
        `Marker Class` = Marker_Class
      )
    
    DT::datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
    )
  })
  
  # ----------------------------
  # BSAarray downloads
  # ----------------------------
  output$download_all_chr_plot <- downloadHandler(
    filename = function() {
      paste0("BSAarray_All_Chromosomes_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(input$analysis_mode == "BSArray", analysis_results())
      png(file, width = 1800, height = 1200, res = 180)
      print(plot_all_chromosomes_shaded(analysis_results()))
      dev.off()
    }
  )
  
  output$download_chr_plot <- downloadHandler(
    filename = function() {
      paste0("BSAarray_", input$chr_select, "_Plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(input$analysis_mode == "BSArray", analysis_results(), input$chr_select)
      
png(file, width = 1600, height = 900, res = 150)
print(
  plot_bsarray_signature_regions(
    results = analysis_results(),
    chromosome = input$chr_select,
    start_mb = NULL,
    end_mb = NULL,
    show_thresholds = TRUE,
    show_regions = TRUE,
    title = paste0(
      "Chromosome ",
      sub("^Gm", "", input$chr_select),
      ": Genomic Landscape Plot"
    )
  )
)
dev.off()
    }
  )
  
  output$download_region_table <- downloadHandler(
    filename = function() {
      paste0("Putative_Regions_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(input$analysis_mode == "BSArray")
      readr::write_csv(region_df_reactive(), file)
    }
  )
  
  output$download_marker_table <- downloadHandler(
    filename = function() {
      paste0("Significant_Markers_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(input$analysis_mode == "BSArray")
      readr::write_csv(marker_df_reactive(), file)
    }
  )
  
  output$download_top5 <- downloadHandler(
    filename = function() {
      paste0("Top5_Per_Region_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(input$analysis_mode == "BSArray")
      readr::write_csv(top5_df_reactive(), file)
    }
  )
  
  output$download_gene_table <- downloadHandler(
    filename = function() {
      paste0("Candidate_Gene_Annotations_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(input$analysis_mode == "BSArray")
      readr::write_csv(gene_df_reactive(), file)
    }
  )
  
  output$download_full_results <- downloadHandler(
    filename = function() {
      paste0("BSAarray_Full_Results_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      req(input$analysis_mode == "BSArray", analysis_results())
      
      region_df <- region_df_reactive()
      sig_markers_df <- marker_df_reactive()
      top5_df <- top5_df_reactive()
      gene_table_df <- gene_export_reactive()
      
      wb <- openxlsx::createWorkbook()
      
      openxlsx::addWorksheet(wb, "Putative_Regions")
      openxlsx::writeData(wb, "Putative_Regions", region_df)
      
      openxlsx::addWorksheet(wb, "Significant_Markers")
      openxlsx::writeData(wb, "Significant_Markers", sig_markers_df)
      
      openxlsx::addWorksheet(wb, "Top5_Per_Region")
      openxlsx::writeData(wb, "Top5_Per_Region", top5_df)
      
      openxlsx::addWorksheet(wb, "Candidate_Genes")
      openxlsx::writeData(wb, "Candidate_Genes", gene_table_df)
      
      header_style <- openxlsx::createStyle(
        textDecoration = "bold",
        halign = "center",
        valign = "center"
      )
      
      wrap_style <- openxlsx::createStyle(
        wrapText = TRUE,
        valign = "top"
      )
      
      for (sheet_name in c("Putative_Regions", "Significant_Markers", "Top5_Per_Region", "Candidate_Genes")) {
        df_obj <- switch(
          sheet_name,
          "Putative_Regions"    = region_df,
          "Significant_Markers" = sig_markers_df,
          "Top5_Per_Region"     = top5_df,
          "Candidate_Genes"     = gene_table_df
        )
        
        openxlsx::addStyle(
          wb, sheet = sheet_name,
          style = header_style,
          rows = 1,
          cols = 1:ncol(df_obj),
          gridExpand = TRUE
        )
        
        openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
      }
      
      openxlsx::setColWidths(wb, "Putative_Regions", cols = 1:ncol(region_df), widths = "auto")
      openxlsx::setColWidths(wb, "Significant_Markers", cols = 1:ncol(sig_markers_df), widths = "auto")
      openxlsx::setColWidths(wb, "Top5_Per_Region", cols = 1:ncol(top5_df), widths = "auto")
      openxlsx::setColWidths(wb, "Candidate_Genes", cols = 1:ncol(gene_table_df), widths = "auto")
      
      openxlsx::addStyle(
        wb, "Candidate_Genes",
        style = wrap_style,
        rows = 1:(nrow(gene_table_df) + 1),
        cols = 1:ncol(gene_table_df),
        gridExpand = TRUE,
        stack = TRUE
      )
      
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # ----------------------------
  # Classic mode downloads
  # ----------------------------
  output$download_classic_informative <- downloadHandler(
    filename = function() {
      paste0("Classic_Informative_SNPs_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(input$analysis_mode == "classic", classic_results())
      
      res <- classic_results()
      
      df <- filter_classic_informative_to_intervals(
        res,
        min_snp = input$classic_min_snp
      ) %>%
        rename_classic_output_samples(res)
      
      readr::write_csv(df, file)
    }
  )
  
  output$download_classic_classified <- downloadHandler(
    filename = function() {
      paste0("Classic_Classified_SNPs_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(input$analysis_mode == "classic", classic_results())
      
      res <- classic_results()
      
      df <- res$Classified_Data %>%
        rename_classic_output_samples(res)
      
      readr::write_csv(df, file)
    }
  )
}

shinyApp(ui = ui, server = server)
