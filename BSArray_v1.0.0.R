# SoyBSArray_v1.8.R

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(IRanges)
library(tibble)

# 1) Import GenomeStudio report
read_genomestudio_report <- function(filepath) {
  all_lines <- readLines(filepath, warn = FALSE, encoding = "UTF-8")
  all_lines <- sub("^\ufeff", "", all_lines)
  all_lines_trim <- trimws(all_lines)
  
  data_start <- grep("^\\[Data\\]$", all_lines_trim, ignore.case = TRUE)
  
  if (length(data_start) == 0) {
    message("First 20 lines of uploaded file:")
    print(utils::head(all_lines_trim, 20))
    stop(
      "Could not find [Data] section in GenomeStudio file. ",
      "Please upload the original GenomeStudio text report containing [Header] and [Data]."
    )
  }
  
  data_lines <- all_lines[(data_start[1] + 1):length(all_lines)]
  
  readr::read_tsv(
    file = I(data_lines),
    col_names = TRUE,
    show_col_types = FALSE
  ) %>%
    dplyr::rename(
      SNP_Name    = `SNP Name`,
      Sample_ID   = `Sample ID`,
      Allele1_Top = `Allele1 - Top`,
      Allele2_Top = `Allele2 - Top`,
      GT_Score    = `GT Score`,
      GC_Score    = `GC Score`,
      Theta       = Theta
    )
}

# 2) Import parental genotype file
read_parent_file <- function(filepath, parent1_name = NULL, parent2_name = NULL) {
  all_lines <- readLines(filepath, warn = FALSE, encoding = "UTF-8")
  all_lines <- sub("^\ufeff", "", all_lines)
  all_lines_trim <- trimws(all_lines)
  
  data_start <- grep("^\\[Data\\]$", all_lines_trim, ignore.case = TRUE)
  
  if (length(data_start) > 0) {
    data_lines <- all_lines[(data_start[1] + 1):length(all_lines)]
    
    raw_df <- readr::read_tsv(
      file = I(data_lines),
      col_names = TRUE,
      show_col_types = FALSE
    )
    
    names(raw_df)[1] <- "SNP_Name"
    available_parents <- setdiff(names(raw_df), "SNP_Name")
    
    if (is.null(parent1_name) || is.null(parent2_name)) {
      stop(
        "Please provide parent1_name and parent2_name. Available parent columns: ",
        paste(available_parents, collapse = ", ")
      )
    }
    
    missing_parents <- setdiff(c(parent1_name, parent2_name), available_parents)
    if (length(missing_parents) > 0) {
      stop("Selected parent(s) not found: ", paste(missing_parents, collapse = ", "))
    }
    
    return(
      raw_df %>%
        dplyr::select(
          SNP_Name,
          R_Parent = all_of(parent1_name),
          S_Parent = all_of(parent2_name)
        )
    )
  }
  
  df <- readr::read_tsv(filepath, show_col_types = FALSE)
  
  names(df) <- trimws(names(df))
  names(df) <- gsub("\\s+", "_", names(df))
  
  if ("Marker" %in% names(df)) {
    names(df)[names(df) == "Marker"] <- "SNP_Name"
  }
  
  if ("S_parent" %in% names(df) && !"S_Parent" %in% names(df)) {
    names(df)[names(df) == "S_parent"] <- "S_Parent"
  }
  
  required_cols <- c("SNP_Name", "R_Parent", "S_Parent")
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop("Parent file is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  df
}

# 3) Helper: get parent sample names
get_parent_sample_names <- function(filepath) {
  all_lines <- readLines(filepath, warn = FALSE, encoding = "UTF-8")
  all_lines <- sub("^\ufeff", "", all_lines)
  all_lines_trim <- trimws(all_lines)
  
  data_start <- grep("^\\[Data\\]$", all_lines_trim, ignore.case = TRUE)
  
  if (length(data_start) > 0) {
    data_lines <- all_lines[(data_start[1] + 1):length(all_lines)]
    
    raw_df <- readr::read_tsv(
      file = I(data_lines),
      col_names = TRUE,
      show_col_types = FALSE
    )
    
    names(raw_df)[1] <- "SNP_Name"
    return(setdiff(names(raw_df), "SNP_Name"))
  }
  
  c("R_Parent", "S_Parent")
}

# 4) Helper: get bulk sample names
get_bulk_sample_names <- function(filepath) {
  data <- read_genomestudio_report(filepath)
  sort(unique(data$Sample_ID))
}

# 5) Estimate Theta-based allele frequency
estimate_AF_theta <- function(data) {
  data %>%
    dplyr::mutate(
      Allele_Freq_A = dplyr::if_else(is.na(Theta), NA_real_, 1 - Theta)
    )
}

# 6) Summarize bulks
summarize_bulks <- function(data, bulk1_name, bulk2_name) {
  if (bulk1_name == bulk2_name) {
    stop("bulk1_name and bulk2_name must be different.")
  }
  
  available_samples <- unique(data$Sample_ID)
  missing_samples <- setdiff(c(bulk1_name, bulk2_name), available_samples)
  
  if (length(missing_samples) > 0) {
    stop("Selected bulk sample(s) not found: ", paste(missing_samples, collapse = ", "))
  }
  
  data %>%
    dplyr::filter(Sample_ID %in% c(bulk1_name, bulk2_name)) %>%
    dplyr::group_by(SNP_Name, Sample_ID) %>%
    dplyr::summarise(
      Mean_AF = mean(Allele_Freq_A, na.rm = TRUE),
      Mean_GC = mean(GC_Score, na.rm = TRUE),
      Mean_GT = mean(GT_Score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Sample_ID = dplyr::case_when(
        Sample_ID == bulk1_name ~ "Bulk_1",
        Sample_ID == bulk2_name ~ "Bulk_2",
        TRUE ~ Sample_ID
      )
    ) %>%
    tidyr::pivot_wider(
      names_from = Sample_ID,
      values_from = c(Mean_AF, Mean_GC, Mean_GT),
      names_sep = "_"
    ) %>%
    dplyr::mutate(
      Bulk_1    = .data$Mean_AF_Bulk_1,
      Bulk_2    = .data$Mean_AF_Bulk_2,
      GC_Bulk_1 = .data$Mean_GC_Bulk_1,
      GC_Bulk_2 = .data$Mean_GC_Bulk_2,
      GT_Bulk_1 = .data$Mean_GT_Bulk_1,
      GT_Bulk_2 = .data$Mean_GT_Bulk_2
    )
}

merge_parents_with_bulks <- function(bulks, parents) {
  dplyr::left_join(bulks, parents, by = "SNP_Name")
}

filter_monomorphic <- function(data) {
  data %>%
    dplyr::filter(!is.na(R_Parent), !is.na(S_Parent)) %>%
    dplyr::filter(R_Parent != S_Parent) %>%
    dplyr::filter(!is.na(Bulk_1), !is.na(Bulk_2)) %>%
    dplyr::filter(Bulk_1 != Bulk_2)
}

calculate_delta_AF_and_weights <- function(data, min_weight = 0.01) {
  data %>%
    dplyr::mutate(
      Delta_AF = Bulk_1 - Bulk_2,
      Abs_Delta_AF = abs(Delta_AF),
      Mean_GC = rowMeans(cbind(GC_Bulk_1, GC_Bulk_2), na.rm = TRUE),
      Mean_GT = rowMeans(cbind(GT_Bulk_1, GT_Bulk_2), na.rm = TRUE),
      Weight_Hybrid = Mean_GT * Mean_GC,
      Weight_Hybrid = dplyr::if_else(
        is.na(Weight_Hybrid),
        min_weight,
        pmax(Weight_Hybrid, min_weight)
      ),
      W_SNP = Abs_Delta_AF * Weight_Hybrid
    )
}

parse_snp_names <- function(data, snp_db) {
  if (missing(snp_db) || is.null(snp_db)) {
    stop("snp_db was not provided.")
  }
  
  snp_lookup <- snp_db %>%
    dplyr::select(SNP_Name, Position_v6) %>%
    dplyr::distinct(SNP_Name, .keep_all = TRUE)
  
  data %>%
    dplyr::filter(stringr::str_detect(SNP_Name, "^Gm\\d{1,2}_")) %>%
    tidyr::separate(
      SNP_Name,
      into = c("Chr", "Position", "Allele1", "Allele2"),
      sep = "_",
      remove = FALSE
    ) %>%
    dplyr::mutate(
      Position = as.numeric(Position),
      Chr = factor(Chr, levels = paste0("Gm", stringr::str_pad(1:20, 2, pad = "0")))
    ) %>%
    dplyr::left_join(snp_lookup, by = "SNP_Name") %>%
    dplyr::mutate(Position = dplyr::coalesce(Position_v6, Position)) %>%
    dplyr::select(-Position_v6)
}

sliding_window_smooth <- function(data, window_size, step_size, min_snp) {
  out <- data %>%
    dplyr::filter(!is.na(Chr), !is.na(Position), !is.na(Delta_AF), !is.na(Abs_Delta_AF)) %>%
    dplyr::group_by(Chr) %>%
    dplyr::arrange(Position, .by_group = TRUE) %>%
    dplyr::group_modify(~ {
      chr_data <- .x
      min_pos <- min(chr_data$Position, na.rm = TRUE)
      max_pos <- max(chr_data$Position, na.rm = TRUE)
      window_starts <- seq(min_pos, max_pos, by = step_size)
      
      tibble::tibble(
        Window_Start = window_starts,
        Window_End = window_starts + window_size,
        Mean_Delta_AF = purrr::map2_dbl(window_starts, window_starts + window_size, ~ {
          idx <- chr_data$Position >= .x & chr_data$Position < .y
          if (sum(idx, na.rm = TRUE) == 0) return(NA_real_)
          mean(chr_data$Delta_AF[idx], na.rm = TRUE)
        }),
        Mean_Abs_Delta_AF = purrr::map2_dbl(window_starts, window_starts + window_size, ~ {
          idx <- chr_data$Position >= .x & chr_data$Position < .y
          if (sum(idx, na.rm = TRUE) == 0) return(NA_real_)
          mean(chr_data$Abs_Delta_AF[idx], na.rm = TRUE)
        }),
        W_stat = purrr::map2_dbl(window_starts, window_starts + window_size, ~ {
          idx <- chr_data$Position >= .x & chr_data$Position < .y
          if (sum(idx, na.rm = TRUE) == 0) return(NA_real_)
          stats::weighted.mean(chr_data$Abs_Delta_AF[idx], w = chr_data$Weight_Hybrid[idx], na.rm = TRUE)
        }),
        Mean_GC = purrr::map2_dbl(window_starts, window_starts + window_size, ~ {
          idx <- chr_data$Position >= .x & chr_data$Position < .y
          if (sum(idx, na.rm = TRUE) == 0) return(NA_real_)
          mean(chr_data$Mean_GC[idx], na.rm = TRUE)
        }),
        SNP_Count = purrr::map2_int(
          window_starts,
          window_starts + window_size,
          ~ sum(chr_data$Position >= .x & chr_data$Position < .y, na.rm = TRUE)
        )
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::filter(SNP_Count >= min_snp)
  
  if (nrow(out) == 0) {
    stop("No valid windows remained after smoothing. Try reducing min_snp.")
  }
  
  out
}

calculate_threshold <- function(data, cutoff_high = 0.99, cutoff_low = 0.95) {
  list(
    W = unname(quantile(data$W_stat, probs = cutoff_high, na.rm = TRUE)),
    W_low = unname(quantile(data$W_stat, probs = cutoff_low, na.rm = TRUE)),
    cutoff = cutoff_high
  )
}

get_shaded_regions <- function(results_df, threshold, stat_col = "W_stat") {
  stat_sym <- rlang::sym(stat_col)
  
  out <- results_df %>%
    dplyr::filter(!is.na(!!stat_sym), !!stat_sym >= threshold) %>%
    dplyr::group_by(Chr) %>%
    dplyr::arrange(Window_Start, .by_group = TRUE) %>%
    dplyr::group_modify(~ {
      merged_ranges <- IRanges::reduce(
        IRanges::IRanges(start = .x$Window_Start, end = .x$Window_End)
      )
      
      tibble::tibble(
        Region_Start = IRanges::start(merged_ranges),
        Region_End = IRanges::end(merged_ranges)
      )
    }) %>%
    dplyr::ungroup()
  
  if (nrow(out) == 0) {
    return(tibble::tibble(Chr = character(), Region_Start = numeric(), Region_End = numeric()))
  }
  
  out
}

extract_significant_markers <- function(raw_data, regions_df, top_n = NULL) {
  empty_tbl <- tibble::tibble(
    Chr = character(),
    Region_Start = numeric(),
    Region_End = numeric(),
    SNP_Name = character(),
    Position = numeric(),
    Delta_AF = numeric(),
    Abs_Delta_AF = numeric(),
    Weight_Hybrid = numeric(),
    W_SNP = numeric(),
    Allele1 = character(),
    Allele2 = character(),
    Bulk_1 = numeric(),
    Bulk_2 = numeric()
  )
  
  if (is.null(regions_df) || nrow(regions_df) == 0) return(empty_tbl)
  
  markers <- raw_data %>%
    dplyr::filter(!is.na(Position), !is.na(Chr), !is.na(W_SNP)) %>%
    dplyr::inner_join(regions_df, by = "Chr", relationship = "many-to-many") %>%
    dplyr::filter(Position >= Region_Start & Position <= Region_End) %>%
    dplyr::distinct(Chr, Region_Start, Region_End, SNP_Name, .keep_all = TRUE) %>%
    dplyr::select(
      Chr, Region_Start, Region_End, SNP_Name, Position,
      Delta_AF, Abs_Delta_AF, Weight_Hybrid, W_SNP,
      Allele1, Allele2, Bulk_1, Bulk_2
    )
  
  if (nrow(markers) == 0) return(empty_tbl)
  
  if (!is.null(top_n)) {
    markers <- markers %>%
      dplyr::group_by(Chr, Region_Start, Region_End) %>%
      dplyr::arrange(desc(Abs_Delta_AF), desc(W_SNP), Position, .by_group = TRUE) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::ungroup()
  } else {
    markers <- markers %>%
      dplyr::arrange(Chr, Region_Start, Region_End, desc(Abs_Delta_AF), Position)
  }
  
  markers
}

plot_all_chromosomes_shaded <- function(results,
                                        title = "Genome-wide distribution of W statistic across soybean chromosomes") {
  plot_data <- results$Results %>%
    dplyr::arrange(Chr, Window_Start) %>%
    dplyr::group_by(Chr) %>%
    dplyr::mutate(
      Window_Start_Mb = Window_Start / 1e6,
      Above_Threshold = W_stat >= results$Threshold$W,
      Plot_Gap = Window_Start - dplyr::lag(Window_Start),
      sig_new_group = dplyr::case_when(
        !Above_Threshold ~ 0L,
        is.na(Plot_Gap) ~ 1L,
        !dplyr::lag(Above_Threshold, default = FALSE) ~ 1L,
        Plot_Gap > 1e6 ~ 1L,
        TRUE ~ 0L
      ),
      Sig_Group = cumsum(sig_new_group),
      Sig_Group = dplyr::if_else(Above_Threshold, Sig_Group, NA_integer_)
    ) %>%
    dplyr::ungroup()
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = Window_Start_Mb, y = W_stat)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = 0, ymax = W_stat, group = Chr),
      fill = "grey80",
      alpha = 0.5
    ) +
    ggplot2::geom_ribbon(
      data = plot_data %>% dplyr::filter(Above_Threshold),
      ggplot2::aes(ymin = 0, ymax = W_stat, group = interaction(Chr, Sig_Group)),
      fill = "#9467bd",
      alpha = 0.25
    ) +
    ggplot2::geom_line(ggplot2::aes(group = Chr), color = "black", linewidth = 1) +
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = results$Threshold$W, color = "thr_99"),
      linetype = "dashed",
      linewidth = 1
    ) +
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = results$Threshold$W_low, color = "thr_95"),
      linetype = "dashed",
      linewidth = 0.7
    ) +
    ggplot2::scale_color_manual(
      name = "Threshold",
      values = c("thr_99" = "#D62728", "thr_95" = "#1f77b4"),
      labels = c("thr_99" = "99%", "thr_95" = "95%")
    ) +
    ggplot2::facet_wrap(~Chr, scales = "free_x") +
    ggplot2::labs(title = title, x = "Genomic Position (Mb)", y = "W statistic") +
    ggplot2::theme_minimal(base_size = 20) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

plot_bsarray_signature_regions <- function(results,
                                           chromosome = NULL,
                                           start_mb = NULL,
                                           end_mb = NULL,
                                           show_thresholds = TRUE,
                                           show_regions = TRUE,
                                           title = NULL) {
  if (is.null(chromosome)) stop("Please provide a chromosome.")
  
  plot_data <- results$Results %>%
    dplyr::filter(Chr == chromosome) %>%
    dplyr::arrange(Window_Start) %>%
    dplyr::mutate(
      Position_Mb = Window_Start / 1e6,
      Above_Threshold = W_stat >= results$Threshold$W,
      sig_gap_bp = Window_Start - dplyr::lag(Window_Start),
      sig_new_group = dplyr::case_when(
        !Above_Threshold ~ 0L,
        is.na(sig_gap_bp) ~ 1L,
        !dplyr::lag(Above_Threshold, default = FALSE) ~ 1L,
        sig_gap_bp > 1e6 ~ 1L,
        TRUE ~ 0L
      ),
      Sig_Group = cumsum(sig_new_group),
      Sig_Group = dplyr::if_else(Above_Threshold, Sig_Group, NA_integer_)
    )
  
  if (!is.null(start_mb) && !is.null(end_mb)) {
    plot_data <- plot_data %>%
      dplyr::filter(Position_Mb >= start_mb, Position_Mb <= end_mb)
  }
  
  if (is.null(title)) title <- paste0(chromosome, ": BSAarray Plot")
  
  threshold_df <- tibble::tibble(
    Threshold = factor(c("99%", "95%"), levels = c("99%", "95%")),
    yint = c(results$Threshold$W, results$Threshold$W_low),
    line_wd = c(1.2, 0.8)
  )
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Position_Mb, y = W_stat)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = 0, ymax = W_stat),
      fill = "grey80",
      alpha = 0.5
    ) +
    ggplot2::geom_ribbon(
      data = plot_data %>% dplyr::filter(Above_Threshold),
      ggplot2::aes(ymin = 0, ymax = W_stat, group = Sig_Group),
      fill = "#9467bd",
      alpha = 0.25
    ) +
    ggplot2::geom_line(color = "black", linewidth = 1.2)
  
  if (show_thresholds) {
    p <- p +
      ggplot2::geom_hline(
        data = threshold_df,
        ggplot2::aes(yintercept = yint, color = Threshold, linewidth = line_wd),
        linetype = "dashed",
        show.legend = TRUE
      ) +
      ggplot2::scale_color_manual(
        name = "Threshold",
        values = c("99%" = "#D62728", "95%" = "#1f77b4")
      ) +
      ggplot2::scale_linewidth_identity()
  }
  
  p +
    ggplot2::labs(title = title, x = "Position (Mb)", y = "W statistic") +
    ggplot2::theme_minimal(base_size = 25) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}

annotate_markers_with_genes <- function(marker_df, gene_db, flank_bp = 10000) {
  if (is.null(marker_df) || !is.data.frame(marker_df) || nrow(marker_df) == 0) {
    return(marker_df)
  }
  
  gene_db_sub <- gene_db %>%
    dplyr::filter(Chromosome %in% marker_df$Chr)
  
  annotated_list <- purrr::map_dfr(seq_len(nrow(marker_df)), function(i) {
    this_chr <- marker_df$Chr[i]
    this_pos <- marker_df$Position[i]
    this_snp <- marker_df$SNP_Name[i]
    
    hits <- gene_db_sub %>%
      dplyr::filter(
        Chromosome == this_chr,
        Start <= this_pos + flank_bp,
        End >= this_pos - flank_bp
      ) %>%
      dplyr::mutate(
        Distance_bp = dplyr::case_when(
          this_pos >= Start & this_pos <= End ~ 0,
          TRUE ~ pmin(abs(this_pos - Start), abs(this_pos - End))
        )
      )
    
    if (nrow(hits) == 0) {
      tibble::tibble(
        SNP_Name = this_snp,
        Nearby_Genes = NA_character_,
        Nearby_Annotations = NA_character_,
        Closest_Gene = NA_character_,
        Closest_Annotation = NA_character_,
        Closest_Distance_kb = NA_real_,
        Gene_Annotation_Display = NA_character_
      )
    } else {
      closest <- hits %>%
        dplyr::arrange(Distance_bp) %>%
        dplyr::slice(1)
      
      pairs <- hits %>%
        dplyr::transmute(Gene_model, Annotation = Annotations_clean) %>%
        dplyr::distinct()
      
      display_text <- pairs %>%
        dplyr::mutate(
          Label = ifelse(
            Gene_model == closest$Gene_model[1],
            paste0("<b>", Gene_model, " — ", Annotation, "</b>"),
            paste0(Gene_model, " — ", Annotation)
          )
        ) %>%
        dplyr::pull(Label) %>%
        paste(collapse = "<br>")
      
      tibble::tibble(
        SNP_Name = this_snp,
        Nearby_Genes = paste(unique(hits$Gene_model), collapse = "; "),
        Nearby_Annotations = paste(unique(hits$Annotations_clean), collapse = "; "),
        Closest_Gene = closest$Gene_model,
        Closest_Annotation = closest$Annotations_clean,
        Closest_Distance_kb = round(closest$Distance_bp / 1000, 2),
        Gene_Annotation_Display = display_text
      )
    }
  }) %>%
    dplyr::distinct(SNP_Name, .keep_all = TRUE)
  
  marker_df %>%
    dplyr::left_join(annotated_list, by = "SNP_Name")
}

BSArray_run_W <- function(bulkfile, parentfile, snp_db,
                          bulk1_name, bulk2_name,
                          parent1_name = NULL, parent2_name = NULL,
                          window_size = 1e6,
                          step_size = 1e6,
                          min_snp = 10,
                          cutoff = 0.99) {
  mydata <- read_genomestudio_report(bulkfile)
  
  parents <- read_parent_file(
    filepath = parentfile,
    parent1_name = parent1_name,
    parent2_name = parent2_name
  )
  
  mydata_AF <- estimate_AF_theta(mydata)
  AF_summary <- summarize_bulks(mydata_AF, bulk1_name, bulk2_name)
  AF_merged <- merge_parents_with_bulks(AF_summary, parents)
  AF_filtered <- filter_monomorphic(AF_merged)
  Stage1_Delta <- calculate_delta_AF_and_weights(AF_filtered)
  
  Stage1_valid <- parse_snp_names(Stage1_Delta, snp_db) %>%
    dplyr::filter(!is.na(Chr))
  
  if (nrow(Stage1_valid) == 0) stop("No valid SNPs remained after filtering and parsing.")
  
  Stage1_smoothed <- sliding_window_smooth(Stage1_valid, window_size, step_size, min_snp)
  Threshold <- calculate_threshold(Stage1_smoothed, cutoff_high = cutoff)
  
  Shaded_Regions <- get_shaded_regions(Stage1_smoothed, Threshold$W, stat_col = "W_stat")
  
  Significant_Markers <- extract_significant_markers(Stage1_valid, Shaded_Regions)
  Top5_Per_Region <- extract_significant_markers(Stage1_valid, Shaded_Regions, top_n = 5)
  
  list(
    Results = Stage1_smoothed,
    Threshold = Threshold,
    Filtered_Data = AF_filtered,
    Parsed_Data = Stage1_valid,
    Shaded_Regions = Shaded_Regions,
    Significant_Markers = Significant_Markers,
    Top5_Per_Region = Top5_Per_Region,
    Primary_Statistic = "W_stat",
    Secondary_Statistic = "Mean_Abs_Delta_AF",
    Bulk_1_Name = bulk1_name,
    Bulk_2_Name = bulk2_name,
    Parent_1_Name = parent1_name,
    Parent_2_Name = parent2_name
  )
}

# Classic mode
read_genomestudio_matrix_report <- function(filepath,
                                            sample1_name = NULL,
                                            sample2_name = NULL,
                                            sample1_label = "Sample_1",
                                            sample2_label = "Sample_2") {
  all_lines <- readLines(filepath, warn = FALSE, encoding = "UTF-8")
  all_lines <- sub("^\ufeff", "", all_lines)
  all_lines_trim <- trimws(all_lines)
  
  data_start <- grep("^\\[Data\\]$", all_lines_trim, ignore.case = TRUE)
  if (length(data_start) == 0) stop("Could not find [Data] section in GenomeStudio matrix file.")
  
  data_lines <- all_lines[(data_start[1] + 1):length(all_lines)]
  
  raw_df <- readr::read_tsv(file = I(data_lines), col_names = TRUE, show_col_types = FALSE)
  
  names(raw_df)[1] <- "SNP_Name"
  available_samples <- setdiff(names(raw_df), "SNP_Name")
  
  missing_samples <- setdiff(c(sample1_name, sample2_name), available_samples)
  if (length(missing_samples) > 0) {
    stop("Selected sample(s) not found: ", paste(missing_samples, collapse = ", "))
  }
  
  raw_df %>%
    dplyr::select(
      SNP_Name,
      !!sample1_label := all_of(sample1_name),
      !!sample2_label := all_of(sample2_name)
    ) %>%
    dplyr::mutate(dplyr::across(-SNP_Name, ~ toupper(trimws(.x))))
}

get_genomestudio_matrix_sample_names <- function(filepath) {
  all_lines <- readLines(filepath, warn = FALSE, encoding = "UTF-8")
  all_lines <- sub("^\ufeff", "", all_lines)
  all_lines_trim <- trimws(all_lines)
  
  data_start <- grep("^\\[Data\\]$", all_lines_trim, ignore.case = TRUE)
  if (length(data_start) == 0) stop("Could not find [Data] section in GenomeStudio matrix file.")
  
  data_lines <- all_lines[(data_start[1] + 1):length(all_lines)]
  
  raw_df <- readr::read_tsv(file = I(data_lines), col_names = TRUE, show_col_types = FALSE)
  names(raw_df)[1] <- "SNP_Name"
  
  setdiff(names(raw_df), "SNP_Name")
}

standardize_genotype <- function(gt) {
  gt <- toupper(trimws(gt))
  
  sapply(gt, function(x) {
    if (is.na(x) || x %in% c("", "-", "--", "NA", "#NAME?")) return(NA_character_)
    if (nchar(x) != 2) return(NA_character_)
    alleles <- strsplit(x, "")[[1]]
    paste(sort(alleles), collapse = "")
  }) %>%
    unname()
}

is_homozygous_gt <- function(gt) {
  !is.na(gt) & nchar(gt) == 2 & substr(gt, 1, 1) == substr(gt, 2, 2)
}

allele_in_genotype <- function(gt, allele) {
  gt <- as.character(gt)
  allele <- as.character(allele)
  
  mapply(
    function(g, a) {
      if (is.na(g) || is.na(a) || g == "" || a == "") return(FALSE)
      grepl(a, g, fixed = TRUE)
    },
    gt,
    allele,
    USE.NAMES = FALSE
  )
}

classify_matrix_snps_strict <- function(data,
                                        matching_model = c(
                                          "strict",
                                          "dominant_favorable",
                                          "recessive_favorable"
                                        )) {
  matching_model <- match.arg(matching_model)
  
  data <- data %>%
    dplyr::mutate(
      R_Bulk   = standardize_genotype(R_Bulk),
      S_Bulk   = standardize_genotype(S_Bulk),
      R_Parent = standardize_genotype(R_Parent),
      S_Parent = standardize_genotype(S_Parent),
      R_parent_allele = substr(R_Parent, 1, 1),
      S_parent_allele = substr(S_Parent, 1, 1)
    )
  
  if (matching_model == "strict") {
    data <- data %>%
      dplyr::mutate(
        R_Bulk_match = R_Bulk == R_Parent,
        S_Bulk_match = S_Bulk == S_Parent
      )
  }
  
  if (matching_model == "dominant_favorable") {
    data <- data %>%
      dplyr::mutate(
        R_Bulk_match = allele_in_genotype(R_Bulk, R_parent_allele),
        S_Bulk_match = S_Bulk == S_Parent
      )
  }
  
  if (matching_model == "recessive_favorable") {
    data <- data %>%
      dplyr::mutate(
        R_Bulk_match = R_Bulk == R_Parent,
        S_Bulk_match = allele_in_genotype(S_Bulk, S_parent_allele)
      )
  }
  
  data %>%
    dplyr::mutate(
      Marker_Class = dplyr::case_when(
        is.na(R_Parent) | is.na(S_Parent) ~ "Missing parent genotype",
        is.na(R_Bulk) | is.na(S_Bulk) ~ "Missing bulk genotype",
        R_Parent == S_Parent ~ "Monomorphic parents",
        
        !is_homozygous_gt(R_Parent) ~ "Heterozygous Parent 1",
        !is_homozygous_gt(S_Parent) ~ "Heterozygous Parent 2",
        
        matching_model == "strict" & !is_homozygous_gt(R_Bulk) ~ "Heterozygous Bulk 1",
        matching_model == "strict" & !is_homozygous_gt(S_Bulk) ~ "Heterozygous Bulk 2",
        
        !S_Bulk_match ~ "Bulk 2 does not match Parent 2",
        !R_Bulk_match ~ "Bulk 1 does not match Parent 1",
        
        matching_model == "dominant_favorable" & !is_homozygous_gt(R_Bulk) ~
          "Informative - heterozygous favorable bulk",
        
        matching_model == "recessive_favorable" & !is_homozygous_gt(S_Bulk) ~
          "Informative - heterozygous contrasting bulk",
        
        TRUE ~ "Informative"
      )
    ) %>%
    dplyr::select(
      -R_parent_allele,
      -S_parent_allele,
      -R_Bulk_match,
      -S_Bulk_match
    )
}

rename_classic_output_samples <- function(df,
                                          bulk1_name,
                                          bulk2_name,
                                          parent1_name,
                                          parent2_name) {
  df %>%
    dplyr::rename(
      !!bulk1_name := R_Bulk,
      !!bulk2_name := S_Bulk,
      !!parent1_name := R_Parent,
      !!parent2_name := S_Parent
    )
}

parse_snp_names_matrix <- function(data) {
  data %>%
    dplyr::filter(stringr::str_detect(SNP_Name, "^Gm\\d{1,2}_")) %>%
    tidyr::separate(
      SNP_Name,
      into = c("Chr", "Position", "Allele1", "Allele2"),
      sep = "_",
      remove = FALSE
    ) %>%
    dplyr::mutate(
      Position = suppressWarnings(as.numeric(Position)),
      Chr = factor(Chr, levels = paste0("Gm", stringr::str_pad(1:20, 2, pad = "0"))),
      Allele1 = toupper(trimws(Allele1)),
      Allele2 = toupper(trimws(Allele2))
    )
}

BSArray_run_classic <- function(bulk_matrix_file,
                                parent_matrix_file,
                                bulk1_name,
                                bulk2_name,
                                parent1_name,
                                parent2_name,
                                snp_db = NULL,
                                matching_model = "strict") {
  bulk_matrix <- read_genomestudio_matrix_report(
    filepath = bulk_matrix_file,
    sample1_name = bulk1_name,
    sample2_name = bulk2_name,
    sample1_label = "R_Bulk",
    sample2_label = "S_Bulk"
  )
  
  parent_matrix <- read_genomestudio_matrix_report(
    filepath = parent_matrix_file,
    sample1_name = parent1_name,
    sample2_name = parent2_name,
    sample1_label = "R_Parent",
    sample2_label = "S_Parent"
  )
  
  merged <- bulk_matrix %>%
    dplyr::left_join(parent_matrix, by = "SNP_Name")
  
  if (!is.null(snp_db)) {
    merged <- merged %>% parse_snp_names(snp_db = snp_db)
  } else {
    merged <- merged %>% parse_snp_names_matrix()
  }
  
  classified <- merged %>%
    classify_matrix_snps_strict(matching_model = matching_model)
  
  informative <- classified %>%
    dplyr::filter(stringr::str_starts(Marker_Class, "Informative")) %>%
    dplyr::arrange(Chr, Position)
  
  chr_summary <- classified %>%
    dplyr::count(Chr, Marker_Class, name = "Count") %>%
    dplyr::arrange(Chr, dplyr::desc(Count))
  
  list(
    Bulk_Matrix = bulk_matrix,
    Parent_Matrix = parent_matrix,
    Merged_Data = merged,
    Classified_Data = classified,
    Informative_SNPs = informative,
    Chromosome_Summary = chr_summary,
    Mode = "classic",
    Matching_Model = matching_model,
    Bulk_1_Name = bulk1_name,
    Bulk_2_Name = bulk2_name,
    Parent_1_Name = parent1_name,
    Parent_2_Name = parent2_name
  )
}

BSArray_run <- function(mode = c("BSArray", "classic"),
                        bulkfile = NULL,
                        parentfile = NULL,
                        bulk_matrix_file = NULL,
                        parent_matrix_file = NULL,
                        snp_db = NULL,
                        bulk1_name = NULL,
                        bulk2_name = NULL,
                        parent1_name = NULL,
                        parent2_name = NULL,
                        window_size = 1e6,
                        step_size = 1e6,
                        min_snp = 10,
                        cutoff = 0.99,
                        classic_matching_model = "strict") {
  mode <- match.arg(mode)
  
  if (mode == "BSArray") {
    return(
      BSArray_run_W(
        bulkfile = bulkfile,
        parentfile = parentfile,
        snp_db = snp_db,
        bulk1_name = bulk1_name,
        bulk2_name = bulk2_name,
        parent1_name = parent1_name,
        parent2_name = parent2_name,
        window_size = window_size,
        step_size = step_size,
        min_snp = min_snp,
        cutoff = cutoff
      )
    )
  }
  
  if (mode == "classic") {
    return(
      BSArray_run_classic(
        bulk_matrix_file = bulk_matrix_file,
        parent_matrix_file = parent_matrix_file,
        bulk1_name = bulk1_name,
        bulk2_name = bulk2_name,
        parent1_name = parent1_name,
        parent2_name = parent2_name,
        snp_db = snp_db,
        matching_model = classic_matching_model
      )
    )
  }
}
