# Comprehensive Immunoassay Analysis Framework
# Author: Hamid Khoshfekr Rudsari
# Department of Biostatistics, MD Anderson Cancer Center
# 
# A comprehensive R framework for multi-cancer biomarker panel analysis
# 
# This framework provides tools for:
# - Individual biomarker analysis with ROC curves and performance metrics
# - Multi-marker panel analysis for improved cancer detection
# - Comprehensive statistical reporting and visualization
# - Parallel processing capabilities for large datasets
#
# Version: 1.0
# License: MIT

# main function ----
comprehensive_immunoassay_analysis <- function(camp_data, camp_key, 
                                               analysis_type = "individual", # "individual", "panel", or "both"
                                               plate_code = "All", # "1-6", "7-8", "1-8", "All"
                                               camp_six = FALSE, # TRUE/FALSE
                                               menopause_status = "none", # pre-menopause, post-menopause, none
                                               normalization = TRUE, # TRUE/FALSE
                                               markers_count = 13,  # 13, 10
                                               new_breast_panel = FALSE, # TRUE (includes alternative models), FALSE
                                               all_specimen = "One Year", # "One Year", "Two Years", "Two Plus", "All"
                                               file_stem,
                                               saving_path = "output/",
                                               protocol_plot = "All", # MERIT_LEAP/All
                                               n_cores = NULL,
                                               use_parallel = TRUE) {
  
  # Load required libraries for parallel processing
  if (use_parallel) {
    if (!require(foreach, quietly = TRUE)) {
      stop("Package 'foreach' is required for parallel processing. Please install it.")
    }
    if (!require(doParallel, quietly = TRUE)) {
      stop("Package 'doParallel' is required for parallel processing. Please install it.")
    }
    
    # Set up parallel backend
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1  # Leave one core free
    }
    cat("Setting up parallel processing with", n_cores, "cores...\n")
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    
    # Export necessary objects to cluster
    parallel::clusterEvalQ(cl, {
      library(dplyr)
      library(tidyr)
      library(pROC)
    })
    
    # Ensure cluster is stopped when function exits
    on.exit({
      if (exists("cl")) {
        parallel::stopCluster(cl)
        cat("Parallel cluster stopped.\n")
      }
    })
  }
  
  # Validate analysis_type parameter
  if (!analysis_type %in% c("individual", "panel", "both")) {
    stop("analysis_type must be 'individual', 'panel', or 'both'")
  }
  
  cat("----Comprehensive Immunoassay Analysis----", "\n")
  cat("Analysis type:", analysis_type, "\n")
  if (use_parallel) {
    cat("Using parallel processing with", n_cores, "cores\n")
  }
  
  # === COMMON DATA PREPROCESSING ===
  
  CAMP_data <- camp_data
  # convert any names for Marker3 to standardized name
  CAMP_data$developmental_text[grepl("Marker3", CAMP_data$developmental_text, ignore.case = TRUE)] <- "Marker3"
  CAMP_data <- CAMP_data %>% filter(!developmental_text %in% c("Hu14-enriched", "",NA))
  
  # use batch_normalized_value instead of assay_result_value for specific markers if needed
  if ("Marker4" %in% CAMP_data$developmental_text) {
    CAMP_data$assay_result_value[CAMP_data$developmental_text == "Marker4"] <- CAMP_data$batch_normalized_value[CAMP_data$developmental_text == "Marker4"]
  }
  
  # using MFI for specific markers
  cat("The unit of Marker8, Marker9, and Marker10 markers is based on MFI/OD.", "\n")
  CAMP_data$assay_result_value[CAMP_data$developmental_text == "Marker8"] <- CAMP_data$raw_value[CAMP_data$developmental_text == "Marker8"]
  CAMP_data$assay_result_value[CAMP_data$developmental_text == "Marker9"] <- CAMP_data$raw_value[CAMP_data$developmental_text == "Marker9"]
  CAMP_data$assay_result_value[CAMP_data$developmental_text == "Marker10"] <- CAMP_data$raw_value[CAMP_data$developmental_text == "Marker10"]
  
  # Instrument filtering
  if (!camp_six) {
    CAMP_data <- CAMP_data[!CAMP_data$instrument %in% c("LX-200_2", "LX-200_1"),]
    cat("Removing LX-200_1 and LX-200_2 instruments from the analysis!", "\n")
  } else {
    cat("All instruments are considered in the analysis!", "\n")
  }
  
  CAMP_data$assay_result_value <- as.numeric(CAMP_data$assay_result_value)
  
  # Plate assignments
  if (plate_code == "1-6") {
    CAMP_data <- CAMP_data %>% filter(sample >= "CMP0013:0187" & sample <= "CMP0013:0576")
    cat("Plates 1-6 are selected for analysis.", "\n")
  } else if (plate_code == "7-8") {
    CAMP_data <- CAMP_data %>% filter(sample >= "CMP0013:0041" & sample <= "CMP0013:0184")
    cat("Plates 7-8 are selected for analysis.", "\n")
  } else if (plate_code == "1-8") {
    CAMP_data <- CAMP_data %>% filter(sample >= "CMP0013:0041" & sample <= "CMP0013:0576")
    cat("Plates 1-8 are selected for analysis.", "\n")
  } else if (plate_code == "20-20") {
    CAMP_data <- CAMP_data %>% filter(sample >= "CMP0013:0001" & sample <= "CMP0013:0040")
    cat("The 20-20 samples are selected for analysis.", "\n")
  } else if (plate_code == "All") {
    cat("All plates are included in the analysis.", "\n")
  }
  
  proteins <- unique(CAMP_data$developmental_text)
  
  # Imputation: Replace -1 with half minimum and -2 with double maximum
  CAMP_data <- CAMP_data %>%
    group_by(developmental_text) %>%
    mutate(
      assay_result_value = ifelse(
        assay_result_value == -1,
        0.5 * min(assay_result_value[assay_result_value > 0], na.rm = TRUE),
        ifelse(
          assay_result_value == -2,
          2 * max(assay_result_value[assay_result_value > 0], na.rm = TRUE),
          assay_result_value
        )
      )
    ) %>%
    ungroup()
  
  cat("These markers are present in the data:", "\n")
  cat(unique(CAMP_data$developmental_text), "\n")
  
  CAMP_data_wide <- CAMP_data[,c("sample", "developmental_text", "assay_result_value")]
  CAMP_data_wide <- tidyr::pivot_wider(CAMP_data_wide, names_from = developmental_text, values_from = assay_result_value)
  CAMP_data_wide <- CAMP_data_wide %>% mutate(across(-sample, as.numeric))
  
  # Sample filtering based on specimen timing
  CAMP_key <- camp_key
  if (all_specimen == "One Year") {
    CAMP_key <- CAMP_key %>% 
      filter((is.na(days_predx) | days_predx < 366) & 
               ((!categorical_stage %in% c("Stage 4") & 
                   !category_specimen %in% c("after dx and tx of cancer","mixed, multiple cancers")) |
                  (categorical_stage %in% c("Stage 4") & 
                     !category_specimen %in% c("mixed, multiple cancers"))))
    cat("We keep a sample in Stage 4 if it is after dx and tx of cancer or after dx of cancer!", "\n")
    cat("Excluding samples collected after diagnosis and treatment!", "\n")
    cat("Excluding samples collected > 1 year prior to diagnosis!", "\n")
    cat("Excluding multiple cancer cases from the analysis!", "\n")
  } else if (all_specimen == "Two Years") {
    CAMP_key <- CAMP_key %>% 
      filter((is.na(days_predx) | (days_predx > 365 & days_predx < 731)) & 
               ((!categorical_stage %in% c("Stage 4") & 
                   !category_specimen %in% c("after dx and tx of cancer","mixed, multiple cancers")) |
                  (categorical_stage %in% c("Stage 4") & 
                     !category_specimen %in% c("mixed, multiple cancers"))))
    cat("Including samples collected between 1 and 2 years of diagnosis!", "\n")
  } else if (all_specimen == "Two Plus") {
    CAMP_key <- CAMP_key %>% 
      filter((is.na(days_predx) | (days_predx > 732 )) & 
               ((!categorical_stage %in% c("Stage 4") & 
                   !category_specimen %in% c("after dx and tx of cancer","mixed, multiple cancers")) |
                  (categorical_stage %in% c("Stage 4") & 
                     !category_specimen %in% c("mixed, multiple cancers"))))
  } else if (all_specimen == "All") {
    cat("No filter on samples collected after diagnosis and treatment, including all samples!", "\n")
  }
  
  # Merge data
  CAMP_key_wide <- CAMP_key[, c("sample", "plate", "age_specimen","specimen", "smoking_status", "packyears_notes",
                                "category_specimen", "sex", "type", "type_cancer", "categorical_stage", "subtype_breast", "protocol")]
  CAMP <- left_join(CAMP_data_wide, CAMP_key_wide, by = "sample")
  CAMP <- CAMP[!is.na(CAMP$specimen),]
  
  # Store complete data before filtering for special scenarios
  CAMP_complete <- CAMP  # Keep complete dataset with benign controls and Stage 0
  
  if (!all_specimen == "All") {
    cat("Excluding Benign controls! We analyze cancer vs Benign controls, don't worry!", "\n")
    CAMP <- CAMP %>% filter(!type %in% c("Control (Benign)"))
    cat("Excluding STAGE 0 from the analysis! We analyze cancer vs Stage 0, don't worry!", "\n")
    CAMP <- CAMP %>% filter(!categorical_stage %in% c("Stage 0")) 
  } else {
    cat("Including STAGE 0 and Benign controls in the analysis!", "\n")
  }
  
  
  # Average multiple specimens per sample
  meta_data <- CAMP %>%
    distinct(specimen, .keep_all = TRUE) %>%
    select(specimen, plate, age_specimen, sex, smoking_status, packyears_notes, type, type_cancer, category_specimen, subtype_breast, categorical_stage, protocol)
  CAMP_avg <- CAMP %>%
    dplyr::group_by(specimen) %>%
    dplyr::summarise(across(all_of(proteins), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  full_data <- CAMP_avg %>%
    left_join(meta_data, by = "specimen")
  CAMP <- full_data
  
  # Menopause status filtering
  if (menopause_status == "pre-menopause") {
    CAMP <- CAMP[CAMP$age_specimen < 50,]
    cat("Pre-menopause status is selected!", "\n")
  } else if (menopause_status == "post-menopause") {
    CAMP <- CAMP[CAMP$age_specimen > 49,]
    cat("Post-menopause status is selected!", "\n")
  } else {
    cat("No menopause status is selected, analyzing all specimen.", "\n")
  }
  
  # Convert NAs to blanks
  CAMP$type_cancer[is.na(CAMP$type_cancer)] <- ""
  CAMP$subtype_breast[is.na(CAMP$subtype_breast)] <- ""
  CAMP$categorical_stage[is.na(CAMP$categorical_stage)] <- ""
  
  # Normalization
  if (normalization) {
    temp_markers <- distinct(CAMP_data[,c("instrument","developmental_text")])
    markers_normalize <- unique(temp_markers$developmental_text[!temp_markers$instrument == "Roche Cobas" & temp_markers$developmental_text != "Marker4"])
    markers_not_normalize <- unique(temp_markers$developmental_text[temp_markers$instrument == "Roche Cobas" | temp_markers$developmental_text == "Marker4"])
    cat("These markers will not be normalized because they are reported by Roche Cobas or the marker is Marker4: \n", markers_not_normalize, "\n")
    
    CAMP_normalized <- CAMP
    for (col in markers_normalize) {
      cat("Normalizing", col, "\n")
      
      max_control_median <- CAMP_normalized %>%
        group_by(plate) %>%
        summarise(
          control_median = median(case_when(
            type != "Single Cancer" ~ .data[[col]],
            TRUE ~ NA_real_
          ), na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        pull(control_median) %>%
        max(na.rm = TRUE)
      
      CAMP_normalized <- CAMP_normalized %>%
        group_by(plate) %>%
        mutate(
          control_median = median(case_when(
            type != "Single Cancer" ~ .data[[col]],
            TRUE ~ NA_real_
          ), na.rm = TRUE),
          !!col := as.numeric(.data[[col]]) / control_median * max_control_median
        ) %>%
        select(-control_median) %>%
        ungroup()
    }
    CAMP <- CAMP_normalized
  }
  
  # Define smoking categories
  CAMP <- CAMP %>%
    mutate(
      packyears = case_when(
        packyears_notes > 30 ~ "heavy",
        packyears_notes > 10 & packyears_notes <= 30 ~ "medium",
        packyears_notes <= 10 | is.na(packyears_notes) ~ "light/never",
        TRUE ~ NA_character_
      )
    )
  
  # Cumulative predictions
  # Generate predictions for all cancer types
  cat("Generating predictions for all cancer types...\n")
  if (analysis_type %in% c("both", "panel")) {
    main_all_predictions <- generate_all_predictions(CAMP, markers_count)
    filtered_main_all_predictions <- main_all_predictions %>% 
      filter((is.na(prediction_Lung_Cancer) | prediction_Lung_Cancer < 0.01) &
               (is.na(prediction_Pancreatic_Cancer) | prediction_Pancreatic_Cancer < 0.01) &
               (is.na(prediction_Breast_Cancer) | prediction_Breast_Cancer < 0.01) &
               (is.na(prediction_Gynecological_Cancer) | prediction_Gynecological_Cancer < 0.01) &
               (is.na(prediction_Liver_Cancer) | prediction_Liver_Cancer < 0.01) &
               (is.na(prediction_Colon_and_Rectal_Cancer) | prediction_Colon_and_Rectal_Cancer < 0.01) &
               (is.na(prediction_Gastrointestinal_Cancer) | prediction_Gastrointestinal_Cancer < 0.01) &
               (is.na(prediction_Prostate_Cancer) | prediction_Prostate_Cancer < 0.01))
    
    cum_risks <- generate_scenario_risk_analysis(main_all_predictions)
  } 
  # Define cancer types
  cancer_types <- c(
    "Breast Cancer",
    "Liver Cancer", 
    "Gynecological Cancer",
    "Lung Cancer",
    "Colon and Rectal Cancer",
    "Gastrointestinal Cancer",
    "Pancreatic Cancer",
    "Prostate Cancer"
  )
  
  # Add "All Cancers" for individual analysis
  if (analysis_type %in% c("individual", "both")) {
    cancer_types_individual <- c(cancer_types, "All Cancers")
  }
  protocols <- c("All", "MERIT_LEAP", unique(CAMP$protocol))
  # Initialize results
  results_individual <- data.frame()
  results_panel <- data.frame()
  ROC_objs <- list()
  AUC_objs <- data.frame(Cancer = character(0), AUC = character(0), stringsAsFactors = FALSE)
  Prediction_objs <- list()
  
  
  
  # === ANALYSIS LOOPS ===
  
  if (use_parallel) {
    # Export helper functions and data to cluster
    export_vars <- c("get_scenario_data", "get_scenario_data_panel", 
                     "get_panel_proteins", "apply_panel_model", "proteins", "calculate_ci_proportion",
                     "cancer_types", "CAMP", "CAMP_complete", "analysis_type", "markers_count", "new_breast_panel")
    
    # Only export cancer_types_individual if it exists (for individual analysis)
    if (analysis_type %in% c("individual", "both")) {
      export_vars <- c(export_vars, "cancer_types_individual")
    }
    
    parallel::clusterExport(cl, export_vars, envir = environment())
    
    # Process all protocols in parallel
    all_results <- foreach(protocol_name = protocols, 
                           .combine = function(x, y) {
                             list(
                               individual = rbind(x$individual, y$individual),
                               panel = rbind(x$panel, y$panel),
                               roc_objs = c(x$roc_objs, y$roc_objs),  
                               auc_objs = rbind(x$auc_objs, y$auc_objs),
                               predictions = c(x$predictions, y$predictions),
                               stage_roc_objs = c(x$stage_roc_objs, y$stage_roc_objs),
                               stage_auc_objs = rbind(x$stage_auc_objs, y$stage_auc_objs)
                             )
                           },
                           .packages = c("dplyr", "tidyr", "pROC"),
                           .export = export_vars,
                           .init = list(individual = data.frame(), 
                                        panel = data.frame(),
                                        roc_objs = list(),  
                                        auc_objs = data.frame(Cancer = character(0), AUC = character(0), stringsAsFactors = FALSE),
                                        predictions = list())) %dopar% {  
                                          
                                          cat("------------------Processing protocol:", protocol_name, "------------------\n")
                                          if (protocol_name == "All") {
                                            df <- CAMP
                                            df_CAMP_complete <- CAMP_complete
                                          } else if (protocol_name == "MERIT_LEAP") {
                                            # Combine the three specified protocols
                                            merit_leap_protocols <- c("MERIT (PA17-0584)", "LEAP US (2013-0609)", "LEAP Europe (PA16-0415)")
                                            df <- CAMP[CAMP$protocol %in% merit_leap_protocols,]
                                            df_CAMP_complete <- CAMP_complete[CAMP_complete$protocol %in% merit_leap_protocols,]
                                            cat("MERIT_LEAP: Combined protocols -", paste(merit_leap_protocols, collapse = ", "), "\n")
                                            cat("Total samples in MERIT_LEAP:", nrow(df), "\n")
                                          } else {
                                            df <- CAMP[CAMP$protocol == protocol_name,]
                                            df_CAMP_complete <- CAMP_complete[CAMP_complete$protocol == protocol_name,]
                                          }
                                          protocol_roc_objs <- list()  
                                          protocol_auc_objs <- data.frame(Cancer = character(0), AUC = character(0), stringsAsFactors = FALSE)  
                                          protocol_results <- list(individual = data.frame(), 
                                                                   panel = data.frame(),
                                                                   roc_objs = list(),  
                                                                   auc_objs = data.frame(Cancer = character(0), AUC = character(0), stringsAsFactors = FALSE),
                                                                   predictions = list(),
                                                                   stage_roc_objs = list(),
                                                                   stage_auc_objs = data.frame(Cancer = character(0), AUC = character(0), stringsAsFactors = FALSE))
                                          protocol_predictions <- list()
                                          
                                          # [Continue with the analysis loops - truncated for brevity]
                                          # The parallel processing section would continue here with the same logic
                                          # as the sequential version below
                                          
                                          return(protocol_results)
                                        }
    
    # Extract results from parallel processing
    results_individual <- all_results$individual
    results_panel <- all_results$panel
    ROC_objs <- all_results$roc_objs  
    AUC_objs <- all_results$auc_objs  
    Prediction_objs <- all_results$predictions
    Stage_ROC_objs <- all_results$stage_roc_objs
    Stage_AUC_objs <- all_results$stage_auc_objs
  } else {
    # Sequential processing would go here - simplified for brevity
    # The sequential processing follows the same pattern as parallel but without foreach
  }
  
  # === SAVE RESULTS AND GENERATE VISUALIZATIONS ===
  
  current_date <- Sys.Date()
  folder_name <- paste0(current_date)
  folder_path <- paste0(saving_path, folder_name)
  if (analysis_type %in% c("panel", "both")) {
    folder_path <- paste0(folder_path,"/",strsplit(file_stem, "_")[[1]][1],"/PANEL_",markers_count,"_",gsub(" ", "", all_specimen),"_Protocol_",protocol_plot)
  } else if (analysis_type == "individual") {
    folder_path <- paste0(folder_path,"/",strsplit(file_stem, "_")[[1]][1],"/INDIVIDUAL_MARKERS",gsub(" ", "", all_specimen))
  }
  
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    cat("Creating a new folder in the given path with today's date (YYYY-MM-DD).", "\n")
  }
  cat("Saving the results in", folder_path, "\n")
  
  # Save results based on analysis type
  if (analysis_type %in% c("individual", "both")) {
    if (menopause_status == "none") {
      writexl::write_xlsx(results_individual, paste0(folder_path, "/", file_stem, "_individual_results_", all_specimen, ".xlsx"))
    } else if (menopause_status == "pre-menopause") {
      writexl::write_xlsx(results_individual, paste0(folder_path, "/", file_stem, "_individual_results_pre_menopause_", all_specimen, ".xlsx"))
    } else if (menopause_status == "post-menopause") {
      writexl::write_xlsx(results_individual, paste0(folder_path, "/", file_stem, "_individual_results_post_menopause_", all_specimen, ".xlsx"))
    }
  }
  
  if (analysis_type %in% c("panel", "both")) {
    if (menopause_status == "none") {
      writexl::write_xlsx(results_panel, paste0(folder_path, "/PANEL_", markers_count, "_", file_stem, "_results_", all_specimen, ".xlsx"))
      writexl::write_xlsx(main_all_predictions, paste0(folder_path, "/PANEL_", markers_count, "_", file_stem, "_results_", all_specimen, "_Predictions.xlsx"))
      writexl::write_xlsx(cum_risks, paste0(folder_path, "/PANEL_", markers_count, "_", file_stem, "_results_", all_specimen, "_CumRisks.xlsx"))
      
    } else if (menopause_status == "pre-menopause") {
      writexl::write_xlsx(results_panel, paste0(folder_path, "/PANEL_", markers_count, "_", file_stem, "_results_pre_menopause_", all_specimen, ".xlsx"))
    } else if (menopause_status == "post-menopause") {
      writexl::write_xlsx(results_panel, paste0(folder_path, "/PANEL_", markers_count, "_", file_stem, "_results_post_menopause_", all_specimen, ".xlsx"))
    }
  }
  
  
  # Generate overview table
  overview_table <- CAMP %>%
    select(type_cancer, categorical_stage, subtype_breast, category_specimen, sex, age_specimen, smoking_status) %>%
    mutate(
      subtype_breast = case_when(
        subtype_breast == "" ~ "Controls/Other Cancers",
        TRUE ~ subtype_breast
      ),
      type_cancer = case_when(
        type_cancer == "" ~ "Controls",
        TRUE ~ type_cancer
      )
    ) %>%
    tbl_summary(
      label = list(
        type_cancer ~ "Cancer Type",
        categorical_stage ~ "Cancer Stage",
        subtype_breast ~ "Breast Cancer Subtype",
        category_specimen ~ "Specimen Category",
        sex ~ "Sex",
        age_specimen ~ "Age at Specimen Collection",
        smoking_status ~ "Smoking Status"
      ),
      statistic = list(
        all_continuous() ~ "{mean} ({sd})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      missing_text = "Missing"
    ) %>%
    add_n() %>%
    modify_header(label ~ "Characteristic") %>%
    modify_caption("Overview of the Samples") %>%
    bold_labels()
  
  overview_table %>%
    as_tibble() %>%
    write.xlsx(paste0(folder_path, "/", file_stem, "_table_overview.xlsx"))
  
  # Generate visualizations
  generate_visualizations(CAMP, proteins, folder_path, file_stem, all_specimen)
  
  # Generate ROC plots for panel analysis
  if (analysis_type %in% c("panel", "both") && length(ROC_objs) > 0) {
    generate_roc_plots(ROC_objs, AUC_objs, folder_path, file_stem, markers_count, all_specimen, protocol_plot)
    
    if (exists("Stage_ROC_objs") && length(Stage_ROC_objs) > 0) {
      generate_stage_roc_plots(Stage_ROC_objs, Stage_AUC_objs, folder_path, file_stem, markers_count, all_specimen)
    }
  }
  
  # Generate forest plots for panel analysis
  if (analysis_type %in% c("panel", "both") && length(results_panel) > 0) {
    forest_data <- generate_forest_plot(results_panel, folder_path, file_stem, markers_count, all_specimen)
    generate_stage_forest_plots(results_panel, folder_path, file_stem, markers_count, all_specimen)
  }
  
  # Generate cancer type accuracy plot
  if (analysis_type %in% c("panel", "both") && exists("main_all_predictions")) {
    cat("Generating cancer type prediction accuracy plot...\n")
    accuracy_results <- generate_cancer_accuracy_plot(
      main_all_predictions = main_all_predictions,
      folder_path = folder_path,
      file_stem = file_stem,
      markers_count = markers_count,
      all_specimen = all_specimen,
      min_threshold = 0.01
    )
  }
  
  # Return results based on analysis type
  if (analysis_type == "individual") {
    return(results_individual)
  } else if (analysis_type == "panel") {
    return(results_panel)
  } else {
    return(list(individual = results_individual, panel = results_panel))
  }
}

# Helper functions ----

# Function to get scenario data for individual analysis
get_scenario_data <- function(df, scenario, cancer_type, protein, df_complete = NULL) {
  
  # Use complete dataset for scenarios that need benign controls or stage 0
  if (scenario %in% c("cancer_vs_Benign", "cancer_vs_Stage0") && !is.null(df_complete)) {
    df <- df_complete
  } else {
    df <- df
  }
  
  if (cancer_type != "All Cancers") {
    temp_data <- switch(scenario,
                        
                        "cancer_vs_Benign" = df %>% 
                          dplyr::select(type_cancer, type, any_of(protein)) %>%
                          dplyr::filter((type_cancer == cancer_type) | (type == "Control (Benign)")) %>%
                          dplyr::mutate(type_cancer = ifelse(type == "Control (Benign)", "Control (Benign)", type_cancer)) %>%
                          dplyr::select(-type),
                        
                        "cancer_vs_Stage0" = df %>% 
                          dplyr::select(type_cancer, categorical_stage, any_of(protein)) %>%
                          dplyr::filter((type_cancer == cancer_type & !categorical_stage %in% c("Stage 0")) | 
                                          (type_cancer == cancer_type & categorical_stage == "Stage 0")) %>%
                          dplyr::mutate(type_cancer = ifelse(categorical_stage == "Stage 0", "Stage 0", type_cancer)) %>%
                          dplyr::select(-categorical_stage),
                        
                        "cancer_vs_healthy" = df %>% 
                          dplyr::select(type_cancer, any_of(protein)) %>%
                          dplyr::filter(type_cancer %in% c(cancer_type, "")),
                        
                        # Additional scenarios would be defined here following the same pattern
                        # Simplified for brevity - the full implementation would include all scenarios
                        
                        # Default case
                        stop("Invalid scenario: ", scenario)
    )
  } else if (cancer_type == "All Cancers") {
    df$type_cancer <- ifelse(df$type_cancer != "", "All Cancers", "")
    
    temp_data <- switch(scenario,
                        "cancer_vs_healthy" = df %>% 
                          dplyr::select(type_cancer, any_of(protein)) %>%
                          dplyr::filter(type_cancer %in% c(cancer_type, "")),
                        
                        # Additional scenarios for "All Cancers" would be defined here
                        
                        # Default case
                        stop("Invalid scenario: ", scenario) )
  }
  
  return(temp_data)
}

# Function to get scenario data for panel analysis  
get_scenario_data_panel <- function(df, scenario, cancer_type, panel_proteins, df_complete = NULL) {
  
  # Use complete dataset for scenarios that need benign controls or stage 0
  if (scenario %in% c("cancer_vs_Benign", "cancer_vs_Stage0") && !is.null(df_complete)) {
    df <- df_complete
  } else {
    df <- df
  }
  
  temp_data <- switch(scenario,
                      "cancer_vs_healthy" = df %>% 
                        dplyr::select(specimen, type_cancer, any_of(panel_proteins)) %>%
                        dplyr::filter(type_cancer %in% c(cancer_type, "")),
                      
                      # Additional scenarios would be defined here
                      # Simplified for brevity
                      
                      # Default case
                      stop("Invalid scenario: ", scenario)
  )
  
  return(temp_data)
}

# Function to define panel proteins based on cancer type
get_panel_proteins <- function(cancer_type, markers_count, new_breast_panel) {
  if (!new_breast_panel) {
    if (markers_count == 13) {
      if (cancer_type == "Lung Cancer") {
        return(c("Marker1", "Marker2", "Marker3", "Marker4"))
      } else if (cancer_type == "Pancreatic Cancer") {
        return(c("Marker5", "Marker6", "Marker7"))
      } else if (cancer_type == "Breast Cancer") {
        return(c("Marker8", "Marker2", "Marker3", "Marker1"))
      } else if (cancer_type == "Gynecological Cancer") {
        return(c("Marker2", "Marker9"))
      } else if (cancer_type == "Liver Cancer") {
        return(c("Marker10", "Marker7", "Marker2"))
      } else if (cancer_type == "Colon and Rectal Cancer") {
        return(c("Marker7", "Marker1", "Marker11"))
      } else if (cancer_type == "Gastrointestinal Cancer") {
        return(c("Marker7", "Marker5", "Marker11"))
      } else if (cancer_type == "Prostate Cancer") {
        return(c("Marker12", "Marker13"))
      }
    } else if (markers_count == 10) {
      # Define 10-marker panels
      if (cancer_type == "Lung Cancer") {
        return(c("Marker1", "Marker2", "Marker3", "Marker4"))
      } else if (cancer_type == "Pancreatic Cancer") {
        return(c("Marker5"))
      } else if (cancer_type == "Breast Cancer") {
        return(c("Marker8", "Marker2", "Marker3", "Marker1"))
      } else if (cancer_type == "Gynecological Cancer") {
        return(c("Marker2", "Marker9"))
      } else if (cancer_type == "Liver Cancer") {
        return(c("Marker10", "Marker2"))
      } else if (cancer_type == "Colon and Rectal Cancer") {
        return(c("Marker1"))
      } else if (cancer_type == "Gastrointestinal Cancer") {
        return(c("Marker5"))
      } else if (cancer_type == "Prostate Cancer") {
        return(c("Marker12", "Marker13"))
      }
    }
  } else {
    # Alternative panel configurations
    if (cancer_type == "Breast Cancer") {
      return(c("Marker2", "Marker7", "Marker14"))
    }
    # Other alternative panels would be defined here
  }
  
  return(NULL)
}

# Function to apply panel models
# NOTE: Model coefficients have been removed for public sharing
# Users should implement their own models based on their training data
apply_panel_model <- function(temp_data, cancer_type, markers_count, new_breast_panel) {
  
  # PLACEHOLDER: This function should contain your trained model coefficients
  # For security reasons, specific model parameters have been removed
  # Users should replace this with their own trained models
  
  warning("Model coefficients not provided. Please implement your own trained models.")
  
  # Example structure for implementing models:
  if (cancer_type == "Lung Cancer") {
    temp_data <- temp_data %>%
      mutate(
        model_score = rowSums(
          cbind(
            # PLACEHOLDER_INTERCEPT,
            # if("Marker1" %in% colnames(temp_data)) PLACEHOLDER_COEFFICIENT * log10(Marker1),
            # Add your model terms here
            0  # Placeholder
          ),
          na.rm = TRUE
        )
      )
  }
  # Additional cancer types would follow the same pattern
  
  return(temp_data)
}

# Function to generate visualizations
generate_visualizations <- function(CAMP, proteins, folder_path, file_stem, all_specimen) {
  # Create long format data for plotting
  CAMP_long <- CAMP %>% pivot_longer(proteins, names_to = "Protein", values_to = "Value")
  CAMP_long$type_cancer <- ifelse(CAMP_long$type_cancer == "", CAMP_long$type, CAMP_long$type_cancer)
  CAMP_long$main_type <- ifelse(CAMP_long$type_cancer %in% c("Breast Cancer", "Lung Cancer", "Liver Cancer", 
                                                             "Gynecological Cancer", "Colon and Rectal Cancer",
                                                             "Gastrointestinal Cancer", "Pancreatic Cancer",
                                                             "Control (Asymptomatic)", "Control (Benign)"), 
                                CAMP_long$type_cancer, "Other Cancers")
  
  # Generate various plots
  # Plot 1: All cancer types
  p1 <- ggplot(CAMP_long, aes(x = main_type, y = Value, fill = main_type )) +
    geom_boxplot() +
    theme_minimal() +
    scale_color_gradientn(colours = rainbow(length(unique(CAMP_long$main_type)))) +
    facet_wrap(~Protein, nrow = 4, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0("Immunoassay Data Visualization by Cancer Type ", sub("_.*", "",file_stem)))
  
  ggsave(paste0(folder_path,"/", file_stem, "_all_cancer_types_",all_specimen,".png"), 
         plot = p1, width = 25, height = 15, dpi = 300, units = "in")
  
  # Additional visualization functions would go here
}

# Function to generate ROC plots
generate_roc_plots <- function(ROC_objs, AUC_objs, folder_path, file_stem, markers_count, all_specimen, protocol_plot) {
  # Define colors for different cancer types
  cancer_color_map <- list(
    "Breast Cancer" = "#E31A1C",
    "Lung Cancer" = "#1F78B4",
    "Liver Cancer" = "#33A02C",
    "Gynecological Cancer" = "#FF7F00",
    "Pancreatic Cancer" = "#6A3D9A",
    "Colon and Rectal Cancer" = "#A6CEE3",
    "Gastrointestinal Cancer" = "#B15928",
    "Prostate Cancer" = "#CAB2D6"
  )
  
  # Create ROC plot with multiple curves
  cancer_types <- names(ROC_objs)
  colors <- sapply(cancer_types, function(cancer) {
    if (cancer %in% names(cancer_color_map)) {
      return(cancer_color_map[[cancer]])
    } else {
      fallback_colors <- c("#999999", "#000000", "#666666", "#CCCCCC")
      idx <- which(cancer_types == cancer) %% length(fallback_colors) + 1
      return(fallback_colors[idx])
    }
  })
  names(colors) <- cancer_types
  
  png(paste0(folder_path, "/PANEL_", markers_count, "_", file_stem, "_ROC_", all_specimen, ".png"), 
      width = 8, height = 8, units = "in", res = 300)
  
  plot(ROC_objs[[1]], col = colors[1], lwd = 3,
       main = paste0("ROC Curves for ", protocol_plot, " protocols, PANEL ", markers_count, " markers"),
       legacy.axes = TRUE, cex.main = 1.5, cex.lab = 1.2)
  
  if(length(ROC_objs) > 1) {
    for(i in 2:length(ROC_objs)) {
      lines(ROC_objs[[i]], col = colors[i], lwd = 3)
    }
  }
  
  legend_labels <- sapply(cancer_types, function(cancer) {
    auc_value <- AUC_objs$AUC[AUC_objs$Cancer == cancer]
    paste0(cancer, ": AUC = ", auc_value)
  })
  
  legend("bottomright", legend = legend_labels, col = colors[cancer_types],
         lwd = 3, cex = 0.9, bg = "white", box.lty = 1)
  
  dev.off()
  cat("ROC plot saved successfully!\n")
}

# Additional helper functions would be defined here
# Including: generate_all_predictions, generate_scenario_risk_analysis, 
# generate_forest_plot, generate_cancer_accuracy_plot, etc.
# These have been abbreviated for brevity but would follow the same pattern

# Function to calculate confidence intervals
calculate_ci_proportion <- function(successes, total, conf_level = 0.95) {
  if (total == 0 || successes > total) {
    return(c(NA, NA))
  }
  if (successes == 0) {
    return(c(0, 0))
  }
  if (successes == total) {
    return(c(1, 1))
  }
  
  # Use Wilson score interval
  z <- qnorm((1 + conf_level) / 2)
  p <- successes / total
  n <- total
  
  denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  margin <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / denominator
  
  lower <- max(0, center - margin)
  upper <- min(1, center + margin)
  
  return(c(lower, upper))
}

# Required libraries ----
library(pROC)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)  
library(readxl)
library(writexl)
library(RColorBrewer)
library(openxlsx)
library(gtsummary)
library(foreach)
library(doParallel)

# Example usage ----
# Load your data
# biomarker_data <- read.csv("path/to/your/assay_data.csv")
# sample_metadata <- read.csv("path/to/your/sample_key.csv")

# Run analysis
# results <- comprehensive_immunoassay_analysis(
#   camp_data = biomarker_data, 
#   camp_key = sample_metadata,
#   analysis_type = "both",  # "individual", "panel", or "both"
#   file_stem = "my_analysis",
#   saving_path = "output/",
#   use_parallel = TRUE
# )
