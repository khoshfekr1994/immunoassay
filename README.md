# Comprehensive Immunoassay Analysis Framework

A comprehensive R framework for analyzing multi-cancer biomarker panel data from immunoassay studies. This tool performs both individual marker analysis and multi-marker panel analysis for cancer detection and classification.

## Features

- **Individual Marker Analysis**: ROC analysis, sensitivity/specificity calculations, and performance metrics for individual biomarkers
- **Panel Analysis**: Multi-marker panel models for improved cancer detection accuracy
- **Comprehensive Scenarios**: Analysis across different patient populations, smoking status, cancer stages, and demographic groups
- **Parallel Processing**: Optional parallel computation for large-scale analyses
- **Visualization**: Automated generation of ROC curves, forest plots, box plots, and performance visualizations
- **Multiple Output Formats**: Excel reports, publication-ready plots, and statistical summaries

## Requirements

### R Dependencies
```r
# Core analysis packages
library(pROC)
library(dplyr)
library(tidyr)
library(ggplot2)

# Additional packages
library(ggpubr)
library(readxl)
library(writexl)
library(RColorBrewer)
library(openxlsx)
library(gtsummary)

# Parallel processing (optional)
library(foreach)
library(doParallel)
```

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/immunoassay-analysis.git
cd immunoassay-analysis
```

2. Install required R packages:
```r
# Install from CRAN
install.packages(c("pROC", "dplyr", "tidyr", "ggplot2", "ggpubr", 
                   "readxl", "writexl", "RColorBrewer", "openxlsx", 
                   "gtsummary", "foreach", "doParallel"))
```

## Usage

### Basic Usage

```r
# Source the main function
source("comprehensive_immunoassay_analysis.R")

# Load your data
biomarker_data <- read.csv("your_assay_data.csv")
sample_metadata <- read.csv("your_sample_key.csv")

# Run analysis
results <- comprehensive_immunoassay_analysis(
  camp_data = biomarker_data,
  camp_key = sample_metadata,
  analysis_type = "both",  # "individual", "panel", or "both"
  file_stem = "your_analysis_name",
  saving_path = "output/",
  use_parallel = TRUE
)
```

### Analysis Types

- **Individual Analysis**: Analyzes each biomarker separately
- **Panel Analysis**: Uses multi-marker models for improved performance
- **Both**: Runs both individual and panel analyses

### Key Parameters

- `analysis_type`: Type of analysis to perform ("individual", "panel", "both")
- `normalization`: Whether to apply plate-based normalization (TRUE/FALSE)
- `markers_count`: Number of markers to include in panel models
- `menopause_status`: Filter by menopause status ("pre-menopause", "post-menopause", "none")
- `use_parallel`: Enable parallel processing for faster computation

## Data Format Requirements

### Biomarker Data
Your assay data should include columns:
- `sample`: Sample identifier
- `developmental_text`: Marker name
- `assay_result_value`: Measured concentration/value
- `instrument`: Assay platform used

### Sample Metadata
Your sample key should include:
- `sample`: Sample identifier (matching biomarker data)
- `specimen`: Specimen identifier
- `type_cancer`: Cancer type classification
- `categorical_stage`: Cancer stage information
- `age_specimen`: Age at specimen collection
- `sex`: Patient sex
- `smoking_status`: Smoking history

## Output

The analysis generates:

### Statistical Reports
- Individual marker performance metrics (Excel)
- Panel model performance results (Excel)
- Comprehensive sample overview tables

### Visualizations
- ROC curves for individual markers and panels
- Forest plots showing performance metrics with confidence intervals
- Box plots by cancer type and demographic groups
- Cancer type prediction accuracy heatmaps

### Performance Metrics
- Area Under the Curve (AUC) with confidence intervals
- Sensitivity and specificity at various thresholds
- Positive and negative predictive values
- Risk-based classification performance

## Advanced Features

### Scenario Analysis
The framework supports multiple comparison scenarios:
- Cancer vs. healthy controls
- Cancer vs. other cancer types
- Stage-specific analyses
- Smoking status subgroups
- Age and sex stratifications

### Parallel Processing
For large datasets, enable parallel processing:
```r
results <- comprehensive_immunoassay_analysis(
  # ... other parameters
  use_parallel = TRUE,
  n_cores = 4  # specify number of cores
)
```

### Custom Panel Models
The framework supports custom multi-marker panel models with configurable:
- Marker combinations for different cancer types
- Model coefficients and transformations
- Risk thresholds and classification rules

## File Structure

```
project/
├── comprehensive_immunoassay_analysis.R  # Main analysis function
├── README.md                            # This file
├── example_data/                        # Example datasets
├── output/                             # Analysis results
└── documentation/                      # Additional documentation
```

## Performance Considerations

- **Memory**: Large datasets may require 8GB+ RAM
- **Processing Time**: Panel analysis can be computationally intensive
- **Parallel Processing**: Recommended for datasets with >1000 samples
- **Storage**: Output files can be substantial for comprehensive analyses

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Author

**Hamid Khoshfekr Rudsari**  
Department of Biostatistics  
MD Anderson Cancer Center

Contact:
khoshfekr1994@gmail.com
hkhoshfekr@mdanderson.org
