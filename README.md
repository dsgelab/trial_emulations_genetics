# Integrating Genetic Data in Target Trial Emulations

This repository contains the supporting code for the manuscript **"Integrating genetic data in target trial emulations improves their design and informs the value of polygenic scores for prognostic and predictive enrichment"**, available on medRxiv with DOI: [10.1101/2024.11.05.24316763](https://doi.org/10.1101/2024.11.05.24316763).

## Overview

This repository provides code to accompany the manuscript, offering transparency on the analysis steps undertaken. The analyses are designed for use with the **FinnGen** cohort and leverage proprietary datasets specific to this study.

## Repository Structure

- `scripts/`: R scripts for data preprocessing, trial emulation, and analysis.
- `scripts_MR/`: R scripts related to Mendelian Randomization (MR) analyses.

## Prerequisites

The code requires access to FinnGen-specific data and is intended for use within the FinnGen research environment.

### Required Software and Packages

- **R (version 4.0 or higher)**
- R packages:
  - `dplyr`
  - `ggplot2`
  - `TwoSampleMR`
  - `data.table`
  - `survival`

## Usage

1. **Data Access:**

   - Approved access to the FinnGen dataset is required.

2. **Running Analyses:**

   - Use the provided scripts to replicate the analyses described in the manuscript.

3. **Mendelian Randomization Analyses:**

   - MR analysis scripts can be found in the `scripts_MR/` folder.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

We acknowledge the FinnGen project and all collaborators who contributed to this research.

## Citation

If you use any part of this repository in your work, please cite the following manuscript:

> German, J. et al. (2024). Integrating genetic data in target trial emulations improves their design and informs the value of polygenic scores for prognostic and predictive enrichment. medRxiv. DOI: [10.1101/2024.11.05.24316763](https://www.medrxiv.org/content/10.1101/2024.11.05.24316763v1.full-text).

## Contact

For questions or further information, please contact the corresponding authors listed in the manuscript.
