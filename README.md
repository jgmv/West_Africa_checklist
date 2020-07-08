# West_Africa_checklist
Data and code for the data analyses used in [Piepenbring, Maciá-Vicente, Codjia, Glatthorn, Kirk, Meswaet, Minter, Olou, Reschke, Schmidt, Yorou (2020) Mapping mycological ignorance – checklists and diversity patterns of fungi known for West Africa. IMA Fungus 10:13](https://doi.org/10.1186/s43008-020-00034-y).

## Contents
### R_process.R
R script with commands to run all the analyses.

### Code/functions.R
Custom R functions called by the script `R_process.R`.

### Data/
Files with raw data necessary for the analyses.
* `checklist.csv`: checklist data, including all fungal records for West Africa and their associated metadata.
* `first_author_activity.csv`: data on the origin of first authors of publications on West African fungi (for Figure&nbsp;3c).
* `GBIF records.csv`: Fungal records for West Africa available in the [Global Biodiversity Information Facility (GBIF)](https://www.gbif.org/) database.
* `IMI_records.csv`: Fungal records for West Africa available in the former International Mycological Institute (IMI) fungal reference collection. 
* `known_spp_numbers.csv`: numbers of fungal species known worldwide for systematic groups of fungi and fungus-like organisms. Data from the [Catalogue of Life](http://www.catalogueoflife.org).
* `known_spp_numbers_selection.csv`: subset of `known_spp_numbers.csv` with a selection of systematic groups (for Table&nbsp;2).

## Usage
Run the script `R_process.R` in R, using as working directory the same folder where all data and code files and folders are stored. 
The following R packages need to be installed: `ggplot2`, `ggthemes`, `Hmisc`, `maps`, `RColorBrewer`, `rnaturalearth`, `sf`, `sp`, `vegan`.
Upon running the script, all results will be stored in a folder named `Output` that is automatically created. 
