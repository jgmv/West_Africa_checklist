### load custom functions
source("Code/functions.R")


### input data
data_with_genera <- parse_raw_checklist(remove_genus = F)
data <- parse_raw_checklist(remove_genus = T)
# output file: 'Output/modified_checklist.csv'


### number of records
nrow(data)
nrow(data_with_genera)


### number of species
length(unique(data$species))
length(unique(data_with_genera$species))


### number of publications
length(unique(data$literature_reference))
length(unique(data_with_genera$literature_reference))


### summary of studies' historical context
plot_study_motivation(data_with_genera)
# output file: 'Output/study_motivation_study.pdf' 

authors_data <- prepare_authors_data(data)
# output file: 'Output/first_author_activity.csv'
# file to be edited with author's origin

plot_author_publications(authors_data)
# output file: 'Output/author_publications.pdf' 

plot_records_per_time(data)
# output file: 'Output/records_x_year.pdf' 


### sumaries of records' ecology
species_ecology(data)
# output file: 'Output/sp_x_ecology.csv' 
# output file: 'Output/sp_x_ecology.pdf' 

lichen_substrata(data)
# output file: 'Output/lichen_substrata.pdf' 
# output file: 'Output/lichen_substrata.csv' 

### diversity analyses
plot_rank_abundance(data, n_labels = 8)
# output file: 'Output/rank_abundance.pdf' 

country <- diversity_per_country(data)
# output files: 'Output/country_summary.csv',
#  'Output/sp_x_ref.pdf' 
#  'Output/estimations_vs_known.pdf' 

plot_accumulation_curves(data)
# output files: 'Output/accumulation_curves_per_record.pdf',
#  'Output/accumulation_curve_per_publication.pdf', and
#  'Output/accumulation_curve_per_publication_WA.pdf'

sp_x_country <- species_per_country(data)
# output file: 'Output/sp_x_country.csv' 

top_rich_genera(data, 15)
# output file: 'Output/sp_x_genus.csv' 


### obtain proportion of species reported vs known
sp_x_taxon <- proportion_known_spp()
# output file: 'Output/sp_x_taxon.csv' 


### comparison of species with GBIF and IMI data
db_comparison(data)
# output file: 'Output/db_comparisons.pdf' 
# output file: 'Output/GBIF_spp.csv'
# output file: 'Output/IMI_spp.csv' 

db_comparison_country(data)
# output file: 'Output/db_comparisons.pdf' 
# output file: 'Output/GBIF_comparison_country.csv'
# output file: 'Output/IMI_comparison_country.csv'


### number of type specimens
type_specimens(data)
# output file: 'Output/types_x_country.csv' 


### number of species with incertae sedis taxa
(is_taxa <- incertae_sedis_taxa(data))

# number of records with incertae sedis order
length(grep("incertae", data$order))

# number of records with incertae sedis family
length(grep("incertae", data$family))

# number of Ascomycota species in dataset
asco_spp <- unique(data[data$division == "Ascomycota", "species"])
length(asco_spp)

# percentage of incertae sedis ascomycota species
100 * is_taxa[1] / length(asco_spp)


### plot map with colored countries
plot_map()
# output file: 'Output/West_Africa_map.pdf' 

plot_richness_map(data)
# output file: 'Output/West_Africa_map2.pdf' 


### species turnover with distance (Mantel correlogram)
mantel_correlogram(data, country)
# output file: 'Output/mantel_correlogram.pdf' 


### end
