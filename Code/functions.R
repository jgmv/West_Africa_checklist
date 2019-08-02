### create folder for data output ----------------------------------------------
if (!dir.exists("Output")) dir.create("Output")


### parse raw checklist --------------------------------------------------------
parse_raw_checklist <- function(file = "Data/checklist.csv", remove_genus = T) {  
  # input data
  data <- read.csv(file, h = T, sep = ";")
  
  # create a column with genus data
  options(warn = -1)
  data$genus <- as.vector(do.call(rbind, strsplit(as.character(data$species),
    " "))[, 1])
  options(warn = 0)
  data <- data[, c(1:6, 31, 7:30)]
  
  # modify column names
  colnames(data) <- c("status", "kingdom", "division", "order", "family",
    "modifier", "genus", "species", "authors", "study_focus", "ecology",
    "substrate", "host_family", "literature_reference_with_comments",
    "literature_reference", "year", "Benin", "Burkina Faso",
    "Gambia", "Ghana", "Guinea", "Guinea-Bissau","Ivory Coast", "Liberia",
    "Mali", "Niger", "Nigeria", "Senegal", "Sierra Leone", "Togo", "comments")
  
  # remove synonyms (select desired option by uncommenting)
  #data <- droplevels(data[data$status == "valid", ])
  data <- droplevels(data[grep("^syn.", data$status, invert = T), ])
  
  # modify publication's year data
  data$year[data$year == "present publi"] <- "2019"
  data$year <- gsub("[^0-9\\.]", "", data$year)
  data$year[data$year == ""] <- NA
  
  # remove trailing spaces in columns
  for(i in ncol(data)) data[, i] <- trimws(data[, i], which = c("both"))
  
  # drop empty levels
  data <- droplevels(data)
  
  # modify country columns
  countries <- colnames(data)[17:30]
  for(i in countries) {
    x <- data[data[, i] == "x", -c(17:30)]
    x$country <- rep(i, nrow(x))
    assign(paste("country", i, sep = "_"), x)
    rm(x, i)
  }
  for(i in ls(pattern = "country_")) {
    if(i == ls(pattern = "country_")[1]) {
      x <- get(i)
    } else {
      x <- rbind(x, get(i))
    }
  }
  data_mod <- x
  
  # modify empty fields
  data_mod[data_mod$ecology == "", "ecology"] <- "unclear"
  data_mod[data_mod$host_family == "", "host_family"] <- NA
  levels(data_mod$substrate) <- c(levels(data_mod$substrate), "unknown")
  data_mod[data_mod$substrate == "", "substrate"] <- "unknown"
  data_mod <- droplevels(data_mod)
  write.table(data_mod, file = "Output/modified_checklist.csv", sep = ";",
    row.names = F)
  
  if(remove_genus) {
    data_mod <- droplevels(data_mod[grep(" sp.$", data_mod$species,
      invert = T), ])
  }    
  return(data_mod)  
}


### country colors -------------------------------------------------------------
col_country <- function(country = NULL, x = data) {  
  require(RColorBrewer)
  
  col_country <- brewer.pal(8, "Set2")
  col_country <- c(col_country, brewer.pal(6, "Set1"))
  names(col_country) <- unique(x$country)
  
  if(!is.null(country)) col_country <- col_country[country]
  
  return(col_country)  
}


### plot map -------------------------------------------------------------------
plot_map <- function(data = data) {  
  require(maps)
  
  # plot map with countries
  pdf("Output/West_Africa_map.pdf", w = 6, h = 6, pointsize = 14)
  map("world", xlim = c(-20, 54), ylim = c(-38, 40), boundary = F, interior = F,
    fill = T, col = "white", lty = 0, wrap = T)
  polygon(c(par("usr")[1], rep(par("usr")[2], 2), rep(par("usr")[1], 2)),
    c(rep(par("usr")[3], 2), rep(par("usr")[4], 2), par("usr")[3]),
    col = gray(0.85), border = NA)
  map("world", xlim = c(-35, 60), ylim = c(-60, 40),
    boundary = T, interior = F, fill = T, col = "white", border = "white",
    lty = 1, wrap = T, add = T)
  for(i in names(col_country())) {
    map("world", i, add = T, fill = T,
      col = col_country(i), boundary = T, interior = F, lty = 1,
      border = "white")
  }
  #for(i in names(col_country())) {
  #  text(get_centroid(i)[1], get_centroid(i)[2], i)
  #}
    
  #map.scale(ratio = F, cex = 0.8)
  legend("bottomright", legend = names(col_country()), ncol = 1,
    fill = col_country(), bty = "n", border = "transparent", cex = 1.3,
    inset = c(-0.1, 0.025), y.intersp = 0.75)
  dev.off()  
}


### proportion of reported vs known --------------------------------------------
proportion_known_spp <- function(file = "Data/knonw_spp_numbers.csv") {  
  # known fungi per country
  known_fungi <- rep(NA, length(unique(data$country)) + 1)
  names(known_fungi) <- c(unique(data$country), "total")
  known_fungi["total"] <- 4673
  known_fungi["Benin"] <- 415
  known_fungi["Burkina Faso"] <- 13
  known_fungi["Gambia"] <- 117
  known_fungi["Ghana"] <- 1530
  known_fungi["Guinea"] <- 786
  known_fungi["Guinea-Bissau"] <- 19
  known_fungi["Ivory Coast"] <- 1125
  known_fungi["Liberia"] <- 159
  known_fungi["Mali"] <- 117
  known_fungi["Niger"] <- 85
  known_fungi["Nigeria"] <- 1325
  known_fungi["Senegal"] <- 299
  known_fungi["Sierra Leone"] <- 1618
  known_fungi["Togo"] <- 543
  
  # read data
  x <- read.csv(file, h = T, sep = ";")
  
  # keep only first word
  x$taxon <- trimws(x$taxon, which = c("both"))
  options(warn = -1)
  x$taxon <- as.vector(do.call(rbind, strsplit(as.character(x$taxon),
                                               " "))[, 1])
  options(warn = 0)
  x[x$taxon == "Incertae", "taxon"] <- "Incertae sedis"
  x[x$taxon == "fungi", "taxon"] <- "Fungi and fungus-like"
  #x <- x[-which(x$taxon == "incertae"), ]
  #x <- x[-which(x$taxon == "fungi"), ]
  
  # use taxon names as rownames
  rownames(x) <- x$taxon
  
  # add columns for each country
  for(i in c(unique(data$country), "total")) {
    i_perc <- paste0(i, "_perc")
    x[, i] <- rep(NA, nrow(x))
    x[, i_perc] <- rep(NA, nrow(x))
    
    # subset dataset per country
    data_sub <- droplevels(data[data$country == i, ])
    if(i == "total") data_sub <- data
    
    # count species per taxon
    for(j in x$taxon) {      
      if(j %in% data$order) {
        n_spp <- as.vector(data_sub[data_sub$order == j, "species"])
        } else {
          n_spp <- as.vector(data_sub[data_sub$division == j, "species"])
        }
      if(j == "Fungi") n_spp <- data_sub[data_sub$kingdom == "Fungi", "species"]
      if(j == "Fungi and fungus-like") n_spp <-
          data_sub[, "species"]
          # data_sub[data_sub$kingdom != "Fungi", "species"] # only fungus-like     
      if(j == "Ascomycota_Incertae_sedis") {
        data_sub_sub <- data_sub[data_sub$division == "Ascomycota", ]
        data_sub_sub <- droplevels(data_sub_sub)        
        is_order <- grep("incertae", data_sub_sub$order, ignore.case = T)
        is_family <- grep("incertae", data_sub_sub$famyly, ignore.case = T)
        sel <- unique(c(is_order, is_family))        
        n_spp <- data_sub_sub[sel, "species"]
      }
      if(j == "Basidiomycota_Incertae_sedis") {
        data_sub_sub <- data_sub[data_sub$division == "Basidiomycota", ]
        data_sub_sub <- droplevels(data_sub_sub)        
        is_order <- grep("incertae", data_sub_sub$order, ignore.case = T)
        is_family <- grep("incertae", data_sub_sub$famyly, ignore.case = T)
        sel <- unique(c(is_order, is_family))        
        n_spp <- data_sub_sub[sel, "species"]
      }
      if(j == "Zygomycota_Incertae_sedis") {
        data_sub_sub <- data_sub[data_sub$division == "Zygomycota", ]
        data_sub_sub <- droplevels(data_sub_sub)        
        is_order <- grep("incertae", data_sub_sub$order, ignore.case = T)
        is_family <- grep("incertae", data_sub_sub$famyly, ignore.case = T)
        sel <- unique(c(is_order, is_family))        
        n_spp <- data_sub_sub[sel, "species"]
      }
      n_spp <- unique(n_spp)
      x[j, i] <- length(n_spp)      
    }    
    x[, i_perc] <- round((x[, i] * 100) / known_fungi[i], 2)
  }  
  # add last row with totals
  #row_total <- c("total", "", round(colSums(x[, 3:ncol(x)], na.rm = T), 0))
  #x <- rbind(x, row_total)  
  write.table(x, file = "Output/sp_x_taxon.csv", sep = ";", col.names = T,
    row.names = F)
  
  # selected taxa
  selected_taxa <- read.csv("Data/known_spp_numbers_selection.csv", sep = ";",
    h = F)
  selected_taxa <- as.vector(selected_taxa$V1)
  selected_taxa <- trimws(selected_taxa, which = c("both"))
  options(warn = -1)
  selected_taxa <- as.vector(do.call(rbind,
    strsplit(as.character(selected_taxa), " "))[, 1])
  options(warn = 0)
  selected_taxa[selected_taxa == "Incertae"] <- "Incertae sedis"
  selected_taxa <- selected_taxa[-which(selected_taxa == "incertae")]
  x_sel <- droplevels(x[rownames(x) %in% selected_taxa, ])
  x_not_sel <- droplevels(x[!(rownames(x) %in% selected_taxa), ])
  write.table(x, file = "Output/sp_x_taxon_selected_taxa.csv", sep = ";",
    col.names = T, row.names = F)  
  return(x)  
}


### plot study motivation ------------------------------------------------------
plot_study_motivation <- function(data) {  
  require(Hmisc)
  
  study_motiv <- table(data$literature_reference, data$study_focus)
  study_motiv[study_motiv > 0] <- 1
  study_motiv[rowSums(study_motiv) > 1, ]

  x <- rep(NA, nrow(study_motiv))
  for(i in 1:length(study_motiv[, "p"])) if(study_motiv[i, "p"]) x[i] <- "p"
  for(i in 1:length(study_motiv[, "m"])) if(study_motiv[i, "m"]) x[i] <- "m"
  for(i in 1:length(study_motiv[, "l"])) if(study_motiv[i, "l"]) x[i] <- "l"
  for(i in 1:length(study_motiv[, "other"])) if(study_motiv[i, "other"]) {
    x[i] <- "other"
  }
  study_motiv <- as.data.frame(cbind(rownames(study_motiv), x))
  colnames(study_motiv) <- c("literature_reference", "study_focus")
  study_motiv$year <- substr(as.vector(study_motiv$literature_reference),
    nchar(as.vector(study_motiv$literature_reference)) - 4,
    nchar(as.vector(study_motiv$literature_reference)))
  study_motiv$year <- gsub("[a-z]", "", study_motiv$year)
  study_motiv$year <- trimws(study_motiv$year)
  study_motiv$year <- as.numeric(study_motiv$year)
  
  pdf("Output/study_motivation_study.pdf", w = 5, h = 4, pointsize = 14)
  par(xpd = F, mar = c(4, 4, 1, 1), tck = -0.025, mgp = c(3, 0.6, 0.6))
  years <- as.numeric(data$year)
  plot(1, 1, type = "n", ylim = c(min(years), max(years)), xlim = c(0.5, 4.5),
    axes = F, ylab = "Year", xlab = NA)
  grid(lty = 1, col = gray(0.9))
  for(i in study_motiv[study_motiv$study_focus == "m", "year"]) {
    j <- jitter(1, amount = 0.3)
    points(j, i, pch = 16, col = alpha("#88AA00", 0.5), cex = 1.2)
    rm(i, j)
  }
  for(i in study_motiv[study_motiv$study_focus == "p", "year"]) {
    j <- jitter(2, amount = 0.3)
    points(j, i, pch = 16, col = alpha("#88AA00", 0.5), cex = 1.2)
    rm(i, j)
  }
  for(i in study_motiv[study_motiv$study_focus == "l", "year"]) {
    j <- jitter(3, amount = 0.3)
    points(j, i, pch = 16, col = alpha("#88AA00", 0.5), cex = 1.2)
    rm(i, j)
  }
  for(i in study_motiv[study_motiv$study_focus == "other", "year"]) {
    j <- jitter(4, amount = 0.3)
    points(j, i, pch = 16, col = alpha("#88AA00", 0.5), cex = 1.2)
    rm(i, j)
  }
  axis(1, at = 1:4, lwd = 0, lwd.tick = 1,
    label = c("Mycological", "Phytopathological", "Lichenology",
    "other/unknown"))
  axis(1, at = c(par("usr")[1], par("usr")[2]), lwd = 1, lwd.tick = 0,
    label = NA)
  axis(2, las = 2, lwd = 0, lwd.tick = 1)
  axis(2, at = c(par("usr")[3], par("usr")[4]), lwd = 1, lwd.tick = 0,
    label = NA)
  legend_marks <- c(1, 5, 10, 15)
  text(1, 1825, sum(study_motiv$study_focus == "m"))
  text(2, 1825, sum(study_motiv$study_focus == "p"))
  text(3, 1825, sum(study_motiv$study_focus == "l"))
  text(4, 1825, sum(study_motiv$study_focus == "other"))
  dev.off()  
}


### plot rank abundance --------------------------------------------------------
plot_rank_abundance <- function(data, n_labels = 5) {  
  require(Hmisc)
  
  tab <- sort(table(data$species) , decreasing = T)
  x <- as.vector(tab)
  names(x) <- names(tab)
  
  # color by ecology
  pt_col <- gray(0.2)
  if(exists("data_ecology")) { # generated by 'secies_ecology()'    
    sp_x_ecology <- sort(table(data_ecology$ecology), decreasing = T)
    sp_x_ecology <- sp_x_ecology[c(1, 4, 2:3, 5:18)]
    unk <- which(names(sp_x_ecology) == "unclear")
    sp_x_ecology <- sp_x_ecology[c(1:(unk - 1),
                                   (unk + 1):length(sp_x_ecology), unk)] 
    col_ecology <- color(length(sp_x_ecology))
    names(col_ecology) <- names(sp_x_ecology)
    
    rownames(data_ecology) <- data_ecology$species
    ecol <- as.character(data_ecology[names(x), "ecology"])
    pt_col <- col_ecology[ecol]    
  }
  
  pdf("Output/rank_abundance.pdf", w = 7, h = 3, pointsize = 14)
  par(xpd = F, mar = c(4, 4, 1, 1), tck = -0.025, mgp = c(3, 0.6, 0.6))
  plot(x, log = "x", axes = F, ylab = "Number of records",
    xlab = "Species rank", type = "n", cex = 1.5)
  grid(lty = 1, col = gray(0.9))
  points(x, pch = 16, col = alpha(pt_col, 0.5), cex = 1.5)
  axis(1, at = c(0, 1, 10, 100, 1000, 5000), las = 1)
  axis(2, at = c(par("usr")[3], 200), lwd.ticks = 0, labels = NA)
  axis(2, at = c(0, 50, 100, 150, 200), lwd = 0, lwd.ticks = 1, las = 1)
  
  if(n_labels > 0) {
    for(i in 1:n_labels) {
      text(i, x[i], names(x)[i], adj = -0.1, font = 3)
    }
  }
  dev.off()
  
  # print number of records per category
  for(i in 1:5) {
    message(paste0("Species with ", i, " counts = ", sum(x == i),
      " (", round((sum(x == i) / length(x)) * 100, 1), " %)"))
  }  
  message(paste0("Species with > 5 counts = ", sum(x > 5),
    " (", round((sum(x > 5) / length(x)) * 100, 1), " %)"))
  
  write.table(x, file = "Output/records_x_species.csv", sep = ";",
    col.names = T, row.names = T)  
}


### calculate diversity per country --------------------------------------------
diversity_per_country <- function(data) {  
  require(vegan)
  require(maps)
  
  # sp_x_country
  sp_x_country <- table(data$country, data$species)
  sp_x_country[sp_x_country > 0] <- 1
  
  # species per publication
  sp_x_ref <- table(data$literature_reference, data$species)
  sp_x_ref[sp_x_ref > 0] <- 1
  write.table(sp_x_ref, file = "Output/sp_x_ref.csv", col.names = NA, sep = ";")
  
  # references per country
  x <- table(data$country, data$literature_reference)
  x[x > 0] <- 1
  ref_x_country <- rowSums(x)
  ref_x_country <- as.data.frame(ref_x_country)
  colnames(ref_x_country) <- "references"
  ref_x_country <<- ref_x_country
  
  # species per country and publication (list)
  sp_x_ref_x_country <- list()
  for(i in unique(data$country)) {
    x <- table(data[data$country == i,
                    "literature_reference"], data[data$country == i, "species"])
    if(sum(x) == 1) {
      sp  <- colnames(x)[colSums(x) > 0]
      ref <- rownames(x)[rowSums(x) > 0]
      x <- x[x > 0]
      names(x) <- sp
    } else {
      x <- x[rowSums(x) > 0, ]
      x <- x[, colSums(x) > 0]
      x[x > 0] <- 1
    }
    x <- x[, grep(" sp.$", colnames(x), invert = T)] # remove genus-level records
    sp_x_ref_x_country[[i]] <- x
    rm(x, i)
  }
  
  # records per country
  country <- table(data$country)
  country <- as.data.frame(country)
  colnames(country) <- c("country", "records")
  rownames(country) <- country$country
  country$n_species <- rowSums(sp_x_country)
  country$n_publications <- ref_x_country$references
  country <- country[, -which(colnames(country) == "country")]
  
  # richness and richness estimators
  for(i in 1:length(names(sp_x_ref_x_country))) {
    x <- specpool(sp_x_ref_x_country[[i]])[, c(2, 3, 7, 8)]
    if(i == 1) {
      est <- x
    } else {
      est <- rbind(est, x)
    }
    rm(x, i)
  }
  country <- cbind(country, est)
  
  # country centroids and areas
  country$x <- rep(NA, nrow(country))
  country$y <- rep(NA, nrow(country))
  for(i in rownames(country)) country[i, c("x", "y")] <- get_centroid(i); rm(i)
  m <- map("world", rownames(country), proj = "bonne", param = 45, plot = F,
           fill = T)
  areas_km <- rep(NA, nrow(country))
  names(areas_km) <- rownames(country)
  for(i in names(areas_km)) areas_km[i] <- area.map(m, i) / 0.3861
  country$sq_km <- areas_km
  
  # number of vascular plants known per country
  country$vascular_plants <- rep(NA, length(unique(data$country)))
  country["Benin", "vascular_plants"] <- 2807
  country["Burkina Faso", "vascular_plants"] <- 2080
  country["Gambia", "vascular_plants"] <- 1760
  country["Ghana", "vascular_plants"] <- 2971
  country["Guinea", "vascular_plants"] <- 2923
  country["Guinea-Bissau", "vascular_plants"] <- 1507
  country["Ivory Coast", "vascular_plants"] <- 3853
  country["Liberia", "vascular_plants"] <- 2403
  country["Mali", "vascular_plants"] <- 1739
  country["Niger", "vascular_plants"] <- 1218
  country["Nigeria", "vascular_plants"] <- 3378
  country["Senegal", "vascular_plants"] <- 2300
  country["Sierra Leone", "vascular_plants"] <- 1883
  country["Togo", "vascular_plants"] <- 3134
  
  # total values
  est_total <- specpool(sp_x_ref)[, c(2, 3, 7, 8)]
  total <- c(sum(country$records), length(unique(data$species)),
    length(unique(data$literature_reference)), as.numeric(est_total),
    NA, NA, sum(country$sq_km), sum(country$vascular_plants, na.rm = T))
  country <- rbind(total, country)
  rownames(country)[1] <- "West Africa"
  country["West Africa", "vascular_plants"] <- 7072
  
  
  # Hawksworth index-based estimation
  country$hawksworth_estimate <- country$vascular_plants * 6
  
  # percentage known vs estimated richness
  country$known_richness_perc <-
    (country$n_species / country$hawksworth_estimate) * 100
  
  write.table(country, file = "Output/country_summary.csv", col.names = NA,
    row.names = T, sep = ";")
  
  # plot barplot
  #bp <- country[, c("n_species", "chao", "hawksworth_estimate")]
  bp <- country[, c("n_species", "hawksworth_estimate")]
  #bp$chao <- bp$chao - bp$n_species
  #bp$hawksworth_estimate <- bp$hawksworth_estimate - bp$chao - bp$n_species
  bp$hawksworth_estimate <- bp$hawksworth_estimate - bp$n_species
  bp$n_species <- (bp$n_species / country$hawksworth_estimate) * 100
  #bp$chao <- (bp$chao / country$hawksworth_estimate) * 100
  bp$hawksworth_estimate <-
    (bp$hawksworth_estimate / country$hawksworth_estimate) * 100
  bp <- as.matrix(bp)
  
  bcol <- c("#2c5aa0", "#d7e3f4")
  pdf("Output/estimations_vs_known.pdf", w = 8, h = 3.5, pointsize = 14)
  par(mar = c(8, 4, 1, 7), xpd = T, tck = -0.025, mgp = c(3, 0.6, 0.6))
  x <- barplot(t(bp), border = NA, ylab = NA, width = c(2, rep(1, 14)),
    space = c(0.2, 0.4, rep(0.2, 13)), col = bcol, axes = F,
    names.arg = rep("", nrow(bp)))
  axis(2, pos = 0, las = 2, padj = 0.25)
  mtext("Proportion (%)", side = 2, line = 2)
  legend("topright", legend = c("Reported", "Hawksworth's"),
    fill = bcol, inset = c(-0.3, 0), bty = "n", border = F)
  text(x = x, y = -10, rownames(bp), xpd = TRUE, srt = 45, adj = 1)
  text(x, bp[, "n_species"] + 10, round(bp[, "n_species"], 1))
  dev.off()  
  return(country)  
}


### plot accumulation curves ---------------------------------------------------
plot_accumulation_curves <- function(data) {  
  require(vegan)
  
  # records per species and country
  records_sp_country <- table(data$country, data$species)
  
  # species per country and publication (list)
  sp_x_ref_x_country <- list()
  for(i in unique(data$country)) {
    x <- table(data[data$country == i, "literature_reference"],
      data[data$country == i, "species"])
    if(sum(x) == 1) {
      sp  <- colnames(x)[colSums(x) > 0]
      ref <- rownames(x)[rowSums(x) > 0]
      x <- x[x > 0]
      names(x) <- sp
    } else {
      x <- x[rowSums(x) > 0, ]
      x <- x[, colSums(x) > 0]
      x[x > 0] <- 1
    }
    sp_x_ref_x_country[[i]] <- x
    rm(x, i)
  }
  
  # plot accumulation curves per record
  records_sp_country_curves <- rarecurve(records_sp_country)
  dev.off()
  pdf("Output/accumulation_curves_per_record.pdf", w = 3.5, h = 3,
    pointsize = 14)
  par(xpd = F, mar = c(4, 4, 1, 1), tck = -0.025, mgp = c(3, 0.6, 0.6))
  plot(0, 0, type = "n", xlim = c(0, 4500), ylim = c(0, 1800), xlab = "Records",
    ylab = "Species", axes = F)
  grid(lty = 1, col = gray(0.9))
  for(i in 1:nrow(records_sp_country)) {
    lines(records_sp_country_curves[[i]], col = col_country(i))
    #text(nrow(sp_x_ref_x_country[[i]]) + 1,
    #ncol(sp_x_ref_x_country[[i]]), i, adj = 0)
    points(length(records_sp_country_curves[[i]]),
      records_sp_country_curves[[i]][length(records_sp_country_curves[[i]])],
      col = col_country(i), pch = 16)
  }
  axis(1, lwd = 0, lwd.tick = 1)
  axis(1, at = c(par("usr")[1], par("usr")[2]), lwd = 1, lwd.tick = 0,
    label = NA)
  axis(2, las = 2, lwd = 0, lwd.tick = 1)
  axis(2, at = c(par("usr")[3], par("usr")[4]), lwd = 1, lwd.tick = 0,
    label = NA)
  #legend("topleft", legend = unique(data$country),
  #  col = col_country, lty = 1)
  dev.off()
  
  # per publication
  pdf("Output/accumulation_curve_per_publication.pdf", w = 3.5, h = 3,
    pointsize = 14)
  par(xpd = F, mar = c(4, 4, 1, 1), las = 1, tck = -0.025,
    mgp = c(3, 0.6, 0.6))
  plot(0, 0, type = "n", xlim = c(0, max(ref_x_country)), ylim = c(0, 1800),
    xlab = "Publications", ylab = "Species", axes = F)
  grid(lty = 1, col = gray(0.9))
  for(i in unique(data$country)) {
    if(sum(sp_x_ref_x_country[[i]]) > 1) {
      plot(specaccum(sp_x_ref_x_country[[i]]), add = 2, ci = 0,
        col = col_country(i))
      points(nrow(sp_x_ref_x_country[[i]]), ncol(sp_x_ref_x_country[[i]]),
        col = col_country(i), pch = 16)
      #text(nrow(sp_x_ref_x_country[[i]]) + 1,
      #  ncol(sp_x_ref_x_country[[i]]), i, adj = 0)
    }
  }
  axis(1, lwd = 0, lwd.tick = 1)
  axis(1, at = c(par("usr")[1], par("usr")[2]), lwd = 1, lwd.tick = 0,
    label = NA)
  axis(2, las = 2, lwd = 0, lwd.tick = 1)
  axis(2, at = c(par("usr")[3], par("usr")[4]), lwd = 1, lwd.tick = 0,
    label = NA)
  #legend("topleft", legend = unique(data$country),
  #  col = col_country, lty = 1)
  dev.off()
  
  # total accumulation per reference
  sp_x_ref <- table(data$literature_reference, data$species)
  sp_x_ref[sp_x_ref > 0] <- 1
  ref_accum <- specaccum(sp_x_ref)
  ref_est <- poolaccum(sp_x_ref)
  
  lcol <- "#2c5aa0"
  pdf("Output/accumulation_curve_per_publication_WA.pdf", w = 3.5, h = 3,
      pointsize = 14)
  par(xpd = F, mar = c(4, 4, 1, 1), las = 1, tck = -0.025, mgp = c(3, 0.6, 0.6))
  plot(0, 0, type = "n", xlim = c(0, max(ref_accum$sites)),
    ylim = c(0, max(summary(ref_est)$chao[,"Chao"])), xlab = "Publications",
    ylab = "Species", axes = F)
  grid(lty = 1, col = gray(0.9))
  plot(ref_accum, add = 2, ci = 0, col = lcol)
  points(max(ref_accum$sites), max(ref_accum$richness), col = lcol, pch = 16)
  lines(summary(ref_est)$boot[,"Bootstrap"], col = lcol, lty = 2)
  points(max(ref_accum$sites), max(summary(ref_est)$boot[,"Bootstrap"]),
    col = lcol, pch = 16)
  lines(summary(ref_est)$chao[,"Chao"], col = lcol, lty = 6)
  points(max(ref_accum$sites), max(summary(ref_est)$chao[,"Chao"]),
    col = lcol, pch = 16)
  axis(1, lwd = 0, lwd.tick = 1, tck = -0.05)
  axis(1, at = c(par("usr")[1], par("usr")[2]), lwd = 1, lwd.tick = 0,
    label = NA)
  axis(2, las = 2, lwd = 0, lwd.tick = 1, tck = -0.05)
  axis(2, at = c(par("usr")[3], par("usr")[4]), lwd = 1, lwd.tick = 0,
    label = NA)
  legend("bottomright", legend = c("Reported", "Bootstrap", "Chao"),
    col = lcol, lty = c(1, 2, 6), bty = "n")
  dev.off()  
}


### matrix of species records per country --------------------------------------
species_per_country <- function(data) {  
  sp_x_country <- table(data$country, data$species)
  sp_x_country <- t(sp_x_country)
  sp_x_country[sp_x_country > 1] <- 1
  sp_x_country <- cbind(sp_x_country, rowSums(sp_x_country))
  colnames(sp_x_country)[ncol(sp_x_country)] <- "n_countries"
  
  #sp_x_country <- cbind(sp_x_country,
  #                      100 * (rowSums(sp_x_country) / ncol(sp_x_country)))
  #colnames(sp_x_country)[ncol(sp_x_country)] <- "perc_countries"
  #sp_x_country <-  sp_x_country[with(sp_x_country, order("perc_countries")), ]
  write.table(sp_x_country, file = "Output/sp_x_country.csv", col.names = NA,
    row.names = T, sep = ";")
  
  hist(sp_x_country[, "n_countries"])
  # print number of records per category
  for(i in 1:length(unique(data$country))) {
    message(paste0("Species in ", i, " countries = ",
      sum(sp_x_country[, "n_countries"] == i),
      " (", round(100 * sum(sp_x_country[, "n_countries"] == i) /
      nrow(sp_x_country), 1), " %)"))
  }  
  return(sp_x_country)  
}


### prepare authors data -------------------------------------------------------
prepare_authors_data <- function(data,
  file = "Data/first_author_activity.csv") {  
  authors <- data[, c("literature_reference", "year")]
  authors$literature_reference <- data$literature_reference
  strsplit(as.character(authors$literature_reference), " et")[[1]]
  x <- sapply(strsplit(as.character(authors$literature_reference), " et"),
    '[', 1)
  x <- sapply(strsplit(x, " &"), '[', 1)
  x <- gsub(' [0-9]+.*', '', x)
  authors$first_author <- x
  authors <- unique(authors)  
  authors_tab <- table(authors$first_author, authors$year)
  authors_data <- matrix(NA, ncol = 3,
    nrow = length(unique(authors$first_author)),
    dimnames = list(unique(authors$first_author), c("n", "first", "last")))
  for(i in rownames(authors_data)) {
    x <- authors_tab[i, ]
    n <- sum(x)
    x <- x[x > 0]
    f <- min(as.numeric(names(x)))
    l <- max(as.numeric(names(x)))
    authors_data[i, ] <- c(n, f, l)
    rm(x, n, f, l)
  }
  authors_data <- as.data.frame(authors_data)
  authors_data_extracted <- authors_data
  
  # input manualy inserted data
  authors_data <- read.csv(file, h = T, sep = ";", row.names = 1)  
  authors_data_extracted$origin <- rep(NA, nrow(authors_data_extracted))
  for(i in rownames(authors_data_extracted)) {
    authors_data_extracted[i, "origin"] <- as.character(authors_data[i,
      "origin"])
  }
  write.table(authors_data_extracted, file = "Output/first_author_activity.csv",
    col.names = NA, row.names = T, sep = ";")  
  return(authors_data)  
}


### plot author publications ---------------------------------------------------
plot_author_publications <- function(authors_data = authors_data) {  
  pdf("Output/author_publications.pdf", w = 5, h = 4, pointsize = 14)
  par(xpd = F, mar = c(4, 4, 1, 1), tck = -0.025, mgp = c(3, 0.6, 0.6))
  plot(1, 1, type = "n",
    ylim = c(min(authors_data$first), max(authors_data$last)), 
    xlim = c(0.5, 3.5), axes = F, ylab = "Year", xlab = NA)
  #grid(nx= NULL, ny = NULL, lty = 1, col = gray(0.9))
  abline(v = 1:3, h = c(1800, 1850, 1900, 1950, 2000), lty = 1, col = gray(0.9))
  for(i in rownames(authors_data)[authors_data$origin == "African"]) {
    j <- jitter(1, amount = 0.3)
    lines(c(j, j), c(authors_data[i, "first"], authors_data[i, "last"]),
      lwd = 2 + authors_data[i, "n"], col = alpha("#2C89A0", 0.5))
    rm(i, j)
  }
  for(i in rownames(authors_data)[authors_data$origin == "European"]) {
    j <- jitter(2, amount = 0.3)
    lines(c(j, j), c(authors_data[i, "first"], authors_data[i, "last"]),
      lwd = 2 + authors_data[i, "n"], col = alpha("#2C89A0", 0.5))
    rm(i, j)
  }
  for(i in rownames(authors_data)[authors_data$origin == "other"]) {
    j <- jitter(3, amount = 0.3)
    lines(c(j, j), c(authors_data[i, "first"], authors_data[i, "last"]),
      lwd = 2 + authors_data[i, "n"], col = alpha("#2C89A0", 0.5))
    rm(i, j)
  }
  axis(1, at = 1:3, lwd = 0, lwd.tick = 1,
    label = c("African", "European", "other/unknown"))
  axis(1, at = c(par("usr")[1], par("usr")[2]), lwd = 1, lwd.tick = 0,
    label = NA)
  axis(2, las = 2, lwd = 0, lwd.tick = 1)
  axis(2, at = c(par("usr")[3], par("usr")[4]), lwd = 1, lwd.tick = 0,
    label = NA)
  legend_marks <- c(1, 5, 10, 15)
  legend("bottomright", legend = legend_marks, lty = 1, lwd = 2 + legend_marks,
    col = alpha("#2C89A0", 0.5), bty = "n", title = "Publications")
  rm(legend_marks)
  dev.off()  
}


### plot records and species per year ------------------------------------------
plot_records_per_time <- function(data) {  
  # species per year
  sp_x_year <- table(data$year, data$species)
  sp_x_year[sp_x_year > 0] <- 1
  write.table(sp_x_year, file = "Output/sp_x_year.csv", col.names = NA,
    sep = ";")
  
  # plot cummulative species per year
  records_x_year <- table(data$year, data$species)
  sp_year <- as.data.frame(cbind(as.numeric(rownames(records_x_year)),
    rowSums(records_x_year)))
  #sp_year <- as.data.frame(cbind(as.numeric(rownames(sp_x_year)),
  #  rowSums(sp_x_year)))
  colnames(sp_year) <- c("year", "species")
  sp_year$cummulative <- cumsum(sp_year$species)
  new_sp <- function(data) {
    x <- cumsum(data)
    x[x > 1] <- 1
    return(x)
  }
  sp_year$new_sp <- rowSums(apply(sp_x_year, 2, function(x) new_sp(x)))
  rm(new_sp)
  
  pdf("Output/records_x_year.pdf", w = 6, h = 3, pointsize = 14)
  par(mar = c(4, 4, 1, 1), tck = -0.025, mgp = c(3, 0.6, 0.6))
  plot(sp_year$cummulative ~ sp_year$year, type = "s", ylab = "n",
    xlab = "Year", axes = F, col = NA)
  grid(lty = 1, col = gray(0.9))
  polygon(c(sp_year$year, sp_year$year[length(sp_year$year)]),
    c(sp_year$cummulative, 0), col = gray(0.75), border = NA)
  polygon(c(sp_year$year, sp_year$year[length(sp_year$year)]),
    c(sp_year$new_sp, 0), col = "#d35f5f", border = NA)
  points(sp_year$species ~ sp_year$year, type = "h", col = gray(0.25))
  axis(1, lwd = 0, lwd.tick = 1)
  axis(1, at = c(par("usr")[1], par("usr")[2]), lwd = 1, lwd.tick = 0,
    label = NA)
  axis(2, las = 2, lwd = 0, lwd.tick = 1)
  axis(2, at = c(par("usr")[3], par("usr")[4]), lwd = 1, lwd.tick = 0,
    label = NA)
  legend("topleft", c("Records", "Species"), fill = c(gray(0.75), "#d35f5f"),
    bty = "n", border = NA)
  dev.off()  
}


### top rich genera ------------------------------------------------------------
top_rich_genera <- function(data, n = 10) { 
  data_genus <- data[, c("genus", "species")]
  data_genus <- unique(data_genus)
  sp_x_genus <- sort(table(data_genus$genus), decreasing = T)
  write.table(sp_x_genus, file = "Output/sp_x_genus.csv", col.names = NA,
    row.names = T, sep = ";")
  output <- head(sp_x_genus, n)
  return(output)   
}


### summary of species ecology -------------------------------------------------
species_ecology <- function(data) {
  data_ecology <- data[, c("ecology", "species")]
  data_ecology <- unique(data_ecology)
  double_ecology <- 
    names(table(data_ecology$species)[table(data_ecology$species) > 1])
  for(i in double_ecology) {
    data_ecology[data_ecology$species == i, "ecology"] <-
      "saprotrophic or parasitic on vascular plants"
  }
  data_ecology <- unique(data_ecology)
  data_ecology <<- droplevels(data_ecology)
  
  sp_x_ecology <- sort(table(data_ecology$ecology), decreasing = T)
  write.table(sp_x_ecology, file = "Output/sp_x_ecology.csv",
    col.names = NA, row.names = T, sep = ";")
  
  sp_x_ecology <- sp_x_ecology[c(1, 4, 2:3, 5:18)]
  unk <- which(names(sp_x_ecology) == "unclear")
  sp_x_ecology <- sp_x_ecology[c(1:(unk - 1),
    (unk + 1):length(sp_x_ecology), unk)] 
  col_ecology <- color(length(sp_x_ecology))
  pdf("Output/sp_x_ecology.pdf", h = 3, w = 6, pointsize = 14)
  par(mar = rep(2, 4), mfrow = c(1, 2))
  pie(sp_x_ecology, labels = NA, col = col_ecology, border = "white")
  par(mar = rep(0, 4))
  plot(1, 1, type = "n", axes = F)
  legend("left", legend = names(sp_x_ecology), bty = "n", fill = col_ecology,
    border = NA, y.intersp = 0.75)
  dev.off()  
}


### summary of lichen substrata ------------------------------------------------
lichen_substrata <- function(data) {
  lichen <- data[data$ecology == "lichenized", "substrate"]
  lichen <- droplevels(lichen)
  lichen <- as.character(lichen)
  for(i in grep(",", lichen)) {
    x <- unlist(strsplit(lichen[i], ", "))
    lichen <- c(lichen, x)
  }
  lichen <- lichen[-grep(",", lichen)]
  lichen <- trimws(lichen, which = c("both"))  
  tab <- sort(table(lichen), decreasing = T)
  unk <- which(names(tab) == "unknown")
  tab <- tab[c(1:(unk - 1), (unk + 1):length(tab), unk)] 
  col_lichen <- color(length(tab), start = 55)
  pdf("Output/lichen_substrata.pdf", h = 3, w = 6, pointsize = 14)
  par(mar = rep(2, 4), mfrow = c(1, 2))
  pie(tab, labels = NA, col = col_lichen, border = "white")
  par(mar = rep(0, 4))
  plot(1, 1, type = "n", axes = F)
  legend("left", legend = names(tab), bty = "n", fill = col_lichen, border = NA,
    y.intersp = 0.75)
  dev.off()  
}


### comparison with GBIF data --------------------------------------------------
gbif_comparison <- function(data, file = "Data/GBIF_records.csv") {
  
  gbif <- read.csv(file, h = T, sep = ";")
  gbif_spp <- unique(gbif$species)
  chckl_spp <- unique(data$species)
  
  a <- sum(!(chckl_spp %in% gbif_spp))
  b <- sum(gbif_spp %in% chckl_spp)
  c <- sum(!(gbif_spp %in% chckl_spp))
  
  message(paste("Species only in checklist", a))
  message(paste("Species only in GBIF", c))
  message(paste("Species in both", b))
  
  pdf("Output/GBIF_comparisons.pdf", w = 4, h = 3, pointsize = 14)
  par(las = 1, mar = c(4, 4, 1, 8), xpd = T, 
    tck = -0.025, mgp = c(3, 0.6, 0.6))
  barplot(rbind(a, b, c), names.arg = "GBIF", col = c(gray(0.9), 1, gray(0.7)),
    ylab = "n species", border = F)
  legend("topright", legend = c("GBIF", "both", "checklist"), border = F,
    fill = c(gray(0.7), 1, gray(0.9)), inset = c(-1.2, 0), bty = "n")
  dev.off()  
  gbif_species <- c(as.character(gbif_spp[gbif_spp %in% chckl_spp]),
    as.character(gbif_spp[!(gbif_spp %in% chckl_spp)]))
  gbif_species <- as.data.frame(gbif_species)
  gbif_species$in_checklist <- c(rep(T,
    length(gbif_spp[gbif_spp %in% chckl_spp])),
    rep(F, length(gbif_spp[!(gbif_spp %in% chckl_spp)])))  
  write.table(gbif_species,
    file = "Output/GBIF_spp.csv", row.names = F, col.names = F)
}


### function to generate distinct colors ---------------------------------------
# adapted from https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
color <- function(n, start = 1, random = F, last_gray = T) {
  require(RColorBrewer)
  
  # 433 colors
  # x <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
  # 74 colors
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  x <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
    rownames(qual_col_pals)))

  result <- x[start:(start + n - 1)]
  if(random) result <- sample(x, n)
  if(last_gray) result[length(result)] <- gray(0.75)
  
  return(result)  
}


### get country centroids ------------------------------------------------------
get_centroid <- function(country) {
  mp <- map("world", country, plot = F)
  x_cent <- mean(c(min(na.omit(mp$x)), max(na.omit(mp$x))))
  y_cent <- mean(c(min(na.omit(mp$y)), max(na.omit(mp$y))))
  centroid <- c(x_cent, y_cent)
  names(centroid) <- c("x", "y")
  return(centroid)
}


### numbers of type specimens --------------------------------------------------
type_specimens <- function(data) {
  types <- droplevels(data[grep("type", data$modifier), ])
  total <- nrow(types)
  types_x_country <- colSums(table(types$species, types$country))
  types_x_country <- c(total, types_x_country)
  names(types_x_country)[1] <- "Total WA"
  write.table(types_x_country, file = "Output/types_x_country.csv",
    col.names = NA, row.names = T, sep = ";")
  return(types_x_country) 
}


### number of incertae sedis groups --------------------------------------------
incertae_sedis_taxa <- function(data) {
  is_spp <- data[unique(c(grep("sedis", data$order),
    grep("sedis", data$family))), "species"]
  is_spp <- droplevels(is_spp)
  message(paste("Number of records with an incertae sedis taxon:",
    length(is_spp))) 
  is_spp <- unique(is_spp)
  message(paste("Number of species with an incertae sedis taxon:",
    length(is_spp)))
  is_div <- droplevels(data[data$species %in% is_spp, c("division", "species")])
  is_div <- unique(is_div)
  is_x_division <- table(is_div$division)
  return(is_x_division)
}

