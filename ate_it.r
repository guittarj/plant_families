# Gastauer, Markus, Meira Neto, and João Augusto Alves. "Updated angiosperm family tree for analyzing phylogenetic diversity and community structure." Acta Botanica Brasilica 31.2 (2017): 191-198.

library(dplyr)
library(phytools)
library(tidytree)
library(collapsibleTree)

phylo <- read.tree('plant_families_Gastauer2017.new')
x <- as_data_frame(phylo)

ate_raw <- read.csv('ate_it.csv', stringsAsFactors = FALSE)
ate <- ate_raw %>%
  mutate_all(trimws) %>%
  mutate(family = tolower(family), species = tolower(species))
fams <- unique(ate$family[ate$family != ''])
if(!all(fams %in% unique(ate$all_families))) 
  stop(paste("Not all families eaten are present in tree. Check your spelling of:", fams[!fams %in% unique(ate$all_families)]))
ate <- ate %>%
  filter(family %in% fams) %>%
  mutate(label = paste0(genus, ' ', species, ' (', common_name, ')'))

x$parent_name <- x$label[match(x$parent, x$node)]
x$parent_name[x$label == 'angiosperms'] <- NA

a <- list()
ups <- for(i in fams) {
  #location
  k <-1
  #name
  j <- x$parent_name[x$label == i]
  a[[i]] <- c()
  while (!is.na(j)) {
    a[[i]][[k]] <- j
    k <- k + 1
    j <- x$parent_name[x$label == j]
  }
  a[[i]] <- a[[i]][!grepl("_to_", a[[i]], fixed = TRUE)]
}

j <- data.frame(parent = NA, child = 'angiosperms')

for(i in names(a)) {
  len <- length(a[[i]])
  tmp <- data.frame(parent = a[[i]],
                    child = c(i, a[[i]][1:len-1]))
  j <- rbind(j, tmp)
}

j <- rbind(
  mutate(j,
    color = ifelse(child %in% ate_raw$all_families, as.character(child), NA),
    notes = ''),
  transmute(ate, 
    color = family, 
    parent = family, 
    child = label,
    notes = notes)) %>%
  mutate(
    color = as.numeric(as.factor(color)),
    color = ifelse(is.na(color), '#FFFFFF',
                   colorspace::rainbow_hcl(length(unique(color[!is.na(color)])))[color]))

#might be possible to include hover-over images using HTML
#j$tooltip <- NA
#j$tooltip[j$family == 'caricaceae'] <- "<img src='http://www.tizag.com/pics/htmlT/sunset.gif'/>"

collapsibleTreeNetwork(j, 
  attribute = 'notes', fill = 'color', 
  collapsed = FALSE, zoomable = FALSE, nodeSize = 'leafCount'
  #,tooltipHtml = "tooltip"
  )

