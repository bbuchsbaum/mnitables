

harvard_oxford <- function(cds) {
  cort_labels <- read.table(system.file("extdata", "cortical_table.txt", 
                                        package = "mnitables"), header=TRUE)
  subcort_labels <- read.table(system.file("extdata", "subcortical_table.txt", 
                                           package = "mnitables"), header=TRUE)
  
  cort_atlas <- 
    neuroim2::read_vol(system.file("extdata", "HarvardOxford-cort-maxprob-thr0-2mm.nii.gz", 
                                   package = "mnitables"))
  
  subcort_atlas <- 
    neuroim2::read_vol(system.file("extdata", "HarvardOxford-sub-maxprob-thr0-2mm.nii.gz", 
                                   package = "mnitables"))
  
  vox2 <- coord_to_grid(cort_atlas, cds)
  
  ltab <- do.call(rbind, lapply(1:nrow(vox2), function(i) {
    v <- vox2[i,]
    
    clab <- cort_atlas[v[1], v[2], v[3]]
    slab <- subcort_atlas[v[1], v[2], v[3]]
    
    label <- if (clab == 0 && slab == 0) {
      "unknown"
    } else if (clab >= 1 || clab <= 48) {
      idx <- which(cort_labels$index == (clab-1))
      as.character(cort_labels$regions[idx])
    } else {
      idx <- which(subcort_labels$index == (slab-1))
      as.character(subcort_labels$regions[idx])
    }
    
    hemi <- if (cds[i,1] <= 0) {
      "left"
    } else {
      "right"
    }
    
    data.frame(x=cds[i,1], y=cds[i,2], z=cds[i,3], Hemi=hemi, Label=label)
  }))
  
  
}


#' @export
create_table <- function(im, threshold=0, local_maxima=TRUE, local_maxima_dist=10, statname="z-stat", ...) {
  ccomp <- conn_comp(im, local_maxima_dist=local_maxima_dist)
  
  local_maxima <- as.data.frame(ccomp$local_maxima)
  vox <- ccomp$local_maxima[,2:4]
  values <- ccomp$local_maxima[,5]
  area <- ccomp$cluster_table$Area[local_maxima$index]
  
  cds <- grid_to_coord(im, vox)
  
  ltab <- harvard_oxford(cds)
  ltab[[statname]] <- values
  ltab$Area <- area
  ltab
}