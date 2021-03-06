
#' @import stringr
#' @import dplyr
harvard_oxford <- function(cds) {
  cort_labels <- read.table(system.file("extdata", "cortical_table.txt", 
                                        package = "mnitables"), header=TRUE,stringsAsFactors=FALSE)
  subcort_labels <- read.table(system.file("extdata", "subcortical_table.txt", 
                                           package = "mnitables"), header=TRUE, stringsAsFactors=FALSE)
  
  subcort_labels <- subcort_labels %>% mutate(regions = stringr::str_trim(regions))
  cort_labels <- cort_labels %>% mutate(regions = stringr::str_trim(regions))
  
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
    } else if (clab >= 1 && clab <= 48) {
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


#' create a labeled coordinate table
#' 
#' 
#' @param im a\code{NeuroVol} instance
#' @param a cluster threshold value to define the mask. 
#' @param local_maxima_dist the distance threshold used to define local_maxima
#' @param statname the name of the statistic associated with the image values
#' @param allow_duplicate_labels whether to allow the same atlas label to appear more than once in the table. 
#'        If \code{FALSE} and the same label appears twice, then only the coordinate with the highest vlaue will be retained.
#' @param ... extra args to pass to \code{conn_comp}
#' @export
#' @import dplyr
create_table <- function(im, threshold=0, local_maxima_dist=10, statname="zstat", allow_duplicate_labels=TRUE, ...) {
  ccomp <- conn_comp(im, local_maxima_dist=local_maxima_dist,...)
  
  local_maxima <- as.data.frame(ccomp$local_maxima)
  vox <- ccomp$local_maxima[,2:4]
  values <- ccomp$local_maxima[,5]
  area <- ccomp$cluster_table$Area[local_maxima$index]
  
  cds <- grid_to_coord(im, vox)
  
  ltab <- harvard_oxford(cds)
  ltab[[statname]] <- values
  ltab$Area <- area
  ltab$Voxels <- ccomp$cluster_table$N[local_maxima$index]
  
  if (!allow_duplicate_labels) {
    col_name <- rlang::sym(statname)
    ltab <- ltab %>% arrange(desc(UQ(col_name)), Label) %>% 
      group_by(Label) %>% filter(row_number() == 1) %>% arrange(desc(Area))
  }
  ltab
}