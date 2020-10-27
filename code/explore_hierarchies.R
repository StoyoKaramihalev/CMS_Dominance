# Packages:
# - Compete: Curley, J.P. 2016, compete: Analyzing Social Hierarchies: R package version 0.1
# - Steepness: Leiva, D & de Vries, H. Steepness: Testing steepness of dominance hierarchies (2014)

# Note: this is an older version of this function which still deals with the master data at the level of individuals. The newer version used in conjunction with the `posr` & `socialbox` packages will work with daily groupwise entries. 

# Note: if an idVars entry is provided that separates individuals within groups, expect
# this to malfunction

require(compete)
require(steepness)

explore_hierarchies <- function(master, 
                                nPerm = 1000,
                                idVars = NULL, 
                                daily = TRUE,
                                DS_only = FALSE,
                                verbose = FALSE
) {
  
  stopifnot("Suspect column is missing from `master`" = 
              "Suspect" %in% colnames(master))
  stopifnot("Chase column is missing from `master`" = 
              "Chases" %in% colnames(master))
  stopifnot("Some of the `idVars` supplied are not present in `master`" = 
              all(idVars %in% colnames(master)))
  
  idVars <- c('MouseID', 'Day', 'GroupNumber', 'GroupType', idVars)
  
  suspect <- master %>% 
    dplyr::select(tidyselect::all_of(idVars), Suspect) 
  suspect_other <- suspect %>%
    dplyr::rename(OtherID = MouseID, SuspectOther = Suspect)
  
  master <- master %>% 
    dplyr::select(tidyselect::all_of(idVars), Chases) %>% 
    tidyr::unnest(cols = c(Chases)) %>% 
    dplyr::rename(wins = Chases) %>% 
    dplyr::left_join(suspect) %>% 
    dplyr::left_join(suspect_other)
  
  master$wins[master$Suspect | master$SuspectOther] <- NA
  # master <- master %>% filter(!master$Suspect & !master$SuspectOther)
  master$wins[master$MouseID == master$OtherID] <- NA
  
  master <- master %>% select(-Suspect, -SuspectOther)
  
  if(daily) {
    if(length(unique(master$Day)) == 1) message("Only one day is present in the data")
    master <- master %>% 
      dplyr::group_by_at(.vars = tidyselect::all_of(idVars[idVars != 'MouseID'])) %>% 
      tidyr::nest() %>% 
      dplyr::mutate(
        mat = purrr::map(data, function(d) {
          d <- d %>% 
            tidyr::pivot_wider(names_from = OtherID, 
                               values_from = wins,
                             names_prefix = "other_")
          m <- d[,2:ncol(d)] %>% as.matrix()
          rownames(m) <- d$MouseID
          return(m)
        }),
        n_subjects = purrr::map_dbl(mat, ~ifelse(is.null(.x), NA, dim(.x)[1])), 
        NAs = purrr::map_dbl(mat, ~ sum(is.na(.x[lower.tri(.x)])))
      )
  }
  
  # Cumulative
  if(!daily) {
    master <- master %>% 
      dplyr::group_by_at(.vars = c(tidyselect::all_of(idVars[!idVars == 'Day']), 'OtherID')) %>% 
      dplyr::summarize(
        wins = mean(wins, na.rm = T)
        ) %>% 
      tidyr::pivot_wider(names_from = OtherID, values_from = wins,
                  names_prefix = "other_") %>% 
      dplyr::group_by_at(.vars = c(tidyselect::all_of(idVars[!idVars %in% c('Day', 'MouseID')]))) %>% 
      nest() %>%
      mutate(
        mat = purrr::map(data, function(d) {
          m <- d %>% 
            dplyr::select(dplyr::contains("other_")) %>% 
            as.matrix()
          rownames(m) <- d$MouseID
          return(m)
          }),
        n_subjects = purrr::map_dbl(mat, ~ifelse(is.null(.x), NA, dim(.x)[1])), 
        NAs = purrr::map_dbl(mat, ~ sum(is.na(.x[lower.tri(.x)])))
      )
    
    # remove mat rows & cols with all NAs? 
    # master <- master %>% mutate(
    #   mat = map(mat, ~.x[!rowSums(is.na(.x)) == 4, !rowSums(is.na(.x)) == 4])
    # )
  }
  
  if(verbose) message("Calculating David's Scores")
  master <- master %>%
    mutate(
      DS = map(mat, function(m) 
        data.frame(
          MouseID = rownames(m),
          DS = compete::ds(m = m, norm = T, type = 'D') 
        )),
      DS_Pij = map(mat, function(m) 
        data.frame(
          MouseID = rownames(m),
          DS_Pij = compete::ds(m = m, norm = T, type = 'P') 
        ))
    )
  
  if(!DS_only) {
    
    if(verbose) message("Calculating Landau's h'")
    master <- master %>%
      mutate(
        Landau = map_dbl(mat, function(m) compete::devries(m, Nperms = nPerm, history = F, plot = F)$`h-modified` ),
        Landau_pval = map_dbl(mat, function(m) compete::devries(m, Nperms = nPerm, history = F, plot = F)$`p-value` )
      )
    
    if(verbose) message("Calculating Phi")
    master <- master %>%
      mutate(
        Phi = map_dbl(mat, function(m) compete::phi(m))
      )
    
    if(verbose) message("Calculating Triangle Transitivity")
    master <- master %>%
      mutate(
        TriangleTransitivity = map2_dbl(mat, n_subjects, function(m, n) 
          if(n > 3) compete::ttri(m)$ttri else NA
        ),
        TriangleTransitivity_pval = map2_dbl(mat, n_subjects, function(m, n) 
          if(n > 3) compete::ttri_test(m, Nperms = nPerm)$pval else NA
        )
      )
    
    if(verbose) message("Calculating Despotism")
    # Fraction of total chases by the dominant animal
    master <- master %>%
      mutate(
        Despotism = map2_dbl(mat, DS, function(m, DS) {
          alpha = as.numeric(as.character(DS$MouseID[which(DS$DS == max(DS$DS, na.rm = T))]))
          Despotism = sum(m[alpha,], na.rm = T) / sum(m, na.rm = T)
        })
      )
    
    if(verbose) message("Calculating Directional Consistency")
    master <- master %>%
      mutate(
        DirectionalConsistency = map_dbl(mat, function(m) compete::dci(m))
      )
    
    if(verbose)  message("Calculating Steepness")
    master <- master %>%
      mutate(
        Steepness = map_dbl(mat,function(m) {
          # Note - requires matrix with zeros off the diag
          diag(m) = 0 
          if(sum(is.na(m)) > 0) { if(verbose) warning('mat contains NAs') 
            return(NA)
          } else {
            steepness::getStp(m)
          }
        })
      )
  }
  
  
  return(master)
}
