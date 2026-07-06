rank_parent_sets <- function(df, number_of_parents = 5, rel_cutoff = 0.0) {
  
  df <- df %>%
    mutate(
      Pmin = pmin(Parent1, Parent2),
      Pmax = pmax(Parent1, Parent2),
      key = paste(Pmin, Pmax, sep = "_")
    )
  
  U_lookup <- setNames(df$Y, df$key)
  GRM_lookup <- setNames(df$K, df$key)
  
  parents <- sort(unique(c(df$Parent1, df$Parent2)))
  sets <- combn(parents, number_of_parents, simplify = FALSE)
  
  evaluate_set <- function(S) {
    pairs <- combn(S, 2, simplify = FALSE)
    keys <- purrr::map_chr(
      pairs,
      ~ paste(pmin(.x[1], .x[2]), pmax(.x[1], .x[2]), sep = "_")
    )
    
    rels <- GRM_lookup[keys]
    if (any(is.na(rels))) return(NULL)
    if (any(rels > rel_cutoff)) return(NULL)
    
    Uvals <- U_lookup[keys]
    if (any(is.na(Uvals))) return(NULL)
    
    tibble::tibble(
      Parents = paste(S, collapse = " / "),
      Mean = mean(Uvals),
      Min = min(Uvals),
      Max = max(Uvals),
      Max_Relationship = max(rels)
    )
  }
  
  results <- purrr::map_dfr(sets, evaluate_set)
  results %>% arrange(desc(Mean))
}