average_dfs <- function(df_list) { 
  
  combined_df <- dplyr::bind_rows(local_importance) |>  
    dplyr::group_by(variable) |> 
    dplyr::summarize(dplyr::across(where(is.numeric), mean, na.rm = TRUE), .groups = 'drop') 
  
  return(combined_df) 
  
}