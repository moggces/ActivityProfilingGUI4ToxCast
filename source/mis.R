# to convert mapping from inp to out
conversion <- function (master, inp, out)
{
  result <- master[, out]
  names(result) <- master[, inp]
  return(result)
}

# split the master file into a list of matrices
split_master_2_matrix <- function(master, props, id='CAS')
{
  id_data <- master[, id]
  result <- lapply(as.list(props), function (x)
  {
    
    n_ori <- str_count(x, '\\.')
    #col_ids <- grepl(paste('^',x, sep=""), colnames(master))
    col_names <- grep(paste('^',x, '\\.', sep=""), colnames(master), value=TRUE)
    col_names <- col_names[str_count(col_names, '\\.') == n_ori + 1]
    
    mat <- master[,col_names]
    rownames(mat) <- id_data
    colnames(mat) <- sub(paste(x, '.', sep=""), "", colnames(mat))
    return(mat)
  }
  )
  names(result) <- props
  return(result)
}

# required: conversion, rename the GSID to Chemical.Name (row), rename assay to common_name (col, if TRUE)
rename_mat_col_row <- function (partial, master, assay_names, input_chemical_name=NULL, rename_assay=TRUE)
{
  chemical_name_ref <- input_chemical_name
  if (is.null(chemical_name_ref)) {
    chemical_name_ref <- conversion(master, inp='GSID', out='Chemical.Name')
  } 

  for (name in names(partial))
  {
    #print(name)
    #print(rownames(partial[[name]]))
    #print(chemical_name_ref[as.character(rownames(partial[[name]])) ])
    rownames(partial[[name]]) <-  chemical_name_ref[as.character(rownames(partial[[name]])) ]
    if (  name != 'struct' & rename_assay )
    {
      
      pathway_ref <- conversion(assay_names, inp='assay', out='common_name')
      colnames(partial[[name]]) <-  pathway_ref[as.character(colnames(partial[[name]])) ]
    }
  }
  
  return(partial)
}

# sort the matrices in a list based on rownames and columnames
sort_matrix <- function (partial)
{
  for (name2 in names(partial))
  {
    if (name2 == 'struct' )
    {
      partial[[name2]] <- partial[[name2]][order(rownames(partial[[name2]])),]
      
    } else
    {
      partial[[name2]] <- partial[[name2]][order(rownames(partial[[name2]])),]
      partial[[name2]] <- partial[[name2]][,order(colnames(partial[[name2]]))]
    }
  }
  return(partial)
}

# make increasing resp in mitotox as active,  dependent on mitotox name: inhibition_MMP
fix_mitotox_reverse <- function(partial, act_mat_names=c('npod', 'nec50', 'nwauc.logit'))
{
  for (name in act_mat_names)
  {
    mitotox_id <- grepl('inhibition_MMP', colnames(partial[[name]]))
    if (sum(mitotox_id) > 0)
    {
      rev_ids <- partial[[name]][, mitotox_id] < 0 &  ! is.na(partial[[name]][, mitotox_id])
      partial[[name]][rev_ids, mitotox_id] <- (partial[[name]][rev_ids, mitotox_id])*-1
    }
  }
  return(partial)
}

# use the information of assay_names (antagonism|inhibition)
filter_activity_by_type <- function(partial, type, thres=NULL, decision=FALSE, act_mat_names=c('modl_acc', 'modl_acb', 'modl_ga', 'modl_ac10'))
{
  #print(partial[['cyto_lower_bnd']][,1, drop=FALSE])
  for (name in act_mat_names)
  {
    
    if (type == 'cyto_lower_bnd' & isTRUE(decision)) 
    {
      
      sel_assays <- assay_names %>% 
        filter(burst_assay != 1 ) %>% #burst_assay
        select(assay_component_endpoint_name) %>% unlist()
      non_cyto_assay_ids <- colnames(partial[[name]]) %in% sel_assays
      ids <- matrix(FALSE, nrow(partial[[name]][, non_cyto_assay_ids]), ncol(partial[[name]][, non_cyto_assay_ids]))
      ids <-  partial[[name]][,non_cyto_assay_ids] <= partial[[type]][,non_cyto_assay_ids]  & ! is.na(partial[[type]][,non_cyto_assay_ids]) & ! is.na(partial[[name]][,non_cyto_assay_ids]) & partial[[name]][,non_cyto_assay_ids] > 0.0001
      partial[[name]][, non_cyto_assay_ids][ids] <- (partial[[name]][, non_cyto_assay_ids][ids])*-1
      
    } else
    {
      ids <- matrix(FALSE, nrow(partial[[name]]), ncol(partial[[name]]))
      if (type %in% c(act_mat_names, 'scaled_emax')) ids <- partial[[type]] < thres & ! is.na(partial[[type]]) & ! is.na(partial[[name]]) & partial[[name]] > 0.0001
      if (type == 'hitc' & isTRUE(decision)) ids <- partial[[type]] != 1 & ! is.na(partial[[type]]) & ! is.na(partial[[name]]) & partial[[name]] > 0.0001
      if (type == 'flags' & isTRUE(decision)) ids <-  partial[[type]] != '' & ! is.na(partial[[type]]) & ! is.na(partial[[name]]) & partial[[name]] > 0.0001  
      partial[[name]][ids] <- (partial[[name]][ids])*-1
    }
  }
  return(partial)
}

# make < 0 as inconclusive (0.0001), make NA in potency as 0 if inactive, (dependent on nwauc.logit matrix)
assign_reverse_na_number <- function (partial, act_mat_names=c('modl_acc', 'modl_acb', 'modl_ga', 'modl_ac10'))
{
  result <- partial
  for (name in act_mat_names)
  {
  
    result[[name]][ partial[['hitc']] == 0 & ! is.na(partial[['hitc']]) ] <- 0
    result[[name]][ result[[name]] < 0 |  is.na(result[[name]])   ] <- 0.0001

  }
  
  return(result)
}

# remove inconclusive label (0.0001) but keep the untested as 0.0001
remove_inconclusive_label <- function (partial, act_mat_names=c('modl_acc', 'modl_acb', 'modl_ga', 'modl_ac10'))
{
  result <- partial
  for (name in names(partial))
  {
    if (name %in% act_mat_names) 
    {result[[name]][ partial[[name]] == 0.0001 & ! is.na(partial[['hitc']]) ] <- 0}
  }
  return(result)
}

remove_nohit_assays <- function (partial, act_mat_names=c('modl_acc', 'modl_acb', 'modl_ga', 'modl_ac10'))
{
  result <- partial
  for (name in names(partial))
  {
    if (name %in% act_mat_names) 
    {
      # version #1
      #ids <- which(! (colSums(partial[[name]] == 0.0001) == nrow(partial[[name]]) | 
      #               colSums(partial[[name]] == 0) == nrow(partial[[name]])))
      #print(ids)
      ids <- which(colSums(partial[[name]] > 0.0001) > 0 )
      if (length(ids) > 0) result[[name]] <- partial[[name]][,ids]
    }
  }
  return(result)
}
