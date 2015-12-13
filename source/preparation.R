# # load the structure fingerprint
# source(paste(getwd(), "/source/customized.R", sep=""), local=TRUE)
# structure_fp_base <- 'U:/Projects/TOX21/Chemical_Curation_from_Ann/20140722/tox21_v5a_leadscope_fp_extend' # tox21_8598_fp tox21_8306_fp
# struct_mat <- load_struc_fp_file(structure_fp_base,NULL) 
# save(struct_mat, file='struct_mat.RData')

# pathg <- 'INVITRODB_V2_SUMMARY'
# 
# ll <- list.files(file.path(getwd(),pathg), '.*_Matrix_151020.csv')
# names(ll) <- sub("_Matrix_151020.csv", "", ll, fixed=TRUE)
# ll <- ll[! names(ll) %in% c("cyto_dist")]  #which(names(ll) %in% c("cyto_dist"))  the format is different
# 
# #mapping <-  read_excel(file.path(getwd(),  'DSSTox_ToxCastRelease_20151019', 'DSSTox_ToxCastRelease_20151019.xlsx'), sheet=1) # doesn't have some IDs 
# mapping <- read_excel(file.path(getwd(),pathg, 'Chemical_Summary_151020.xlsx'), sheet=1) # add right CAS column
# mapping <- mapping %>% mutate(GSID=str_c("gsid_", mapping$chid))
# gsid <- str_c("gsid_", mapping$chid) %>% setNames(mapping$code) # chid is relevant to GSID
# 
# lla <- sapply(ll, function (x) {
#   csv <- tbl_df(read.csv(file.path(getwd(),pathg, x)))
#   try({
#     rownames(csv) <- gsid[csv$X]
#     csv <- select(csv, -X)
#   })
#   return(csv)
# }
# , simplify = FALSE, USE.NAMES = TRUE )
# 
# # make sure columns are identical
# sapply(lla, function (x) identical(rownames(lla$fitc), rownames(x)))
# sapply(lla, function (x) identical(colnames(lla$fitc), colnames(x)))
# 
# toxcast_acts <- lla
# activities <- lapply(toxcast_acts, function (x) as.data.frame(x)) %>% setNames(names(toxcast_acts))
# save(toxcast_acts, file="toxcast_acts_102015.RData")
# 
# # load the flags
# flags_file <- 'AllResults_flags_151020.csv'
# flags <- tbl_df(read.csv(file.path('U:/Projects/TOXCAST/102815/INVITRODB_V2_SUMMARY', flags_file)))
# flags <- flags %>% mutate(GSID = str_c("gsid_", chid)) %>% group_by(GSID,spid, aenm) %>% summarize(flag_c = str_c(flag, collapse="|"))
# 
# spids <- toxcast_acts[['spid']]
# flags_spids <- spids %>% add_rownames("GSID") %>% gather(aenm, spid, -GSID) %>% 
#   left_join(flags) %>% mutate(flag_c=ifelse(is.na(flag_c), "", flag_c))
# 
# flags_spids[, "GSID"] <- ordered(flags_spids$GSID, levels=rownames(toxcast_acts[[1]]))
# flags_spids[, "aenm"] <- ordered(flags_spids$aenm, levels=colnames(toxcast_acts[[1]]))
# 
# result <- flags_spids %>% select(-spid) %>% spread(aenm, flag_c, drop=FALSE, fill="") 
# rownames(result) <- as.character(result$GSID)
# toxcast_acts[["flags"]] <- as.data.frame(select(result, -GSID))
# 
# 
# # toxcast data files (simplified)
# load("U:/Projects/TOXCAST/102815/toxcast_acts_102015.RData")
# include <- c('hitc', 'modl_ac10', 'modl_acb', 'modl_acc', 'modl_ga', 'modl_tp', 'flags')
# activities <- toxcast_acts[include]
# activities <- lapply(names(activities), function (x) {
#   if (x %in% c('modl_ac10', 'modl_acb', 'modl_acc', 'modl_ga'))
#   {
#      activities[[x]] <- (activities[[x]] - 6)*-1
#   }
#   return(activities[[x]])
# }) %>% setNames(names(activities))
# 
# # cyto lower bound # there are 367 GSID has no lower bnd value
# cyto_dist_file <- 'U:/Projects/TOXCAST/102815/INVITRODB_V2_SUMMARY/cyto_dist_Matrix_151020.csv'
# cyto_dist <- read.csv(cyto_dist_file) %>% 
#   mutate(GSID = paste0("gsid_", chid), log_lower_bnd=log10(lower_bnd_um/1000000)*-1) %>%
#   select(GSID, log_lower_bnd)
# cyto_dist_v <- data.frame(GSID=rownames(activities[[1]])) %>% left_join(cyto_dist)
# lower_bnd <- as.data.frame(apply(activities[[1]], 2, function(x) x <- cyto_dist_v$log_lower_bnd))
# rownames(lower_bnd) <- rownames(activities[[1]])
# activities[['cyto_lower_bnd']] <- lower_bnd
# 
# # scaled resp
# assay_quality_file <- 'U:/Projects/TOXCAST/102815/INVITRODB_V2_SUMMARY/Assay_Quality_Summary_Stats_151020.csv'
# assay_q <- read.csv(assay_quality_file) %>% 
#   arrange(aenm) %>%
#   select(aenm, coff) 
# assay_q_v <- data.frame(aenm=colnames(activities[[1]])) %>% left_join(assay_q)
# coff <- as.data.frame(t(apply(t(activities[[1]]), 2, function(x) x <- assay_q_v$coff))) 
# colnames(coff) <- colnames(activities[[1]])
# activities[['scaled_emax']] <- activities$modl_tp/coff
# activities <- lapply(activities, function (x) as.data.frame(x)) %>% setNames(names(activities))
# save(activities, file='activities.RData')