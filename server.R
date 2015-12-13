
# shiny note: 1) it can't discriminate R vs. r in the r script file 
# 2) the renderTable has trouble with encoding issue (cannot recognize ppp file iconv -t UTF-8 -f ISO-8859-1)
# looks like the new read.table can automatically select best encoding
# 3) in shiny server, once you delete a file but replace a file with same name. somwhow don't know how to refresh its
# but if you update .R file, you can refresh to get new functions

# chemical_loader() out: list(id, ?(nwauc.logit or npod or nec50 or unknown))
# matrix_subsetter() out: list(activities or ?(nwauc.logit or npod or nec50 or unknown), struct)
# activity_filter() out: same as above 
# matrix_editor() out: list(nwauc.logit, npod, nec50,wauc.logit, struct, cv_mark)
# heatmap_para_generator() out: list(dcols, drows, annotation, annt_colors, act=act, struct, cv)
# select_plot()

# todo:
# 1. download potency plot
# 2. broaden the "unknown" color scheme

library(shiny)
library(plyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(dplyr)
#library(tidyr)
library(stringr)
options(stringsAsFactors = FALSE)


#Sys.setlocale(locale="C")
#setwd("~/ShinyApps/profiling/")
source(paste(getwd(), "/source/customized.R", sep=""), local=TRUE)
source(paste(getwd(), "/source/pheatmap_display_number.R", sep=""), local=TRUE)
source(paste(getwd(), "/source/get.R", sep=""), local=TRUE)
source(paste(getwd(), "/source/load.R", sep=""), local=TRUE)
source(paste(getwd(), "/source/mis.R", sep=""), local=TRUE)
#environment(pheatmap_new_label) <- environment(pheatmap) pheatmap v. < 1.0

# load assay related parameters
logit_para_file <- './data/Assay_Summary_151020.csv'
assay_names <- tbl_df(read.csv(logit_para_file)) # global, dataframe output

# load chemical information (will include purity later)
profile_file <- './data/tox21_compound_id.txt' #colunm name has to be GSID # v5a3
master <- load_profile(profile_file) # global, dataframe output

# load the activities (all data) and the structure fp matrix
activities_rdata <- './data/activities.RData'
struct_mat_rdata <- './data/struct_mat.RData'
load(activities_rdata, verbose=TRUE) # global, matrix output, activities
load(struct_mat_rdata, verbose=TRUE) # global, matrix output, struct_mat

# remove the structures with low purity
#struct_mat <- struct_mat[rownames(struct_mat) %in% rownames(activities[[1]]),]
# very weird!! this line causes no error frozen on shiny public server

# heatmap settings
# the negative direction breaks won't capture wauc with very small values
wauc_breaks <- c( -1, -0.75, -0.5, -0.25, -0.1, -0.02, 0, 0.0001, 0.1, 0.25, 0.5, 0.75, 1) # upper is filled , lower is empty 
wauc_colors <-  c("#053061" ,"#2166AC" ,"#4393C3" ,"#92C5DE", "#D1E5F0", "#F7F7F7", "gray", "#FDDBC7" ,"#F4A582" ,"#D6604D" ,"#B2182B", "#67001F"  ) #RdBu
wauc_leg_breaks <- c(-1, -0.75, -0.5, -0.25,  0,   0.25, 0.5, 0.75, 1 )
wauc_leg_labels <- c("-1", "-0.75", "-0.5", "-0.25",  "0", "0.25", "0.5", "0.75", "1")

potency_breaks <- c(-0.02, 0, 0.0001, 4, 4.5, 5, 7.5, 9, 12)
potency_colors <- c("#F5F5F5", "gray", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30", "chartreuse") #BrBG
potency_leg_breaks <- c( 0,  4, 4.5, 5, 7.5, 9,12 )
potency_leg_labels <- c( "inactive",  "100uM", "30uM", "10uM", "0.3uM", "1nM", "1pM")

# potency_breaks <- c(-10, -9, -7.5, -5, -4.5, -4, -0.02, 0, 0.0001, 4, 4.5, 5, 7.5, 9, 10)
# potency_colors <- c("darkorange","#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "gray", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30", "chartreuse") #BrBG
# potency_leg_breaks <- c(-10, -9, -7.5, -5, -4.5, -4,  0,  4, 4.5, 5, 7.5, 9,10 )
# potency_leg_labels <- c("-10", "-9", "-7.5", "-5", "-4.5", "-4",  "0",  "4", "4.5", "5", "7.5", "9", "10")


shinyServer(function(input, output) {

# chemical_loader()
  chemical_loader <- reactive({
    
    result <- NULL
    path <- NULL
    
    # input file
    inFile <- input$file1
    # input textarea
    textdata <- input$cmpds
    
    if (! is.null(inFile)) { path <- inFile$datapath; filen <- inFile$name }
    if (textdata != '' ) result <- load_text_2_df(textdata)
    if (! is.null(path)) result <- load_data_matrix(path, filen) # as long as path or file has something it will override
   
    return(result)

  })

# matrix_subsetter()
  matrix_subsetter <- reactive({
    partial <- NULL
    reg_assay <- input$reg_assay 
    inv_sel_assay <- input$inv_sel_assay 
    reg_gene <- input$reg_gene 
    inv_sel_gene <- input$inv_sel_gene
    reg_source <- paste0(input$assay_source, collapse = "|")
    #rename_assay <- FALSE # use the assay_names df (no for ToxCat)
    
    # get all chemical information
    id_info <- chemical_loader()
    if (is.null(id_info)) return(NULL)
    #chem_id_df <- get_lookup_list(id_info[['id']], master)
    #ip <- subset(chem_id_df, ! is.na(StructureID), select=c(CAS, Cluster))
    
    # the basic identifies , GSID + Cluster
    #ip <- subset(chem_id_df, GSID != '' & CAS != '', select=c(GSID, Cluster))
    
    # collect all the matrices and store in full (list)
    
    full <- list()
    full <- activities 
    
    # if it is a data matrix input, only CAS ID is allowd
#    input_chemical_name <- NULL
    if (length(id_info) > 1) # for loading the data matrix function
    {
      full <- id_info[! grepl('id', names(id_info))]
      chemical_name_ref <- conversion(master, inp='CAS', out='GSID')
      #rownames(full[[1]]) <- chemical_name_ref[as.character(rownames(full[[1]]))]
      avail_name <- chemical_name_ref[as.character(rownames(full[[1]]))]
      full[[1]] <- full[[1]][! is.na(avail_name), ]
      rownames(full[[1]]) <- avail_name[!is.na(avail_name)]
      
#      if (! is.null(id_info[['id']]$input_Chemical.Name)) {
#        input_chemical_name <- conversion(join(id_info[['id']], master), inp='GSID', out='input_Chemical.Name')
#      }
#      rename_assay <- FALSE
    }
    
    # the structure fingerprint matrix
    full[['struct']] <- struct_mat
    
    # subset the matrices by chemicals
    #partial <- get_input_chemical_mat(ip, full)
    partial <- full
    
    # rename the assays & chemicals
#    partial <- rename_mat_col_row(partial,  master, assay_names, input_chemical_name, rename_assay=FALSE)
    
    # subset the matrices by assay names
    partial <- get_assay_mat(partial, sel=reg_source, inv=FALSE, type='source')
    partial <- get_assay_mat(partial, sel=reg_assay, inv=inv_sel_assay, type='assay' )
    partial <- get_assay_mat(partial, sel=reg_gene, inv=inv_sel_gene, type='gene' )
    
    # sort the matrix
#    partial <- sort_matrix(partial)
    
    return(partial)
  })

# activity_filter()  
  activity_filter <- reactive({
    
    # load all the activity filter parameters
    
    partial <- matrix_subsetter()
    if (is.null(partial)) return(NULL)
    
    # if it is data matrix input, don't change 
    if (length(partial) == 2) return(partial)
    
    sresp_thres <- input$sresp_thres
    acc_thres <- ifelse(is.na(input$acc_thres), 0, log10(input$acc_thres/1000000)*-1)
    acb_thres <- ifelse(is.na(input$acb_thres), 0, log10(input$acb_thres/1000000)*-1)
    ga_thres <- ifelse(is.na(input$ga_thres), 0, log10(input$ga_thres/1000000)*-1)
    ac10_thres <- ifelse(is.na(input$ac10_thres), 0, log10(input$ac10_thres/1000000)*-1)
    #isstrong <- input$isstrong
    hitfilter <- TRUE
    cytofilter <- input$cytofilter
    flagfilter <- input$flagfilter

    act_mat_names <- c('modl_acc', 'modl_acb', 'modl_ga', 'modl_ac10')
    # reverse direction of mitotox could be meaningful
    #partial <- fix_mitotox_reverse(partial,act_mat_names=act_mat_names )
    
    # filtering
    if (sresp_thres > 1) partial <- filter_activity_by_type(partial, 'scaled_emax', sresp_thres, act_mat_names=act_mat_names)
    if (acc_thres > 0) partial <- filter_activity_by_type(partial, 'modl_acc', acc_thres,act_mat_names=act_mat_names)
    if (acb_thres > 0) partial <- filter_activity_by_type(partial, 'modl_acb', acb_thres,act_mat_names=act_mat_names)
    if (ga_thres > 0) partial <- filter_activity_by_type(partial, 'modl_ga', ga_thres,act_mat_names=act_mat_names)
    if (ac10_thres > 0) partial <- filter_activity_by_type(partial, 'modl_ac10', ac10_thres,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'hitc', thres=NULL, decision=hitfilter,act_mat_names=act_mat_names)
    partial <- filter_activity_by_type(partial, 'cyto_lower_bnd', thres=NULL, decision=cytofilter,act_mat_names=act_mat_names)
    # it has to be the end
    partial <- filter_activity_by_type(partial, 'flags', thres=NULL, decision=flagfilter,act_mat_names=act_mat_names)    
    
    return(partial)
  })
  
  matrix_subsetter_chemical <- reactive({
    # get all chemical information
    id_info <- chemical_loader()
    if (is.null(id_info)) return(NULL)
    chem_id_df <- get_lookup_list(id_info[['id']], master)
    # the basic identifies , GSID + Cluster
    ip <- subset(chem_id_df, GSID != '' & CAS != '', select=c(GSID, Cluster))
    
    full <- activity_filter()
    
    # subset the matrices by chemicals
    partial <- get_input_chemical_mat(ip, full)
    return(partial)
    
  })

# matrix_editor()
  matrix_editor <- reactive({
    
    noincon_label <- input$noinconlab #inconclusive label
    rm_nohit_assay <- input$rm_nohit_assay
    act_mat_names <- c('modl_acc', 'modl_acb', 'modl_ga', 'modl_ac10')
    
#    partial <- activity_filter()
    partial <- matrix_subsetter_chemical()
    # if it is data matrix input, skip 
    if (length(partial) == 2) return(partial)
    
    # create CV marks (it doesn't matter which activity data was chose)
#     cv_mark <- get_cv_mark_mat(partial[['flags']], partial[['hitc']], partial[['modl_acc']])
#     partial[['cv_mark']] <- cv_mark
    
    # make activities matrix as 0.0001
    partial <- assign_reverse_na_number(partial, act_mat_names=act_mat_names)
    
    
    # remove inconclusive label (but keep the untested ones)
    if (noincon_label) partial <- remove_inconclusive_label(partial, act_mat_names=act_mat_names)
    
    # remove complete no hit assays (non-tested included)
    if (rm_nohit_assay) partial <- remove_nohit_assays(partial, act_mat_names=act_mat_names)

#    acts <- partial[c( act_mat_names, 'struct', 'cv_mark')]
    acts <- partial[c( act_mat_names, 'struct', 'flags')] # the flags already have different colnames to act_mat_names
    return(acts)
  })
  
#heatmap_para_generator()  
  heatmap_para_generator <- reactive({
    sort_meth <- input$sort_method
    activity_type <- ''
    
    # get all chemical information
    input_chemical_name <- NULL
    chem_id_df <- get_lookup_list(chemical_loader()[['id']], master)
    if (! is.null(chem_id_df$input_Chemical.Name)) {
      input_chemical_name <- conversion(chem_id_df, inp='GSID', out='input_Chemical.Name')
    }
    
    
    
    # the basic identifies , GSID + Cluster
    ip <- subset(chem_id_df, GSID != '' & CAS != '', select=c(GSID, Cluster))
    
    # the cleaned matrices
    dt <- matrix_editor()
    
    # get the partial and change name after changing name sort matrix
    dt <- get_input_chemical_mat(ip, dt)
    
    # if the input is data matrix, creat a blank CV matrix
    if (length(dt) == 2 ) 
    {
      activity_type <- names(dt)[1]
      act <- dt[[1]]
      cv <-  matrix("", nrow(act), ncol(act), dimnames=dimnames(act))
    } else
    {
     
      activity_type <- input$acttype
      act <- dt[[activity_type]]
      # create CV marks (it doesn't matter which activity data was chose)
      cv_mark <- get_cv_mark_mat(dt[['flags']],  act)
      #cv <- dt[['cv_mark']]
      cv <- cv_mark
      
    }
    
    # struct matrix
    struct <- dt[['struct']]
    
    # first, cluster the chemicals
    dcols <- dist(struct, method = "binary") ## chemicals
    
    # very, very cumbersome functions. better to split, merge dt + activity_type
    annotation <- get_heatmap_annotation(dcols, ip, master, input_chemical_name=input_chemical_name, dmat=dt, actType=activity_type) #data.frame output
    annt_colors <- get_heatmap_annotation_color(annotation,  actType=activity_type)
    
    # cluster compounds by various methods
    if (sort_meth == 'actclust')
    {
      dcols <- dist(act, method = "euclidean") ## chemicals by assays
      
    } else if (sort_meth == 'toxscore' )
    {
      tox_order <- rownames(annotation)[order(annotation$toxScore)]
      act <- act[tox_order, ]
      cv <- cv[tox_order, ]
    }
    # cluster assays by similarity 
    drows <- dist(t(act) , method = "euclidean") ## assays
    return(list(dcols=dcols, drows=drows, annotation=annotation, annt_colors=annt_colors, act=act, struct=struct, cv=cv))
    
  })

  
    select_plot <- reactive({
    showDendrogram <- input$showdendro
    keepsize <- input$keepsize
    sort_meth <- input$sort_method
    fsize <- input$fontsize


    color <- potency_colors
    breaks <- potency_breaks
    leg_labels <- potency_leg_labels
    leg_breaks <- potency_leg_breaks

    
    if (! is.null(chemical_loader()) )
    {
      # note pheatmap input has to have the same order!!!
      paras <- heatmap_para_generator()
      act <- paras[['act']]
      cv <- paras[['cv']]
      dcols <- paras[['dcols']]
      drows <- paras[['drows']]
      annotation <- paras[['annotation']]
      annt_colors <- paras[['annt_colors']]
      
      if (! showDendrogram)
      {
        if (sort_meth != 'toxscore')
        {
          #pheatmap v. < 1.0
          #p <- pheatmap_new_label(t(act), t(cv), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks,breaks=breaks, color=color, display_numbers=TRUE, clustering_distance_rows = drows, clustering_distance_cols = dcols,  clustering_method = "average")
          p <- pheatmap(t(act),  fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks,breaks=breaks, color=color, display_numbers=t(cv), clustering_distance_rows = drows, clustering_distance_cols = dcols,  clustering_method = "average")
        } else
        {
          #pheatmap v. < 1.0
          #p <- pheatmap_new_label(t(act), t(cv), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks, breaks=breaks, color=color, display_numbers=TRUE, clustering_distance_rows = drows, cluster_cols = FALSE, clustering_method = "average")
          p <- pheatmap(t(act),  fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks, breaks=breaks, color=color, display_numbers=t(cv), clustering_distance_rows = drows, cluster_cols = FALSE, clustering_method = "average")
        }
      } else if (sort_meth != 'toxscore' )
      {
        p <- plot(hclust(dcols, method="average"), hang=-1)
      }
    }
    return(p)
  })

  output$contents <- renderDataTable({
    if ( ! is.null(chemical_loader()) ) get_lookup_list(chemical_loader()[['id']], master)
  })

  output$dd <- renderDataTable({
    
    return(as.data.frame(matrix_editor()[['modl_acc']]))
    
#     paras <- heatmap_para_generator() #heatmap_para_generator
#     act <- paras[['act']]
#     annotation <- paras[['annotation']]
#     
#     id_info <- chemical_loader()
#     id_data <- master
#     isUpload <- FALSE
#     if(length(id_info) > 1) {
#       id_data <- id_info[['id']]
#       isUpload <- TRUE
#     }
#     result <- get_output_df(act, annotation, id_data, isUpload)
#     
#     return(result)
    
    # for testing
#       paras <- heatmap_para_generator()
       #return(paras[['struct']][,1:10])
#       return(data.frame(rownames(paras[['act']])))
       #return(data.frame(paras[['annotation']]))
    
  })
  
  output$assay_info <- renderDataTable({
  
    col_n <- c('assay_source_name','timepoint_hr','organism','tissue','cell_format','cell_short_name',
               'assay_format_type', 'assay_component_target_desc', 'assay_design_type', 
               'biological_process_target','detection_technology', 'assay_component_endpoint_name',
               'normalized_data_type','burst_assay','key_positive_control','signal_direction',
               'intended_target_type','intended_target_family', 'assay_component_endpoint_desc',
               'intended_target_entrez_gene_id','intended_target_gene_name','intended_target_gene_symbol')
    result <- assay_names[, colnames(assay_names) %in% col_n]
    result <- select(result, assay_source_name, assay_component_endpoint_name, signal_direction, key_positive_control,
                     intended_target_family, intended_target_type, 
                     intended_target_gene_symbol, intended_target_gene_name,intended_target_entrez_gene_id,
                     biological_process_target,detection_technology, assay_design_type,normalized_data_type,
                     timepoint_hr,organism,tissue,cell_format,cell_short_name,burst_assay,assay_component_target_desc,
                     assay_component_endpoint_desc
                     )
    #result <- assay_names
    return(result)
    
  })

  getVarWidth <- reactive({
    ncmpd <- 0
    keepsize <- input$keepsize
    if ( ! is.null(chemical_loader()) & ! keepsize) 
    {
      chem_id_df <- get_lookup_list(chemical_loader()[['id']], master)
      ip <- subset(chem_id_df, GSID != '' & CAS != '', select=c(GSID, Cluster))
      ncmpd <- nrow(ip)
    }
    if (ncmpd < 40)
    {
      return(1200)
    } else
    {
      return(ncmpd*30)
    }
  })
  
  getVarHeight <- reactive({
    nassay <- 0
    keepsize <- input$keepsize
    if ( ! is.null(heatmap_para_generator()) & ! keepsize) 
    {
      paras <- heatmap_para_generator()
      nassay <- ncol(paras$act)
    }
    if (nassay < 40)
    {
      return(1000)
    } else
    {
      return(nassay*30)
    }
  })
  
  output$profiling <- renderPlot({
    select_plot()
 
  }, width=getVarWidth, height=getVarHeight)
  
  output$box <- renderPlot({
    
    fsize <- input$fontsize
    sort_meth <- input$sort_method
    
    p <- NULL
    
      activity_type <- input$acttype
     
        paras <- heatmap_para_generator()
        act <- paras[['act']]
        annotation <- paras[['annotation']]
        dcols <- paras[['dcols']]
        
        id_info <- chemical_loader()
        id_data <- master
        isUpload <- FALSE
        if(length(id_info) > 1) {
          id_data <- id_info[['id']]
          isUpload <- TRUE
        }
        result <- get_output_df(act, annotation, id_data, isUpload)
        
        p <- get_pod_boxplot(result, fontsize=fsize, sortby=sort_meth, dcols=dcols, global_para=assay_names)
    
    
    if (! is.null(p)) print(p)
    
  },  width=getVarWidth, height=getVarHeight)

  output$downloadData <-  downloadHandler(
    filename = function() {
      
        paste( input$acttype, '.txt', sep='')
      
       },
    content = function(file) {
      paras <- heatmap_para_generator()
      act <- paras[['act']]
      annotation <- paras[['annotation']]
      result <- get_output_df(act, annotation, master)
      #result <- get_published_data_only_commonname(result, assay_names)  # to remove unpublished data
      write.table(result, file, row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE, append=FALSE)
    }
  )

output$downloadPlot <- downloadHandler(
       filename = function() { 
        
           paste( input$acttype, '.pdf', sep='')
        
         },
       content = function(file) {
         #png(file, width=9, height=6.5, units="in", res=600)
         pdf(file, width=9, height=6.5)
         select_plot2()
         dev.off()
       }
)

select_plot2 <- function () {
  showDendrogram <- input$showdendro
  keepsize <- input$keepsize
  sort_meth <- input$sort_method
  fsize <- input$fontsize

      color <- potency_colors
      breaks <- potency_breaks
      leg_labels <- potency_leg_labels
      leg_breaks <- potency_leg_breaks
  
  
  if (! is.null(chemical_loader()) )
  {
    # note pheatmap input has to have the same order!!!
    paras <- heatmap_para_generator()
    act <- paras[['act']]
    cv <- paras[['cv']]
    dcols <- paras[['dcols']]
    drows <- paras[['drows']]
    annotation <- paras[['annotation']]
    annt_colors <- paras[['annt_colors']]
    
    if (! showDendrogram)
    {
     if (sort_meth != 'toxscore')
      {
        #p <- pheatmap_new_label(t(act), t(cv), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks,breaks=breaks, color=color, display_numbers=TRUE, clustering_distance_rows = drows, clustering_distance_cols = dcols,  clustering_method = "average")
        p <- pheatmap(t(act),  fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks,breaks=breaks, color=color, display_numbers=t(cv), clustering_distance_rows = drows, clustering_distance_cols = dcols,  clustering_method = "average")
      } else
      {
        #p <- pheatmap_new_label(t(act), t(cv), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks, breaks=breaks, color=color, display_numbers=TRUE, clustering_distance_rows = drows, cluster_cols = FALSE, clustering_method = "average")
        p <- pheatmap(t(act), fontsize=fsize,annotation=annotation,annotation_colors=annt_colors,legend_labels=leg_labels,legend_breaks=leg_breaks, breaks=breaks, color=color, display_numbers=t(cv), clustering_distance_rows = drows, cluster_cols = FALSE, clustering_method = "average")
      }
    } else if (sort_meth != 'toxscore' )
    {
      p <- plot(hclust(dcols, method="average"), hang=-1)
    }
  }
  return(p)
}
  

    
  
  
  
})
