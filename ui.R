library(shiny)
library(markdown)


shinyUI(pageWithSidebar(
  
  # application title
  headerPanel("ToxCast Activity Profiler"),
  
  sidebarPanel(
    
    tags$head(tags$meta(`http-equiv`="pragma", content="no-cache"), 
              tags$meta(`http-equiv`="Cache-control", content="no-cache, no-store")),
    
    # profiling options
    #h4('Profiling'),
    
    radioButtons("acttype", "Activity type",
          list("ACC" = "modl_acc",  
               "ACB" = "modl_acb", 
               "AC50" = "modl_ga", 
               "AC10"= 'modl_ac10'), selected="modl_acc"),
 
    tags$hr(),
    
    h4('Activity filtering'),
      
    tags$br(),
    sliderInput("sresp_thres", 
                  "scaled resp threshold", min=1, max=10, value=1, step=1),
    tags$br(),
    #sliderInput("acc_thres", 
    #              "ACC threshold",  min=3, max=10, value=3, step=0.5),
    numericInput("acc_thres", label = "ACC threshold (uM)", value = NA),
    
    tags$br(),
    #sliderInput("acb_thres", 
    #            "ACB threshold",  min=3, max=10, value=3, step=0.5),
    numericInput("acb_thres", label = "ACB threshold (uM)", value = NA),
    tags$br(),
    #sliderInput("ga_thres", 
    #              "AC50 threshold", min=3, max=10, value=3, step=0.5),
    numericInput("ga_thres", label = "AC50 threshold (uM)", value = NA),
    tags$br(),
    #sliderInput("ac10_thres", 
    #            "AC10 threshold", min=3, max=10, value=3, step=0.5),
    numericInput("ac10_thres", label = "AC10 threshold (uM)", value = NA),
    #tags$br(),
    #checkboxInput("hitfilter", "apply hit call filter", TRUE),
    
    tags$br(),
    checkboxInput("cytofilter", "apply cytotoxicity filter", FALSE),
      
    tags$br(),
    checkboxInput("flagfilter", "apply flag (from curve fitting) filter", FALSE),
    
    tags$br(),
    checkboxInput("noinconlab", "make inconclusive (compounds with cytotoxic or flag label) as inactive", TRUE),
      
    tags$hr(),
    
    # arrange the compounds
    h4('Sort compounds by ...'),
    
    tags$br(),
    radioButtons("sort_method", "Method:",
                 list('structure similarity' = 'chemclust',
                      'activity similarity' = 'actclust',
                      'toxicity score (only for activity)' = 'toxscore')),
    
    tags$hr(),
    
    # filter the assays
    h4('Filter assays'),
    
    tags$br(),
  
    checkboxInput("rm_nohit_assay", "remove no-hit (non-tested included) assays", TRUE),
  
    tags$br(),
    h5('by assay source'),
    selectizeInput('assay_source', 'sources', choices = list(
      Primary = c('ACEA','APR','ATG','BSK','NVS','OT', 'CLD'),
      Others = c('CEETOX','TANGUAY','NHEERL_PADILLA','NCCT_SIMMONS', 'TOX21')
      ), 
#      options = list(
#        placeholder = 'Please select multiple options below',
#        onInitialize = I('function() { this.setValue(""); }')),
      multiple = TRUE,
      selected=c('ACEA','APR','ATG','BSK','NVS','OT', 'CLD')
      ),
    
    tags$br(),
    h5('by assay name'),
    textInput('reg_assay', 'assays (regular expression)', ''),
    checkboxInput("inv_sel_assay", "invert your selection", FALSE),
  
    tags$br(),
    h5('by gene'),
    textInput('reg_gene', 'genes (case insensitive)', ''),
    #checkboxInput("inv_sel_gene", "invert your selection", FALSE),
    
    tags$hr(),
    
    # input control
    h4('Input'),
    tags$textarea(id="cmpds", rows=3, cols=1, ""),
    helpText("Note: copy & paste from excel file with two columns: CAS & Cluster"),
    
    tags$hr(),
    
    # miscellaneous functions
    h4('Others'),

    checkboxInput("showdendro", "show compound similarity dendrogram", FALSE),
    checkboxInput("keepsize", "keep heatmap size one page", FALSE),
    
    tags$br(),
    
    # fontsize
    sliderInput("fontsize", 
                "fontsize", min = 2,max = 28, value = 16, step=2),
    helpText("tip: lower the fontsize when saving the plot"),
  
    # output functions
    br(),
    downloadButton('downloadData', 'Download Activities'),
    downloadButton('downloadPlot', 'Save Plot'),
    downloadButton('downloadEnrich', 'Download Enrichment analysis')
    
  ),
  mainPanel(
    tabsetPanel(
     
      tabPanel( 'Input chemicals', dataTableOutput('contents')),
      tabPanel( "Profile", plotOutput("profiling", height=1000, width="500%")), # i think the height don't affect 1000
      #tabPanel( "Potency boxplot", plotOutput("box",  height=1000, width="500%")),
      tabPanel( 'Activity data', dataTableOutput('dd')),
      tabPanel('Enrichment analysis', dataTableOutput('enrich')),
      tabPanel( 'Assays', dataTableOutput('assay_info')),
      tabPanel('About', includeHTML("README.html"))
    )
  )
))

