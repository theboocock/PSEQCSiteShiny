
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggvis)



shinyUI(fluidPage(
  
  # Application title
  titlePanel("Site QC Resquencing"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      #   uiOutput("fileInput"),
      conditionalPanel(condition = "input.conditionedPanels == 'QC Table'"
                       ,checkboxGroupInput('show_vars',"Columns in Site QC table to show:",
                                           names(qc_table),selected=names(qc_table))
      ),
      conditionalPanel(condition ="input.conditionedPanels == 'QC Plots Urate'",
                       numericInput("hwe_thresh",label="HWE threshold for plots",min=0,max=1,step=0.025,value=0.01),
                       checkboxGroupInput("filters","Filters to exclude: ", unique(as.character(qc_table$filter))),
                      numericInput("gqm",min=0,max=99,label="GQM -Min genotype Mean for Cases and Controls",value=0),
                      
                       numericInput("max_prop_ref",min=0,max=1,step=0.025,value=1,label="SAMPLE_REF_AB - Maximum mean proportion of alternate alleles in homozygous reference call"),
                       numericInput("min_prop_alt",min=0,max=1,step=0.025,value=0,label="SAMPLE_HOM_AB - minimum mean proportion of alternate alleles in homozygous alternate call "),
                       numericInput("low_het_ab_prop",min=0,max=0.5,step=0.025,value=0,
                                    label="Min SAMPLE_HET_AB - Minimum mean AB ratio for case and control heterozygotes"),
                       numericInput("high_het_ab_prop",min=0.5,max=1.0,step=0.025,value=1,
                                    label="Max SAMPLE_HET_AB - Maximum mean AB ration for case and control heterozygotes"),
                       numericInput("prop_20_80",label="PROP_DEVIANT_20_80 - Maximum proportion of a site with A/B ration between 0-0.2 and 0.8-1.0",min=0,
                                    max=1,step=0.025,value=1),
                       numericInput("prop_30_70",label="PROP_DEVIANT_30_70 - Maximum proportion of a site with A/B ratio between 0-0.3 and 0.7-1.0",
                                    min = 0, max=1,step=0.025, value =1),
                       numericInput("mean_difference",label="Maximum absolute value of the differences between casen difference",
                                    min = 0, max = 200, step = 1, value = 200),
                       numericInput("number_samples_gq_10",label="Maximum number of sites with genotype quality <= 10",
                                    min = 0,max=1000,step=1,value=0),
                       downloadButton("downloadSites","Download Final Sites")
      ),
      
      conditionalPanel(condition ="input.conditionedPanels == 'All Site QC'",
                       numericInput("n_gqm",min=0,max=99,label="Min genotype Mean",value=0),
                       numericInput("n_low_het_ab_prop",min=0,max=0.5,step=0.025,value=0,
                                    label="Minimum mean AB ratio for heterozygotes"),
                       numericInput("n_high_het_ab_prop",min=0.5,max=1.0,step=0.025,value=1,
                                    label="Maximum mean AB ration for heterozygotes"),
                       numericInput("n_prop_20_80",label="Maximum proportion of a site with A/B ration between 0-0.2 and 0.8-1.0",min=0,
                                    max=1,step=0.025,value=1),
                       numericInput("n_prop_30_70",label="Maximum proportion of a site with A/B ratio between 0-0.3 and 0.7-1.0",
                                    min = 0, max=1,step=0.025, value =1)
      )
    ),
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        #  tabPanel("Load Data", uiOutput("load_data")),
        tabPanel('QC Table', dataTableOutput("qc_table")),
        tabPanel('QC Plots Urate',
                 tableOutput("variant_summary"),
                 tableOutput("ts_tv_table"),
                 tableOutput("singleton_table"),
                 h4("Percentarge of sites passing HWE >= 0.01 filter"),
                 tableOutput('hwe_table'),
                 fluidRow(
                   column(6,
                          plotOutput('hwe_cases_plot')),
                   column(6,
                          plotOutput('hwe_controls_plot'))),
                 plotOutput("allelic_balance_controls"),
                 plotOutput("allelic_balance_cases"),
                 plotOutput("overall_allelic_balance"),
                 plotOutput('qqplot'),
                 h4("Linear model containing predictors for the pvalues"),
                 tableOutput('linear_model')
        ),
        
        tabPanel('All Site QC',
                 tableOutput("variant_summary2"),
                 tableOutput("ts_tv_table2"),
                 plotOutput("hwe_plot")),
        id="conditionedPanels")
    )
  )
)
)

