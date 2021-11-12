#library("sf")
library(plotly)
library(shinydashboard)
library(DESeq2)
library(reshape)
de_sig_res <- read.csv('/home/sohn/projects/transcriptomics_browser/MEND_AD_vs_CO/de_res_sig.csv', header=T, stringsAsFactors = F, check.names = F, sep=",")
modelName <- c("Sex+Age+Status", "Sex+Age+PMI+Status", "Sex+Age+PMI+Cohort+Status","Sex+Age+Cohort+Status",
               "Sex+Age+Braak", "Sex+Age+PMI+Braak", "Sex+Age+PMI+Cohort+Braak", "Sex+Age+Cohort+Braak",
               "Sex+Age+CDR", "Sex+Age+PMI+CDR", "Sex+Age+PMI+Cohort+CDR", "Sex+Age+Cohort+CDR")
function(request) {
  sidebar <- dashboardSidebar(width = 240,
    hr(),
    tags$head(tags$style(HTML('.logo {
                              background-color: #a51417 !important;
                              }
                              .navbar {
                              background-color: #a51417 !important;
                              }
                              .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              border-left-color: #a51417;
                              }
                              .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover {
                              border-left-color: #a51417;
                              }  
                              .box.box-solid.box-primary>.box-header {
                              color:#000000;
                              background:#c8c8c8
                              }
                              
                             .nav-tabs-custom .nav-tabs li.active:hover a, .nav-tabs-custom .nav-tabs li.active a {
                              background-color: transparent;
                              border-color:  #000000;
                              }

                              .nav-tabs-custom .nav-tabs li.active {
                              border-top-color: #FFF;
                              }
                              
                              .col-sm-3 {
                                /* width: 25%; */
                              }
                              
                              .col-sm-4 {
                                /* width: 33.33333333%; */
                              }
                              
                              .col-sm-6 {
                                /* width: 25%; */
                              }
                              .col-sm-12 {
                                /* width: 100%; */
                              }
                              
                              .box.box-solid.box-primary{
                              border-bottom-color:#c8c8c8;
                              border-left-color:#c8c8c8;
                              border-right-color:#c8c8c8;
                              border-top-color:#c8c8c8;
                              }'))),
    
    sidebarMenu(id="tabs",
                menuItem("Home", icon = icon("home"), tabName = "home", selected=TRUE),
                #menuItem("Studies", tabName="plot", icon=icon("list")),
                #menuItem("Table import", tabName = "table", icon=icon("table")),
                menuItem("Studies",  icon = icon("list"),
                      menuItem("Differential Expresssion", icon = icon("angle-right"),
                      menuItem("AD_vs_CO", tabName="Mend_AD_vs_CO", icon = icon("angle-right")
                                          # menuSubItem("Model 3: Sex + Age + BraakTau + Status", tabName = "SexAgeBrack"),
                                          # menuSubItem("Model 4: Sex + Age + PMI + Status", tabName = "SexAgePMIStatus"),
                                          # menuSubItem("Model 5: Sex + Age + CDR", tabName = "SexAgeCDR")
                                           #menuSubItem("ADAD_vs_CO", tabName="Mend_ADAD_vs_CO"),
                                           #menuSubItem("TREM2_vs_CO", tabName="Mend_TREM2_vs_CO"),
                                           #menuSubItem("FTD_vs_CO", tabName="Mend_FTD_vs_CO"),
                                           #menuSubItem("GRN_vs_CO", tabName="Mend_GRN_vs_CO"),
                                           #menuSubItem("C9orf72_vs_CO", tabName="Mend_C9orf72_vs_CO"),
                                           #menuSubItem("MAPT_vs_CO", tabName="Mend_MAPT_vs_CO")
                                         )
                                  
                                  #menuItem("TIN", icon = icon("angle-right"), tabName = "TIN"),
                                  #menuItem("circRNA", icon = icon("angle-right"), tabName = "circ")
                         )
                         #menuSubItem("SUNSHINE", icon = icon("list")),
                         #menuSubItem("MSBB", icon = icon("list")),
                         #menuSubItem("Mayo", icon = icon("list")),
                         #menuSubItem("CommonMind", icon = icon("list"))
                                     
                ),
                #menuItem("ReadMe", tabName = "readme", icon=icon("mortar-board")),
                menuItem("About", tabName = "about", icon = icon("question"))
    )
    
  )
 
  
## set the body structure
body <- dashboardBody(
    tabItems(
        tabItem(tabName = "home",
                fluidPage(
                  tags$iframe(src = './firstPgae.html', 
                              width = '100%', height = '800px',
                              frameborder = 0, scrolling = 'auto'
                  )
                )
        ),
      
        tabItem(tabName = "MenDesc",
                fluidPage(
                  tags$iframe(src = './MenDesc.html', 
                              width = '100%', height = '800px',
                              frameborder = 0, scrolling = 'auto'
                  )
                )
       ),
      tabItem(tabName = "Mend_AD_vs_CO",
              fluidRow(
                box(width = 8, collapsible = FALSE,
                    fluidRow(
                      column(width = 3,
                             selectInput("modelDE", "Choose DE Model", choices = modelName, selected = "Sex_Age_Status"),
                      )
                    ),
                    column(width = 3,
                           selectInput("selectedGeneBiotype", "Choose gene biotype", choices = c("All", unique(de_sig_res$GeneBiotype)), selected = "All"),
                    ),
                    column(width = 3,
                           sliderInput("FCThresholdSlider", "Select a fold change threshold",  min = 0, max = 10, value = NA),
                           textInput("FCThresholdText", "", value=0),
                           useShinyalert(),
                    ),
                    column(width = 4,
                           sliderInput("pvalSlider", "Select pvalue threshold",  min = 0, max = 0.10, value = NA),
                           textInput("pvalText", "", value=0),
                    ),
                    title = "Thresholds", solidHeader = TRUE, status = "primary"
                )),
              fluidRow(
                width = 6,
                box(width = 8, collapsible = FALSE, 
                    #sampleModuleUI("AD_vs_CO"),
                    DT::dataTableOutput('mytable', width = 700),
                    title = "Differential Expression Data", solidHeader = TRUE, status = "primary"
                )
              ),
              fluidRow(
                box(width = 4, collapsible = FALSE,
                    plotOutput('gene_plot', width = 500),
                    title = "Gene Expression", solidHeader = TRUE, status = "primary"
                ),
                box(width = 4, collapsible = FALSE,
                    plotOutput("vol_plot", width = 340), 
                    #plotVolModuleUI("vol_plot"),
                    title = "Volcano Plot", solidHeader = TRUE, status = "primary"),
                
              )
     )
     # tabItem(tabName = "SexAgeStatus",
     #         fluidRow(
     #           box(width = 8, collapsible = FALSE,
     #               column(width = 3,
     #                      selectInput("selectedGeneBiotype1", "Choose gene biotype", choices = c("All", unique(de_sig_res$GeneBiotype)), selected = "All"),
     #               ),
     #               column(width = 3,
     #                      sliderInput("FCThreshold1", "Select a fold change threshold",  min = 0, max = 10, value = 0),
     #               ),
     #               column(width = 4,
     #                      sliderInput("SexAgePMIpvalSlider", "Select pvalue threshold",  min = 0, max = 0.10, value = 0.10),
     #                      textInput("SexAgeStatuspvalText", "", value=0.1),
     #               ),
     #               title = "Thresholds", solidHeader = TRUE, status = "primary"
     #           )),
     #         fluidRow(
     #           width = 6,
     #           box(width = 8, collapsible = FALSE, 
     #               #sampleModuleUI("AD_vs_CO"),
     #               DT::dataTableOutput('mytable1', width = 700),
     #               title = "Differential Expression Data", solidHeader = TRUE, status = "primary"
     #           )
     #         ),
     #         fluidRow(
     #           box(width = 4, collapsible = FALSE,
     #               plotOutput('gene_plot1', width = 280),
     #               title = "Gene Expression", solidHeader = TRUE, status = "primary"
     #           ),
     #           box(width = 4, collapsible = FALSE,
     #               plotOutput("vol_plot1", width = 340), 
     #               #plotVolModuleUI("vol_plot"),
     #               title = "Volcano Plot", solidHeader = TRUE, status = "primary"),
     #           
     #         )
     # )
      
    ) 
  )

## put the dashboard and sidebar together 
dashboardPage(
    dashboardHeader(title = "Transcriptomics Browser"),
    sidebar,
    body
)
  
  
}
  
