# Load required libraries
library(shiny)
library(bs4Dash)
library(ggrepel)
library(geomtextpath)
library(shinyWidgets)
library(DT)
library(fresh)

# setwd("/Users/Cherry/2025/HDS CDT")
#setwd("/Users/orayasrim/Documents/GRAM_shiny_workshop/economic_evaluation_tutorial-main")

# Source the external script containing the decision tree model
source("practical_cea_decisiontree.R")


# Define UI for the application
ui <- dashboardPage(
  freshTheme = create_theme(
    bs4dash_status(primary = "#0077b6", danger = "#BF616A")),
  dark = FALSE,
  title = "CEA",
  help = NULL,
  header = dashboardHeader(
    title = dashboardBrand(
      title = "Cost-Effectiveness Analysis",
      color = "primary"
      #color = "primary"
    ),
    skin = "light"
  ),
  sidebar = dashboardSidebar(
    skin = "light",
    status = "primary",
    collapsed = FALSE,
    minified = TRUE,
    sidebarMenu(
      id = "sidebar",
      menuItem(
        "Dashboard",
        tabName = "dashboard",
        icon = icon("dashboard")
      ),
      menuItem(
        "Infection Parameters",
        tabName = "infection_params",
        icon = icon("virus")
      ),
      # menuItem(
      #   "Culture Positive Parameters",
      #   tabName = "culture_positive",
      #   icon = icon("flask")
      # ),
      menuItem(
        "Culture Negative Parameters",
        tabName = "culture_negative",
        icon = icon("flask")
      ),
      menuItem(
        "No Culture Parameters",
        tabName = "culture_no",
        icon = icon("flask")
      ),
      menuItem(
        "Economic Parameters",
        tabName = "economic_params",
        icon = icon("money-bill")
      ),
      menuItem(
        "About",
        tabName = "about",
        icon = icon("info-circle")
      ),
      sliderInput("N", "Patient Population", value = 1000, min = 100, max = 10000, step = 100),
      sliderInput("will_pay_perDALY", 
                  label = "Willingness to Pay Per DALY averted ($/DALY)",
                  min = 0,
                  max = 5000,
                  value = 1500,
                  step = 20)
    )
  ),
  body = dashboardBody(
    tabItems(
      # Dashboard Tab
      tabItem(
        tabName = "dashboard",
        fluidRow(
          column(
            width = 12,
            bs4Card(
              title = "Overview",
              status = "primary",
              width = 12,
              solidHeader = TRUE,
              collapsible = TRUE,
              closable = FALSE,
              p("This application evaluates the cost-effectiveness of implementing an intervention (example with maintaining an active microbiology laboratory for blood culture testing in healthcare settings). The model uses a decision tree approach to compare strategies with and without active microbiology services."),
              p("Use the sidebar and the 'Bloodstream Infections Parameters' (BSI) tab below to navigate through the different parameter input sections, then run the model by clicking the 'Compute' button to see the results."),
              p("Navigate to the 'Kenya Example: Bloodstream Infections Parameters' tab below to set the proportion of BSI each top pathogen contributes to BSIs and the proportion each pathogen is susceptible to 1st line treatment. "),
              # actionButton("run", "Run Analysis", class = "btn-lg btn-success")
              actionButton(
                inputId = "run",
                class = "btn-primary",
                width = "300px", style = "color: white; background-color: #228B22",
                "Compute"
              )
            )
          )
        ),
        fluidRow(
          column(
            width = 12,
            bs4TabCard(
              title = "",
              elevation = 2,
              status = "primary",
              width = 12,
              solidHeader = FALSE,
              collapsible = TRUE,
              maximizable = TRUE,
              tabPanel(
                title = "Cost-Effectiveness Plane",
                br(),br(),
                
                "The example shown on the cost-effectivness plane is comparing having an microbiology service system vs no microbiology service system. This incremental cost-effectiveness plane is showing 
                disability-adjusted life years (DALYs) averted versus incremental costs per value of hospitalised adult patients (default value is set to 1000) from maintaining a microbiology laboratory compared to no microbiological testing 
                when ceftriaxone and gentamicin is used as empirical treatment. In the upper right quadrant, if the point is below the 
                dotted diagonal line, then the active microbiology service would be considered cost-effective for a willingness to pay to 
                avert one DALY set on the left tab (default value set to $1500/DALY Averted), while any point above this line would not be considered 
                cost-effective at this threshold. You can adjust the different parameters in the tabs on the left to see how this changes the cost-effectivness of this example.",
                br(),br(),
                "Adjust the limit of the plots using the sliders below.",
                br(),br(),
                sliderInput(inputId = "SLIDERX", label = "X-axis limit", min = 0, max = 10000, value = 2000, step = 100 ),
                sliderInput(inputId = "SLIDERY", label = "Y-axis limit", min = 0, max = 100000000, value = 8000000, step = 1000),
                
                br(),br(),
                plotOutput("plotCEP", height = "600px")
              ),
              tabPanel(
                title = "Kenya Example: Bloodstream Infections Parameters",
                br(),
                "In this example, we will look at the top 7 pathogens causing bloodstream infections (BSI) in Kenya. The top pathogens from are Klebsiella pneumoniae, Staphylococcus aureus, Pseudomonas aeruginosa,
                Streptococcus pneumoniae, Neisseria meningitidis, Acinetobacter baumannii, and Escherichia coli. Adjust the sliders to specify the proportion
                of BSI they cause and the proportion suseptible to recommended first line treatment.", br(),br(),

                
                bs4Card(
                  title = "Klebsiella pneumoniae",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  fluidRow(
                    column(
                      width = 6,
                      sliderInput("kleb_prop_bsi", "Proportion of suspected BSI caused by this pathogen", value = 0.05, min = 0, max = 0.10, step = 0.01)#,
                    ) ,
                    column(
                      width = 6,
                       sliderInput("kleb_prop_line", "Proportion Suseptible to Recommended 1st line treatment", value = 0.95, min = 0, max = 1, step = 0.05)
                    )
                  )
                ),
                bs4Card(
                  title = "Staphylococcus aureus",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  fluidRow(
                    column(
                      width = 6,
                      sliderInput("SA_prop_bsi", "Proportion of suspected BSI caused by this pathogen", value = 0.05, min = 0, max = 0.10, step = 0.01)#,
                    ) ,
                    column(
                      width = 6,
                      sliderInput("SA_prop_line", "Proportion Suseptible to Recommended 1st line treatment", value = 0.95, min = 0, max = 1, step = 0.05)
                    )
                  )
                ),
                bs4Card(
                  title = "Pseudomonas aeruginosa",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  fluidRow(
                    column(
                      width = 6,
                      sliderInput("PA_prop_bsi", "Proportion of suspected BSI caused by this pathogen", value = 0.05, min = 0, max = 0.10, step = 0.01)#,
                    ) ,
                    column(
                      width = 6,
                      sliderInput("PA_prop_line", "Proportion Suseptible to Recommended 1st line treatment", value = 0.95, min = 0, max = 1, step = 0.05)
                    )
                  )
                ),
                bs4Card(
                  title = "Streptococcus pneumoniae",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  fluidRow(
                    column(
                      width = 6,
                      sliderInput("SP_prop_bsi", "Proportion of suspected BSI caused by this pathogen", value = 0.05, min = 0, max = 0.10, step = 0.01)#,
                    ) ,
                    column(
                      width = 6,
                      sliderInput("SP_prop_line", "Proportion Suseptible to Recommended 1st line treatment", value = 0.95, min = 0, max = 1, step = 0.05)
                    )
                  )
                ),
                bs4Card(
                  title = "Neisseria meningitidis",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  fluidRow(
                    column(
                      width = 6,
                      sliderInput("NM_prop_bsi", "Proportion of suspected BSI caused by this pathogen", value = 0.05, min = 0, max = 0.10, step = 0.01)#,
                    ) ,
                    column(
                      width = 6,
                      sliderInput("NM_prop_line", "Proportion Suseptible to Recommended 1st line treatment", value = 0.95, min = 0, max = 1, step = 0.05)
                    )
                  )
                ),
                bs4Card(
                  title = "Acinetobacter baumannii",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  fluidRow(
                    column(
                      width = 6,
                      sliderInput("AB_prop_bsi", "Proportion of suspected BSI caused by this pathogen", value = 0.05, min = 0, max = 0.10, step = 0.01)#,
                    ) ,
                    column(
                      width = 6,
                      sliderInput("AB_prop_line", "Proportion Suseptible to Recommended 1st line treatment", value = 0.95, min = 0, max = 1, step = 0.05)
                    )
                  )
                ),
                bs4Card(
                  title = "Escherichia coli",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  fluidRow(
                    column(
                      width = 6,
                      sliderInput("EC_prop_bsi", "Proportion of suspected BSI caused by this pathogen", value = 0.05, min = 0, max = 0.10, step = 0.01)#,
                    ) ,
                    column(
                      width = 6,
                      sliderInput("EC_prop_line", "Proportion Suseptible to Recommended 1st line treatment", value = 0.95, min = 0, max = 1, step = 0.05)
                    )
                  )
                )
              )
            )
          )
        )
      ),
      
      # Infection Parameters Tab
      tabItem(
        tabName = "infection_params",
        fluidRow(
          column(
            width = 12,
            bs4Card(
              title = "Probability of Infection",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              fluidRow(
                column(
                  width = 6,
                  numericInput("p.bact", "Probability of Bacterial Infection", value = 0.5, min = 0, max = 1, step = 0.05)#,
                  #valueBoxOutput("prob_no_bact", width = 12)
                ) #,
                # column(
                #   width = 6,
                #    numericInput("N", "Population Size", value = 1000, min = 100, max = 10000, step = 100)
                # )
              )
            )
          )
        ),
        fluidRow(
          column(
            width = 12,
            bs4Card(
              title = "Culture Results",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              fluidRow(
                column(
                  width = 6,
                  numericInput("p.bact_culture_pos", "Probability Culture Positive when Bacterial Infection Present", value = 0.5, min = 0, max = 1, step = 0.05)
                ),
                column(
                  width = 6,
                  numericInput("p.nobact_culture_pos", "Probability Culture Positive when No Bacterial Infection", value = 0.5, min = 0, max = 1, step = 0.05)
                )
              )
            )
          )
        )
      ),
      
      # # Culture Parameters Tab
      # tabItem(
      #   tabName = "culture_positive",
      #   fluidRow(
      #     column(
      #       width = 12,
      #       bs4Card(
      #         title = "Culture Postive Branch Parameters",
      #         status = "primary",
      #         solidHeader = TRUE,
      #         width = 12,
      #         collapsible = TRUE,
      #         fluidRow(
      #           column(
      #             width = 6,
      #             numericInput("p.CefGen_culpos_stepup", "Probability of stepping up from ceftriaxone/gentamicin patient with culture positive result", value = 0.5, min = 0, max = 1, step = 0.05)
      #           ),
      #           column(
      #             width = 6,
      #             valueBoxOutput("prob_culpos_nochange", width = 12)
      #           )
      #         )
      #       )
      #     )
      #   )
      # ),
      
      # Clinical Response Parameters Tab
      tabItem(
        tabName = "culture_negative",
        fluidRow(
          bs4Card(
            title = "Culture Negatives Branch Parameters",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            numericInput("p.CefGen_bact_culneg_deter", "Probability of clinical deterioration", value = 0.5, min = 0, max = 1, step = 0.05),
            numericInput("p.CefGen_bact_culneg_impro", "Probability of clinical improvement", value = 0.5, min = 0, max = 1, step = 0.05),
            valueBoxOutput("prob_bact_culneg_nocha", width = 12)
          ),
          bs4Card(
            title = "Culture Negative with No Bacterial Infection",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            numericInput("p.CefGen_nobact_culneg_deter", "Probability of clinical deterioration", value = 0.5, min = 0, max = 1, step = 0.05),
            numericInput("p.CefGen_nobact_culneg_impro", "Probability of clinical improvement", value = 0.5, min = 0, max = 1, step = 0.05),
            valueBoxOutput("prob_nobact_culneg_nocha", width = 12)
          )
        ),
        fluidRow(
          bs4Card(
            title = "Culture Negative Treatment Decisions",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            column(
              width = 4,
              numericInput("p.CefGen_culneg_deter_stepup_nochange", "Probability of stepping up (deterioration)", value = 0.5, min = 0, max = 1, step = 0.05)
            ),
            column(
              width = 4,
              numericInput("p.CefGen_culneg_impro_stepdown_nochange", "Probability of stepping down (improvement)", value = 0.5, min = 0, max = 1, step = 0.05)
            ),
            column(
              width = 4,
              numericInput("p.CefGen_culneg_nocha_stepup_nochange", "Probability of stepping up (no change)", value = 0.5, min = 0, max = 1, step = 0.05)
            )
          )
        ),
        fluidRow(
          column(
            width = 12,
            bs4Card(
              title = "Treatment Modification Probabilities",
              status = "info",
              solidHeader = TRUE,
              width = 12,
              collapsible = TRUE,
              column(
                width = 12,
                numericInput("p.culneg_deter_stepup_nochange", "Probability of stepping up from the current antibiotic treatment when clincial condition deteriorates", value = 0.5, min = 0, max = 1, step = 0.05),
                numericInput("p.culneg_impro_stepdown_nochange", "Probability of stepping down from the current antibiotic treatment when clincial condition improves", value = 0.5, min = 0, max = 1, step = 0.05),
                numericInput("p.culneg_nocha_stepup_nochange", "Probability of stepping up from the current antibiotic treatment when clincial condition remains unchanged", value = 0.5, min = 0, max = 1, step = 0.05)
              )
            )
          )
        )
      ),
  tabItem(
        tabName = "culture_no",
        fluidRow(
          # bs4Card(
          #   title = "No Culture with Bacterial Infection",
          #   status = "primary",
          #   solidHeader = TRUE,
          #   width = 6,
          #   collapsible = TRUE
          #   #numericInput("p.CefGen_bact_nocul_deter", "Probability of clinical deterioration", value = 0.5, min = 0, max = 1, step = 0.05),
          #   #numericInput("p.CefGen_bact_nocul_impro", "Probability of clinical improvement", value = 0.5, min = 0, max = 1, step = 0.05),
          #   #valueBoxOutput("prob_bact_nocul_nocha", width = 12)
          # ),
      bs4Card(
        title = "No Culture with No Bacterial Infection",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        collapsible = TRUE,
        numericInput("p.CefGen_nobact_nocul_deter", "Probability of clinical deterioration", value = 0.5, min = 0, max = 1, step = 0.05),
        numericInput("p.CefGen_nobact_nocul_impro", "Probability of clinical improvement", value = 0.5, min = 0, max = 1, step = 0.05),
        valueBoxOutput("prob_nobact_nocul_nocha", width = 12)
      )),
      fluidRow(
        bs4Card(
          title = "No Culture Treatment Decisions",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          column(
            width = 4,
            numericInput("p.CefGen_nocul_deter_stepup_nochange", "Probability of patient stepping up with non-bacterial infection, whose condition deteriorates under ceftriaxone/gentamicin treatment", value = 0.012, min = 0, max = 1, step = 0.01)
          ),
          column(
            width = 4,
            numericInput("p.CefGen_nocul_impro_stepdown_nochange", "Probability of stepping down with non-bacterial infection, whose condition improves under ceftriaxone/gentamicin treatment", value = 0.013, min = 0, max = 1, step = 0.01)
          ),
          column(
            width = 4,
            #numericInput("p.CefGen_nocul_nocha_stepup_nochange", "Probability of stepping up with non-bacterial infection, whose condition has no change under ceftriaxone/gentamicin treatment", value = 0.5, min = 0, max = 1, step = 0.05)
          )
        )
      )),
      
      # Economic Parameters Tab
      tabItem(
        tabName = "economic_params",
        fluidRow(
          bs4Card(
            title = "Cost Parameters",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            fluidRow(
              column(
                width = 6,
                numericInput("cost_henced_microlab", "Additional cost of maintaining microbiology laboratory per sample collected for microbiology testing ($)", value = 34.54, min = 0, step = 10),
                numericInput("cost_CefGen", "Cost of 1-day ceftriaxone/gentamicin ($)", value = 1.76, min = 0, step = 0.1),
                numericInput("cost_Merope", "Cost of 1-day meropenem/vancomycin ($)", value = 35.16, min = 0, step = 0.1)
              ),
              column(
                width = 6,
                numericInput("cost_oralatb", "Cost of oral antibiotics ($)", value = 0.14, min = 0, step = 0.01),
                numericInput("cost_cultpos", "Cost of processing culture positive sample ($)", value = 60, min = 0, step = 1),
                numericInput("cost_cultneg", "Cost of processing culture negative sample ($)", value = 10, min = 0, step = 1)
              )
            )
          )
        ),
        fluidRow(
          bs4Card(
            title = "Treatment Duration Parameters",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            fluidRow(
              column(
                width = 6,
                numericInput("length_targeted", "Length of Targeted Treatment (days)", value = 7, min = 1, max = 14, step = 1)
              )
            )
          )
        )
      ),
      
      # About Tab
      tabItem(
        tabName = "about",
        bs4Card(
          title = "About This Application",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = FALSE,
          h4("Cost-Effectiveness Analysis of Microbiology Laboratory Services"),
          p("This application was developed to conduct show the cost-effectiveness of maintaining active microbiology laboratory services for blood culture testing."),
          p("The model uses a decision tree approach to compare outcomes between two strategies:"),
          tags$ul(
            tags$li("Active microbiology laboratory services"),
            tags$li("No active microbiology laboratory services")
          ),
          p("The model considers various parameters including:"),
          tags$ul(
            tags$li("Probability of bacterial infection"),
            tags$li("Clinical response to treatment"),
            tags$li("Treatment decision pathways"),
            tags$li("Costs of laboratory services and treatments"),
            tags$li("Health outcomes (mortality and disability-adjusted life years)")
          ),
          p("Results are presented in the form of a cost-effectiveness plane, showing the incremental costs against DALYs averted per 1,000 patients.")
        )
      )
    )
  ),
  controlbar = dashboardControlbar(
    skin = "light",
    pinned = FALSE,
    collapsed = TRUE,
    div(
      class = "p-3",
      h4("Analysis Settings"),
      actionButton(
        "reset", 
        "Reset All Parameters", 
        icon = icon("refresh"), 
        status = "danger"
      )
    )
  ),
  footer = dashboardFooter(
    left = "Cost-Effectiveness Analysis of Microbiology Lab Services",
    right = "Health Economics Analysis Tool"
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Calculate derived probabilities for display in valueboxes
  output$prob_no_bact <- renderValueBox({
    valueBox(
      value = 1 - input$p.bact,
      subtitle = "Probability of No Bacterial Infection",
      color = "primary",
      icon = icon("microbes")
    )
  })
  
  output$prob_culpos_nochange <- renderValueBox({
    valueBox(
      subtitle = "Probability of no change in ceftriaxone/gentamicin for patient with culture positive result",
      value = 1 - input$p.CefGen_culpos_stepup,
      color = "primary",
      icon = icon("pills")
    )
  })
  
  output$prob_bact_culneg_nocha <- renderValueBox({
    value <- 1 - input$p.CefGen_bact_culneg_deter - input$p.CefGen_bact_culneg_impro
    color <- if(value < 0) "danger" else "primary"
    
    valueBox(
      value = max(0, value),
      subtitle = "Probability of no clinical change",
      color = color,
      icon = icon("equals")
    )
  })
  
  output$prob_nobact_culneg_nocha <- renderValueBox({
    value <- 1 - input$p.CefGen_nobact_culneg_deter - input$p.CefGen_nobact_culneg_impro
    color <- if(value < 0) "danger" else "primary"
    
    valueBox(
      value = max(0, value),
      subtitle = "Probability of no clinical change",
      color = color,
      icon = icon("equals")
    )
  })
  
  # output$prob_bact_nocul_nocha <- renderValueBox({
  #   value <- 1 - input$p.CefGen_bact_nocul_deter - input$p.CefGen_bact_nocul_impro
  #   color <- if(value < 0) "danger" else "primary"
  #   
  #   valueBox(
  #     value = max(0, value),
  #     subtitle = "Probability of no clinical change",
  #     color = color,
  #     icon = icon("equals")
  #   )
  # })
  
  output$prob_nobact_nocul_nocha <- renderValueBox({
    value <- 1 - input$p.CefGen_nobact_nocul_deter - input$p.CefGen_nobact_nocul_impro
    color <- if(value < 0) "danger" else "primary"
    
    valueBox(
      value = max(0, value),
      subtitle = "Probability of no clinical change",
      color = color,
      icon = icon("equals")
    )
  })
  
  # Reset button functionality
  observeEvent(input$reset, {
    # Reset all input values to their defaults
    updateNumericInput(session, "p.bact", value = 0.5)
    updateNumericInput(session, "p.bact_culture_pos", value = 0.5)
    updateNumericInput(session, "p.nobact_culture_pos", value = 0.5)
    #updateNumericInput(session, "p.CefGen_culpos_stepup", value = 0.5)
    updateNumericInput(session, "p.CefGen_bact_culneg_deter", value = 0.5)
    updateNumericInput(session, "p.CefGen_bact_culneg_impro", value = 0.5)
    updateNumericInput(session, "p.CefGen_nobact_culneg_deter", value = 0.5)
    updateNumericInput(session, "p.CefGen_nobact_culneg_impro", value = 0.5)
    updateNumericInput(session, "p.culneg_deter_stepup_nochange", value = 0.5)
    updateNumericInput(session, "p.culneg_impro_stepdown_nochange", value = 0.5)
    updateNumericInput(session, "p.culneg_nocha_stepup_nochange", value = 0.5)
    #updateNumericInput(session, "p.CefGen_bact_nocul_deter", value = 0.5)
    #updateNumericInput(session, "p.CefGen_bact_nocul_impro", value = 0.5)
    updateNumericInput(session, "p.CefGen_nobact_nocul_deter", value = 0.5)
    updateNumericInput(session, "p.CefGen_nobact_nocul_impro", value = 0.5)
    updateNumericInput(session, "p.CefGen_culneg_deter_stepup_nochange", value = 0.5)
    updateNumericInput(session, "p.CefGen_culneg_impro_stepdown_nochange", value = 0.5)
    updateNumericInput(session, "p.CefGen_culneg_nocha_stepup_nochange", value = 0.5)
    updateNumericInput(session, "p.CefGen_nocul_deter_stepup_nochange", value = 0.012)
    updateNumericInput(session, "p.CefGen_nocul_impro_stepdown_nochange", value = 0.013)
   # updateNumericInput(session, "p.CefGen_nocul_nocha_stepup_nochange", value = 0.5)
    updateNumericInput(session, "cost_henced_microlab", value = 300)
    updateNumericInput(session, "cost_CefGen", value = 1.76)
    updateNumericInput(session, "cost_Merope", value = 35.16)
    updateNumericInput(session, "cost_oralatb", value = 0.14)
    updateNumericInput(session, "cost_cultpos", value = 60)
    updateNumericInput(session, "cost_cultneg", value = 10)
    updateNumericInput(session, "length_targeted", value = 7)
    updateNumericInput(session, "N", value = 1000)
    updateSliderInput(session, "will_pay_perDALY", value = 1500)
    updateSliderInput(session, "SA_prop_bsi", value = 0.05)
    updateSliderInput(session, "PA_prop_bsi", value = 0.05)
    updateSliderInput(session, "SP_prop_bsi", value = 0.05)
    updateSliderInput(session, "NM_prop_bsi", value = 0.05)
    updateSliderInput(session, "AB_prop_bsi", value = 0.05)
    updateSliderInput(session, "EC_prop_bsi", value = 0.05)
    updateSliderInput(session, "kleb_prop_bsi", value = 0.05)
    updateSliderInput(session, "SA_prop_line", value = 0.05)
    updateSliderInput(session, "PA_prop_line", value = 0.05)
    updateSliderInput(session, "SP_prop_line", value = 0.05)
    updateSliderInput(session, "NM_prop_line", value = 0.05)
    updateSliderInput(session, "AB_prop_line", value = 0.05)
    updateSliderInput(session, "EC_prop_line", value = 0.05)
    updateSliderInput(session, "kleb_prop_line", value = 0.05)
    updateSliderInput(session, "SLIDERX", value = 2000)
    updateSliderInput(session, "SLIDERY", value = 8000000)
    
    # Show notification to user
    showNotification("All parameters have been reset to default values", type = "warning")
  })
  
  # Reactive expression to run the model when the "Run Model" button is clicked
  model_results <- eventReactive(input$run, {
    # Show a progress notification
    #showNotification("Running model calculations...", type = "message", duration = NULL, id = "calc_progress")
    
    # Create a list of input parameters
    params <- list(
      p.bact = input$p.bact,
      p.nobact = 1-input$p.bact,
      p.nobact_culture_pos = input$p.nobact_culture_pos,
      p.bact_culture_pos = input$p.bact_culture_pos,
      
       ## Branches of culture positives
      #removed for shiny app
      #p.CefGen_culpos_stepdown = input$p.CefGen_culpos_stepdown
      
      p.CefGen_culpos_nochange = (input$kleb_prop_bsi + input$SA_prop_bsi + input$PA_prop_bsi + input$SP_prop_bsi + 
                                    input$NM_prop_bsi + input$AB_prop_bsi + input$EC_prop_bsi)*(input$kleb_prop_line + input$SA_prop_line + input$PA_prop_line + input$SP_prop_line + 
                                                                                                  input$NM_prop_line + input$AB_prop_line + input$EC_prop_line),
      
      #here we edit the step up fxn to be the sum of the 
      p.CefGen_culpos_stepup = 1- (input$kleb_prop_bsi + input$SA_prop_bsi + input$PA_prop_bsi + input$SP_prop_bsi + 
                                     input$NM_prop_bsi + input$AB_prop_bsi + input$EC_prop_bsi)*(input$kleb_prop_line + input$SA_prop_line + input$PA_prop_line + input$SP_prop_line + 
                                                                                                   input$NM_prop_line + input$AB_prop_line + input$EC_prop_line), 
      #p.CefGen_culpos_stepup = input$p.CefGen_culpos_stepup,
      #p.CefGen_culpos_nochange = 1- input$p.CefGen_culpos_stepup, 
      
      ## Branches of culture negatives 
      p.CefGen_bact_culneg_deter = input$p.CefGen_bact_culneg_deter,
      
      #setting this to the same setup value for now 
      # p.CefGen_bact_culneg_deter = 1- (input$kleb_prop_bsi + input$SA_prop_bsi + input$PA_prop_bsi + input$SP_prop_bsi + 
      #       input$NM_prop_bsi + input$AB_prop_bsi + input$EC_prop_bsi)*(input$kleb_prop_line + input$SA_prop_line + input$PA_prop_line + input$SP_prop_line + 
      #                                                                     input$NM_prop_line + input$AB_prop_line + input$EC_prop_line),
      p.CefGen_bact_culneg_impro = input$p.CefGen_bact_culneg_impro,
      p.CefGen_bact_culneg_nocha = 1 - input$p.CefGen_bact_culneg_deter - input$p.CefGen_bact_culneg_impro,
      
      #setting this to the same setup value for now 
      # p.CefGen_nobact_culneg_deter = 1- (input$kleb_prop_bsi + input$SA_prop_bsi + input$PA_prop_bsi + input$SP_prop_bsi + 
      #                                    input$NM_prop_bsi + input$AB_prop_bsi + input$EC_prop_bsi)*(input$kleb_prop_line + input$SA_prop_line + input$PA_prop_line + input$SP_prop_line + 
      #                                                                                                  input$NM_prop_line + input$AB_prop_line + input$EC_prop_line),
      # 
      p.CefGen_nobact_culneg_deter = input$p.CefGen_nobact_culneg_deter,
      p.CefGen_nobact_culneg_impro = input$p.CefGen_nobact_culneg_impro,
      p.CefGen_nobact_culneg_nocha = 1 - input$p.CefGen_nobact_culneg_deter - input$p.CefGen_nobact_culneg_impro,
      
      p.CefGen_culneg_deter_stepup = input$p.CefGen_culneg_deter_stepup_nochange,
      p.CefGen_culneg_deter_nochange = 1-input$p.CefGen_culneg_deter_stepup_nochange,
      p.CefGen_culneg_impro_stepdown = input$p.CefGen_culneg_impro_stepdown_nochange,
      p.CefGen_culneg_impro_nochange = 1-input$p.CefGen_culneg_impro_stepdown_nochange,
      p.CefGen_culneg_nocha_stepup = input$p.CefGen_culneg_nocha_stepup_nochange,
      p.CefGen_culneg_nocha_nochange = 1-input$p.CefGen_culneg_nocha_stepup_nochange,
      
      p.culneg_deter_stepup = input$p.culneg_deter_stepup_nochange,
      p.culneg_deter_nochange = 1-input$p.culneg_deter_stepup_nochange,
      p.culneg_impro_stepdown = input$p.culneg_impro_stepdown_nochange,
      p.culneg_impro_nochange = 1-input$p.culneg_impro_stepdown_nochange,
      p.culneg_nocha_stepup = input$p.culneg_nocha_stepup_nochange,
      p.culneg_nocha_nochange = 1-input$p.culneg_nocha_stepup_nochange,
      
      ## Branches of no culture performed
      p.CefGen_bact_nocul_deter = 1- (input$kleb_prop_bsi + input$SA_prop_bsi + input$PA_prop_bsi + input$SP_prop_bsi + 
                                        input$NM_prop_bsi + input$AB_prop_bsi + input$EC_prop_bsi)*(input$kleb_prop_line + input$SA_prop_line + input$PA_prop_line + input$SP_prop_line +  input$NM_prop_line + input$AB_prop_line + input$EC_prop_line),
      
      
      p.CefGen_bact_nocul_impro = (input$kleb_prop_bsi + input$SA_prop_bsi + input$PA_prop_bsi + input$SP_prop_bsi + 
                                     input$NM_prop_bsi + input$AB_prop_bsi + input$EC_prop_bsi)*(input$kleb_prop_line + input$SA_prop_line + input$PA_prop_line + input$SP_prop_line +  input$NM_prop_line + input$AB_prop_line + input$EC_prop_line),
      
      #p.CefGen_bact_nocul_deter = input$p.CefGen_bact_nocul_deter,
      #p.CefGen_bact_nocul_impro = input$p.CefGen_bact_nocul_impro,
      #p.CefGen_bact_nocul_nocha = 1-input$p.CefGen_bact_nocul_deter-input$p.CefGen_bact_nocul_impro,
      p.CefGen_nobact_nocul_deter = input$p.CefGen_nobact_nocul_deter,
      p.CefGen_nobact_nocul_impro = input$p.CefGen_nobact_nocul_impro,
      #p.CefGen_nobact_nocul_nocha = 1-input$p.CefGen_nobact_nocul_deter-input$p.CefGen_nobact_nocul_impro,
      
      p.CefGen_nocul_deter_stepup = input$p.CefGen_nocul_deter_stepup_nochange,
      p.CefGen_nocul_deter_nochange = 1-input$p.CefGen_nocul_deter_stepup_nochange,
      
      p.CefGen_nocul_impro_stepdown = input$p.CefGen_nocul_impro_stepdown_nochange,
      p.CefGen_nocul_impro_nochange = 1-input$p.CefGen_nocul_impro_stepdown_nochange,
      
      #p.CefGen_nocul_nocha_stepup = input$p.CefGen_nocul_nocha_stepup_nochange,
      #p.CefGen_nocul_nocha_nochange = 1-input$p.CefGen_nocul_nocha_stepup_nochange,
      
      cost_henced_microlab = input$cost_henced_microlab*2,
      cost_CefGen = input$cost_CefGen,
      cost_Merope = input$cost_Merope,
      cost_oralatb = input$cost_oralatb,
      
      cost_cultpos = input$cost_cultpos,
      cost_cultneg = input$cost_cultneg,
      cost_careunchan = 158.82,
      cost_caredetrio = 203.19,
      cost_careimprov = 14.44,
      
      los_bact_culpos_stepdown = 7,
      los_bact_culpos_stepup = 7*1.5,
      los_bact_culpos_nochange = 7,
      
      ## Culture negative in-hospital mortality based on clinical condition
      los_bact_culneg_stepup_deter = 10.5*2,
      los_bact_culneg_nochange_deter = 10.5*1.5,
      los_bact_culneg_stepdown_impro = 7,
      los_bact_culneg_nochange_impro = 7,
      los_bact_culneg_stepup_uncha = 7,
      los_bact_culneg_nochange_uncha = 7,
      
      ## No microbiology testing performed: in-hospital mortality based on clinical condition
      los_bact_nocul_stepup_deter = 10.5*1.25,
      los_bact_nocul_nochange_deter = 10.5*2,
      los_bact_nocul_stepdown_impro = 10.5*1.5,
      los_bact_nocul_nochange_impro = 10.5,
      los_bact_nocul_stepup_uncha = 13.125,
      los_bact_nocul_nochange_uncha = 10.5*2*1.5,
      
      # Patient with no bacterial infection
      ## Culture positive: in-hospital losality after step-down, step-up, or no change
      los_nobact_culpos_stepdown = 7/2,
      los_nobact_culpos_stepup = 7*1.5/2,
      los_nobact_culpos_nochange = 7/2,
      
      ## Culture negative: in-hospital losality based on clinical condition (for culture negatives)
      los_nobact_culneg_stepup_deter = 10.5*2/2,
      los_nobact_culneg_nochange_deter = 10.5*1.5/2,
      los_nobact_culneg_stepdown_impro = 7/2,
      los_nobact_culneg_nochange_impro = 7/2,
      los_nobact_culneg_stepup_uncha = 7/2,
      los_nobact_culneg_nochange_uncha = 7/2,
      
      ## No microbiology testing performed: in-hospital losality based on clinical condition
      los_nobact_nocul_stepup_deter = 10.5*2/2,
      los_nobact_nocul_nochange_deter = 10.5*1.5/2,
      los_nobact_nocul_stepdown_impro = 7/2,
      los_nobact_nocul_nochange_impro = 7/2,
      los_nobact_nocul_stepup_uncha = 7/2,
      los_nobact_nocul_nochange_uncha = 7/2,
      
      # Mortality of each path ----
      # Patient with bacterial infection
      ## Culture positive: in-hospital mortality after step-down, step-up, or no change
      mort_bact_culpos_stepdown = 0.2,
      mort_bact_culpos_stepup = 0.2*1.25,
      mort_bact_culpos_nochange = 0.2*1.25,
      ## Culture negative in-hospital mortality based on clinical condition
      mort_bact_culneg_stepup_deter = 0.2*1.25,
      mort_bact_culneg_nochange_deter = 0.2*1.25*1.5,
      mort_bact_culneg_stepdown_impro = 0.2, # mort_bact_culneg_stepup_deter/3,
      mort_bact_culneg_nochange_impro = 0.2, # mort_bact_culneg_stepup_deter/3,
      mort_bact_culneg_stepup_uncha = 0.2*1.25, # mort_bact_culneg_stepup_deter/3,
      mort_bact_culneg_nochange_uncha = 0.2*1.25, #mort_bact_culneg_stepup_deter/2,
      ## No microbiology testing performed: in-hospital mortality based on clinical condition
      mort_bact_nocul_stepup_deter = 0.2*1.25*1.25,
      mort_bact_nocul_nochange_deter = 0.2*1.25*2,
      mort_bact_nocul_stepdown_impro = 0.2*1.25*1.5,
      mort_bact_nocul_nochange_impro = 0.2*1.25*1.25,
      mort_bact_nocul_stepup_uncha = 0.2*1.25*1.25,
      mort_bact_nocul_nochange_uncha = 0.2*1.25*1.5,
      
      # Patient with no bacterial infection
      ## Culture positive: in-hospital mortality after step-down, step-up, or no change
      mort_nobact_culpos_stepdown = 0.2/2,
      mort_nobact_culpos_stepup = 0.2*1.25/2,
      mort_nobact_culpos_nochange = 0.2*1.25/2,
      ## Culture negative: in-hospital mortality based on clinical condition (for culture negatives)
      mort_nobact_culneg_stepup_deter = 0.2*1.25/2,
      mort_nobact_culneg_nochange_deter = 0.2*1.25*1.5/2,
      mort_nobact_culneg_stepdown_impro = 0.2/2,
      mort_nobact_culneg_nochange_impro = 0.2/2,
      mort_nobact_culneg_stepup_uncha =  0.2*1.25/2,
      mort_nobact_culneg_nochange_uncha = 0.2*1.25/2,
      ## No microbiology testing performed: in-hospital mortality based on clinical condition
      mort_nobact_nocul_stepup_deter = 0.2*1.25/2,
      mort_nobact_nocul_nochange_deter = 0.2*1.25*1.5/2,
      mort_nobact_nocul_stepdown_impro =  0.2/2,
      mort_nobact_nocul_nochange_impro = 0.2/2,
      mort_nobact_nocul_stepup_uncha = 0.2*1.25/2,
      mort_nobact_nocul_nochange_uncha = 0.2*1.25/2,
      
      # DALY of each path ----
      # DALY = YLL + YLD
      # YLL is calculated as the number of deaths (n) x the standard life expectancy at age of death (L1). 
      # This measures the reduction in life expectancy.
      # Source: https://www.who.int/data/gho/data/indicators/indicator-details/GHO/gho-ghe-life-tables-by-country
      # Patient with bacterial infection
      ## Culture positive: in-hospital mortality after step-down, step-up, or no change
      year_lost = 11,
      yll_bact_culpos_stepdown = 0.2*11,       # expectation of life at age 55-59 years
      yll_bact_culpos_stepup = 0.2*1.25*11,         
      yll_bact_culpos_nochange = 0.2*1.25*11,   
      ## Culture negative in-hospital mortality based on clinical condition
      yll_bact_culneg_stepup_deter = 0.2*1.25*11,                
      yll_bact_culneg_nochange_deter = 0.2*1.25*1.5*11, 
      yll_bact_culneg_stepdown_impro = 0.2*11,     
      yll_bact_culneg_nochange_impro = 0.2*11, 
      yll_bact_culneg_stepup_uncha = 0.2*1.25*11,   
      yll_bact_culneg_nochange_uncha = 0.2*1.25*11,   
      ## No microbiology testing performed: in-hospital mortality based on clinical condition
      yll_bact_nocul_stepup_deter = 0.2*1.25*1.25*11,                
      yll_bact_nocul_nochange_deter = 0.2*1.25*2*11, 
      yll_bact_nocul_stepdown_impro = 0.2*1.25*1.5*11,     
      yll_bact_nocul_nochange_impro = 0.2*1.25*1.25*11, 
      yll_bact_nocul_stepup_uncha = 0.2*1.25*1.25*11,   
      yll_bact_nocul_nochange_uncha = 0.2*1.25*1.5*11,   
      
      # Patient with no bacterial infection
      ## Culture positive: in-hospital mortality after step-down, step-up, or no change
      yll_nobact_culpos_stepdown = 0.2/2*11,       # expectation of life at age 55-59 years
      yll_nobact_culpos_stepup = 0.2*1.25/2*11,         
      yll_nobact_culpos_nochange = 0.2*1.25/2*11,   
      ## Culture negative: in-hospital mortality based on clinical condition (for culture negatives)
      yll_nobact_culneg_stepup_deter = 0.2*1.25/2*11,                
      yll_nobact_culneg_nochange_deter = 0.2*1.25*1.5/2*11, 
      yll_nobact_culneg_stepdown_impro = 0.2/2*11,     
      yll_nobact_culneg_nochange_impro = 0.2/2*11, 
      yll_nobact_culneg_stepup_uncha = 0.2*1.25/2*11,   
      yll_nobact_culneg_nochange_uncha = 0.2*1.25/2*11,   
      ## No microbiology testing performed: in-hospital mortality based on clinical condition
      yll_nobact_nocul_stepup_deter = 0.2*1.25/2*11,                
      yll_nobact_nocul_nochange_deter = 0.2*1.25*1.5/2*11, 
      yll_nobact_nocul_stepdown_impro = 0.2/2*11,     
      yll_nobact_nocul_nochange_impro = 0.2/2*11, 
      yll_nobact_nocul_stepup_uncha =  0.2*1.25/2*11,   
      yll_nobact_nocul_nochange_uncha = 0.2*1.25/2*11,
      
      # Patient with bacterial infection
      ## Culture positive: in-hospital length of hospitalisation after step-down, step-up, or no change
      daly_weight = 0.33,
      yld_bact_culpos_stepdown = 7*0.33,       # expectation of life at age 55-59 years
      yld_bact_culpos_stepup = 7*1.5*0.33,         
      yld_bact_culpos_nochange = 7*0.33,   
      ## Culture negative in-hospital length of hospitalisation based on clinical condition
      yld_bact_culneg_stepup_deter = 10.5*2*0.33,                
      yld_bact_culneg_nochange_deter = 10.5*1.5*0.33, 
      yld_bact_culneg_stepdown_impro = 7*0.33,     
      yld_bact_culneg_nochange_impro = 7*0.33, 
      yld_bact_culneg_stepup_uncha = 7*0.33,   
      yld_bact_culneg_nochange_uncha = 7*0.33,   
      ## No microbiology testing performed: in-hospital length of hospitalisation based on clinical condition
      yld_bact_nocul_stepup_deter = 10.5*1.25*0.33,                
      yld_bact_nocul_nochange_deter = 10.5*2*0.33, 
      yld_bact_nocul_stepdown_impro = 10.5*1.5*0.33,     
      yld_bact_nocul_nochange_impro = 10.5*0.33, 
      yld_bact_nocul_stepup_uncha = 13.125*0.33,   
      yld_bact_nocul_nochange_uncha = 10.5*2*1.5*0.33,   
      
      # Patient with no bacterial infection
      ## Culture positive: in-hospital length of hospitalisation after step-down, step-up, or no change
      yld_nobact_culpos_stepdown = 7/2*0.33,       # expectation of life at age 55-59 years
      yld_nobact_culpos_stepup = 7*1.5/2*0.33,         
      yld_nobact_culpos_nochange = 7/2*0.33,   
      ## Culture negative: in-hospital length of hospitalisation based on clinical condition (for culture negatives)
      yld_nobact_culneg_stepup_deter = 10.5*2/2*0.33,                
      yld_nobact_culneg_nochange_deter = 10.5*1.5/2*0.33, 
      yld_nobact_culneg_stepdown_impro = 7/2*0.33,     
      yld_nobact_culneg_nochange_impro = 7/2*0.33, 
      yld_nobact_culneg_stepup_uncha = 7/2*0.33,   
      yld_nobact_culneg_nochange_uncha = 7/2*0.33,   
      ## No microbiology testing performed: in-hospital length of hospitalisation based on clinical condition
      yld_nobact_nocul_stepup_deter = 10.5*2/2*0.33,                
      yld_nobact_nocul_nochange_deter = 10.5*1.5/2*0.33, 
      yld_nobact_nocul_stepdown_impro = 7/2*0.33,     
      yld_nobact_nocul_nochange_impro = 7/2*0.33, 
      yld_nobact_nocul_stepup_uncha = 7/2*0.33,   
      yld_nobact_nocul_nochange_uncha = 7/2*0.33,
      
      length_targeted = input$length_targeted,
      length_empiric = 3,
      #num_bottle = input$num_bottle,
      N = input$N
      # Add other parameters as needed
      
    )
    
    # Call the dec_tree function with the input parameters
    results <- dec_tree(params)
    
    # Return the results
    results
    
    # Results has the ic/ie -> we want the ic as the y axis, ie x axis to make the ice( slope is icer)
    
  })
  
  # Render the model results
  output$results <- renderPrint({
    model_results()
  })
  
  output$plotCEP <-   renderPlot({
    
    # this should give the results produced by the decision tree model 
    data <- as.data.frame(t(model_results()))
    colnames(data) <- c("diff_im", "diff_cost", "diff_daly", "icer",
                        "mort.noactivelab_CefGen", "mort.activelab_CefGen",
                        "cost.noactivelab_CefGen", "cost.activelab_CefGen",
                        "daly.noactivelab_CefGen", "daly.activelab_CefGen")
    ## Plot CEA
    ggplot(data) +
      geom_point(aes(x=diff_daly/1000, y=diff_cost), size = 20, colour = "red",alpha = 0.8, shape = 21,stroke = 8) + 
      geom_label_repel(
        label="Example Intervention",
        data = data,
        x = data$diff_daly/1000, y=data$diff_cost, #we divide the diff_daly by 1000 since we want it per 1000 people
        color = "black", fill= "lightgrey", box.padding = 4,
        nudge_x = 0, nudge_y = 20
      ) + labs(
        y="Incremental Costs ($)",
        x="DALY Averted /1,000 patients",
        title = paste0("Cost-Effectiveness Plane - Intervention Incremental Cost: ", data$diff_cost,"($),  DALY Averted/1000 Patients: ", data$diff_daly/1000))+
      # xlim(-2000, 2000) +
      # ylim(-8000000, 8000000) +
      xlim(-(input$SLIDERX), input$SLIDERX) +
      ylim(-(input$SLIDERY), input$SLIDERY) +
      theme_minimal() +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=12)) +
      #scale_y_continuous(labels = comma, limits = c(round(min(dta_ceplane_CefGen2$ic), 0)*1.25, 100000)) +
      geom_hline(yintercept=0, colour = "darkgrey") +
      geom_vline(xintercept=0, colour = "darkgrey") +
      geom_textabline(slope = input$will_pay_perDALY, intercept = -100, label = "Cost-Effectiveness Threshold",
                      color = "blue", hjust = 0.8, vjust = -0.3,linetype= 2, linewidth = 1.5, size = 5) 
    
  })
}

# Run the application
shinyApp(ui = ui, server = server)