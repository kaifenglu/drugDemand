library(shiny)
library(shinyMatrix)
library(shinyFeedback)
library(shinyjs, warn.conflicts = FALSE)
library(shinybusy)
library(readxl)
library(writexl)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(purrr)
library(prompter)
library(ggplot2)
library(plotly, warn.conflicts = FALSE)
library(eventPred)
library(drugDemand)


# conditional panels for treatment allocation
f_treatment_allocation <- function(i) {
  conditionalPanel(
    condition = paste0("input.k == ", i),

    shinyMatrix::matrixInput(
      paste0("treatment_allocation_", i),
      label = tags$span(
        "Treatment allocation",
        tags$span(icon(name = "question-circle")) %>%
          add_prompt(message = "in a randomization block",
                     position = "right")),

      value = matrix(rep(1,i), ncol = 1,
                     dimnames = list(paste("Treatment", 1:i), "Size")),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE, editableNames=TRUE),
      cols = list(names=TRUE, extend=FALSE))
  )
}


f_exponential_survival <- function(i) {
  conditionalPanel(
    condition = paste("input.event_prior == 'Exponential' && input.k ==", i),

    shinyMatrix::matrixInput(
      paste0("exponential_survival_", i),
      label = "Hazard rate for each treatment",
      value = matrix(rep(0.0030, i), nrow = 1,
                     dimnames = list(NULL, paste("Treatment", 1:i))),
      inputClass = "numeric",
      rows = list(names=FALSE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}


f_weibull_survival <- function(i) {
  conditionalPanel(
    condition = paste("input.event_prior == 'Weibull' && input.k ==", i),

    shinyMatrix::matrixInput(
      paste0("weibull_survival_", i),
      label = "Weibull parameters",
      value = matrix(rep(c(1.42, 392), i), nrow = 2, byrow = FALSE,
                     dimnames = list(c("Shape", "Scale"),
                                     paste("Treatment", 1:i))),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}


f_llogis_survival <- function(i) {
  conditionalPanel(
    condition = paste("input.event_prior == 'Log-logistic' && input.k ==", i),

    shinyMatrix::matrixInput(
      paste0("llogis_survival_", i),
      label = "Log-logistic parameters",
      value = matrix(rep(c(5.4, 1), i), nrow = 2, byrow = FALSE,
                     dimnames = list(c("Location on log scale",
                                       "Scale on log scale"),
                                     paste("Treatment", 1:i))),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}


f_lnorm_survival <- function(i) {
  conditionalPanel(
    condition = paste("input.event_prior == 'Log-normal' && input.k ==", i),

    shinyMatrix::matrixInput(
      paste0("lnorm_survival_", i),
      label = "Log-normal parameters",
      value = matrix(rep(c(5.4, 1), i), nrow = 2, byrow = FALSE,
                     dimnames = list(c("Mean on log scale",
                                       "SD on log scale"),
                                     paste("Treatment", 1:i))),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}


f_piecewise_exponential_survival <- function(i) {
  conditionalPanel(
    condition = paste(
      "input.event_prior == 'Piecewise exponential' && input.k ==", i),

    shinyMatrix::matrixInput(
      paste0("piecewise_exponential_survival_", i),
      label = "Hazard rate by time interval for each treatment",
      value = matrix(c(0, rep(0.0030, i)), nrow = 1,
                     dimnames = list(
                       "Interval 1",
                       c("Starting time", paste("Treatment", 1:i)))),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    ),

    actionButton(paste0("add_piecewise_exponential_survival_", i),
                 label=NULL, icon=icon("plus")),
    actionButton(paste0("del_piecewise_exponential_survival_", i),
                 label=NULL, icon=icon("minus"))
  )
}


f_drug_description <- function(j) {
  conditionalPanel(
    condition = paste("input.l ==", j),

    shinyMatrix::matrixInput(
      paste0("drug_description_", j),
      label = "Drug names and dose units",
      value = matrix(c(paste("Drug", seq_len(j)), rep("kit", j)),
                     nrow = j, ncol = 2,
                     dimnames = list(NULL, c("Drug Name", "Dose Unit"))),
      inputClass = "character",
      rows = list(names=FALSE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE))
  )
}


f_treatment_by_drug <- function(i, j) {
  conditionalPanel(
    condition = paste("input.k ==", i, "&&", "input.l ==", j),

    shinyMatrix::matrixInput(
      paste0("treatment_by_drug_", i, "_", j),
      label = "Drugs contained in each treatment",
      value = matrix(1, nrow = i, ncol = j,
                     dimnames = list(paste("Treatment", 1:i),
                                     paste("Drug", 1:j))),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}


f_dosing_schedule <- function(j) {
  conditionalPanel(
    condition = paste("input.l ==", j),

    shinyMatrix::matrixInput(
      paste0("dosing_schedule_", j),
      label = "Dosing schedule",
      value = matrix(rep(c(21, 2, 10000), each = j),
                     nrow = j, ncol = 3,
                     dimnames = list(paste("Drug", 1:j),
                                     c("Days per Cycle", "Dose per Cycle",
                                       "Number of Cycles"))),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE))
  )
}


# cumulative dose for duration x and a drug with dosing schedule (w, d, N)
# here x can be a vector
f_cum_dose <- function(x, w, d, N) {
  m = length(w)
  u = c(0, cumsum(w*N))
  i = pmin(findInterval(x, u), m)
  z = pmin(x, u[m+1]) - u[i]
  n = pmin(floor(z/w[i]) + 1, N[i])
  v = c(0, cumsum(d*N))
  v[i] + n*d[i]
}


observed_enrollment_event_data_panel <- tabPanel(
  title = "Enrollment and Event",
  value = "enrollment_event_data_panel",

  htmlOutput("dates"),
  verbatimTextOutput("statistics"),

  plotlyOutput("cum_accrual_plot"),

  conditionalPanel(
    condition = "input.stage != 'Real-time after enrollment completion'",
    plotlyOutput("daily_accrual_plot")
  ),

  conditionalPanel(
    condition = "input.to_predict == 'Enrollment and event' ||
    input.stage == 'Real-time after enrollment completion'",

    plotlyOutput("event_km_plot")
  ),

  conditionalPanel(
    condition = "input.stage != 'Design stage'",
    dataTableOutput("input_df"))
)


observed_dosing_data_panel <- tabPanel(
  title = "Drug Dispensing",
  value = "dosing_data_panel",

  htmlOutput("dates_copy"),
  verbatimTextOutput("cum_dose_t0"),

  plotlyOutput("cum_dispense_plot"),
  plotlyOutput("bar_t0_plot"),
  plotlyOutput("bar_ti_plot"),
  plotlyOutput("bar_di_plot"),

  conditionalPanel(
    condition = "input.stage != 'Design stage'",
    dataTableOutput("input_visitview"))
)


observedPanel <- tabPanel(
  title = "Observed Data",
  value = "observed_data_panel",

  tabsetPanel(
    id = "observed_data_tabset",

    observed_enrollment_event_data_panel,
    observed_dosing_data_panel
  )
)


enrollmentPanel <- tabPanel(
  title = "Enrollment Model",
  value = "enroll_model_panel",

  conditionalPanel(
    condition = "input.stage == 'Design stage'",

    fluidRow(
      column(6, radioButtons(
        "enroll_prior",
        label = "Which enrollment model to use?",
        choices = c("Poisson",
                    "Time-decay",
                    "Piecewise Poisson"),
        selected = "Piecewise Poisson",
        inline = FALSE)
      ),

      column(6,
             conditionalPanel(
               condition = "input.enroll_prior == 'Poisson'",

               numericInput(
                 "poisson_rate",
                 label = "Daily enrollment rate",
                 value = 1,
                 min = 0, max = 100, step = 1)
             ),

             conditionalPanel(
               condition = "input.enroll_prior == 'Time-decay'",

               fluidRow(
                 column(6, numericInput(
                   "mu",
                   label = "Base rate, mu",
                   value = 1.5,
                   min = 0, max = 100, step = 1)
                 ),

                 column(6, numericInput(
                   "delta",
                   label = "Decay rate, delta",
                   value = 2,
                   min = 0, max = 100, step = 1)
                 )
               )
             ),

             conditionalPanel(
               condition = "input.enroll_prior == 'Piecewise Poisson'",

               shinyMatrix::matrixInput(
                 "piecewise_poisson_rate",
                 label = "Daily enrollment rate by time interval",
                 value = matrix(c(0,1), ncol = 2,
                                dimnames = list("Interval 1",
                                                c("Starting time",
                                                  "Enrollment rate"))),
                 inputClass = "numeric",
                 rows = list(names=TRUE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)),

               actionButton("add_piecewise_poisson_rate",
                            label=NULL, icon=icon("plus")),
               actionButton("del_piecewise_poisson_rate",
                            label=NULL, icon=icon("minus"))
             )
      )
    )
  ),

  conditionalPanel(
    condition = "input.stage != 'Design stage'",

    fluidRow(
      column(6, radioButtons(
        "enroll_model",
        label = "Which enrollment model to use?",
        choices = c("Poisson",
                    "Time-decay",
                    "B-spline",
                    "Piecewise Poisson"),
        selected = "B-spline",
        inline = FALSE)
      ),

      column(6,
             conditionalPanel(
               condition = "input.enroll_model == 'B-spline'",

               numericInput(
                 "nknots",
                 label = "How many inner knots to use?",
                 value = 0,
                 min = 0, max = 10, step = 1),

               numericInput(
                 "lags",
                 label = paste("How many days before the last enrollment",
                               "date to average",
                               "the enrollment rate over for prediction?"),
                 value = 30,
                 min = 0, max = 365, step = 1)
             ),

             conditionalPanel(
               condition = "input.enroll_model == 'Piecewise Poisson'",

               shinyMatrix::matrixInput(
                 "accrualTime",
                 label = "What is the starting time of each time interval?",
                 value = matrix(0, ncol = 1,
                                dimnames = list("Interval 1",
                                                "Starting time")),
                 inputClass = "numeric",
                 rows = list(names=TRUE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)),

               actionButton("add_accrualTime",
                            label=NULL, icon=icon("plus")),
               actionButton("del_accrualTime",
                            label=NULL, icon=icon("minus"))
             )
      )
    ),

    plotlyOutput("enroll_fit")
  )
)


eventPanel <- tabPanel(
  title = "Event Model",
  value = "event_model_panel",

  conditionalPanel(
    condition = "input.stage == 'Design stage'",

    fluidRow(
      column(4, radioButtons(
        "event_prior",
        label = "Which time-to-event model to use?",
        choices = c("Exponential",
                    "Weibull",
                    "Log-logistic",
                    "Log-normal",
                    "Piecewise exponential"),
        selected = "Piecewise exponential",
        inline = FALSE)
      ),

      column(8,
             lapply(1:6, f_exponential_survival),
             lapply(1:6, f_weibull_survival),
             lapply(1:6, f_llogis_survival),
             lapply(1:6, f_lnorm_survival),
             lapply(1:6, f_piecewise_exponential_survival)
      )
    )
  ),


  conditionalPanel(
    condition = "input.stage != 'Design stage'",

    fluidRow(
      column(6, radioButtons(
        "event_model",
        label = "Which time-to-event model to use?",
        choices = c("Exponential",
                    "Weibull",
                    "Log-logistic",
                    "Log-normal",
                    "Piecewise exponential",
                    "Model averaging",
                    "Spline"),
        selected = "Model averaging",
        inline = FALSE)
      ),

      column(6,
             conditionalPanel(
               condition = "input.event_model == 'Piecewise exponential'",

               shinyMatrix::matrixInput(
                 "piecewiseSurvivalTime",
                 label = "What is the starting time of each time interval?",
                 value = matrix(0, ncol = 1,
                                dimnames = list("Interval 1",
                                                "Starting time")),
                 inputClass = "numeric",
                 rows = list(names=TRUE, extend=FALSE),
                 cols = list(names=TRUE, extend=FALSE)),

               actionButton("add_piecewiseSurvivalTime",
                            label=NULL, icon=icon("plus")),
               actionButton("del_piecewiseSurvivalTime",
                            label=NULL, icon=icon("minus"))
             ),

             conditionalPanel(
               condition = "input.event_model == 'Spline'",

               numericInput(
                 "spline_k",
                 label = "How many inner knots to use?",
                 value = 1,
                 min = 0, max = 10, step = 1),

               radioButtons(
                 "spline_scale",
                 label = "Which scale to model as a spline function?",
                 choices = c("hazard", "odds", "normal"),
                 selected = "hazard",
                 inline = TRUE)
             )
      )
    ),

    uiOutput("event_fit_ic"),
    uiOutput("event_fit")
  )
)


dosingSchedulePanel <- tabPanel(
  title = "Dosing schedule",
  value = "dosing_schedule_panel",

  selectInput(
    "l", label = "Total number of drugs",
    choices = seq_len(12), selected = 2),

  lapply(1:12, f_drug_description),

  lapply(c(outer(1:6, 1:12, function(i,j) 1000*i+j)),
         function(k) f_treatment_by_drug(k %/% 1000, k %% 1000)),

  lapply(1:12, f_dosing_schedule)
)


k0Panel <- tabPanel(
  title = "k0 model",
  value = "k0_model_panel",

  radioButtons(
    "model_k0",
    label = paste("Which model to use for the number of skipped visits",
                  "between randomization and the first drug dispensing",
                  "visit?"),
    choices = c("Poisson",
                "Zero-inflated Poisson",
                "Negative binomial"),
    selected = "Negative binomial",
    inline = FALSE),

  uiOutput("k0_fit")
)


t0Panel <- tabPanel(
  title = "T0 model",
  value = "t0_model_panel",

  radioButtons(
    "model_t0",
    label = paste("Which model to use for the time between randomization",
                  "and the first drug dispensing visit when there is",
                  "no visit skipping?"),
    choices = c("Exponential",
                "Weibull",
                "Log-logistic",
                "Log-normal"),
    selected = "Log-logistic",
    inline = FALSE),

  uiOutput("t0_fit")
)


t1Panel <- tabPanel(
  title = "T1 model",
  value = "t1_model_panel",

  radioButtons(
    "model_t1",
    label = paste("Which model to use for the time between randomization",
                  "and the first drug dispensing visit when there is",
                  "visit skipping?"),
    choices = c("Least squares",
                "Least absolute deviations"),
    selected = "Least squares",
    inline = FALSE),

  uiOutput("t1_fit")
)


kiPanel <- tabPanel(
  title = "ki model",
  value = "ki_model_panel",

  radioButtons(
    "model_ki",
    label = paste("Which model to use for the number of skipped visits",
                  "between two consecutive drug dispensing visits?"),
    choices = c("Poisson",
                "Zero-inflated Poisson",
                "Negative binomial"),
    selected = "Negative binomial",
    inline = FALSE),

  uiOutput("ki_fit")
)


tiPanel <- tabPanel(
  title = "Ti model",
  value = "ti_model_panel",

  radioButtons(
    "model_ti",
    label = paste("Which model to use for the time between two",
                  "consecutive drug dispensing visits?"),
    choices = c("Least squares",
                "Least absolute deviations"),
    selected = "Least absolute deviations",
    inline = FALSE),

  uiOutput("ti_fit")
)


diPanel <- tabPanel(
  title = "di model",
  value = "di_model_panel",

  radioButtons(
    "model_di",
    label = paste("Which model to use for the doses dispensed at drug",
                  "dispensing visits?"),
    choices = c("Linear model",
                "Linear mixed-effects model"),
    selected = "Linear mixed-effects model",
    inline = FALSE),

  uiOutput("di_fit_ic"),
  uiOutput("di_fit")
)


dosingPanel <- tabPanel(
  title = "Dosing Model",
  value = "dosing_model_panel",

  tabsetPanel(
    id = "dosing_tabset",

    dosingSchedulePanel,
    k0Panel,
    t0Panel,
    t1Panel,
    kiPanel,
    tiPanel,
    diPanel
  )
)


eventPredictPanel <- tabPanel(
  title = "Event Prediction",
  value = "event_prediction_panel",

  uiOutput("pred_date"),
  uiOutput("pred_plot"),

  downloadButton("downloadEventSummaryData", "Download summary data"),
  downloadButton("downloadEventSubjectData", "Download subject data")
)


dosingPredictPanel <- tabPanel(
  title = "Dosing Prediction",
  value = "dosing_prediction_panel",

  uiOutput("dosing_plot"),
  downloadButton("downloadDosingSummaryData", "Download summary data"),
  downloadButton("downloadDosingSubjectData",
                 label = tags$span(
                   "Download subject data",
                   tags$span(icon(name = "question-circle")) %>%
                     add_prompt(message = "for the first simulation only",
                                position = "right"))
                 )
)


# reduced style fileInput
fileInputNoExtra<-function(inputId, label, multiple = FALSE, accept = NULL,
                           width = NULL, buttonLabel = "Browse...",
                           placeholder = "No file selected"){

  restoredValue <- restoreInput(id = inputId, default = NULL)
  if (!is.null(restoredValue) && !is.data.frame(restoredValue)) {
    warning("Restored value for ", inputId, " has incorrect format.")
    restoredValue <- NULL
  }
  if (!is.null(restoredValue)) {
    restoredValue <- toJSON(restoredValue, strict_atomic = FALSE)
  }
  inputTag <- tags$input(id = inputId, name = inputId, type = "file",
                         style = "display: none;",
                         `data-restore` = restoredValue)
  if (multiple)
    inputTag$attribs$multiple <- "multiple"
  if (length(accept) > 0)
    inputTag$attribs$accept <- paste(accept, collapse = ",")

  tags$label(
    class = "input-group-btn",
    type="button",
    style=if (!is.null(width))
      paste0("width: ", validateCssUnit(width), ";",
             "padding-right: 5px; padding-bottom: 0px;
             display:inline-block;"),

    span(class = "btn btn-default btn-file",type="button",
         buttonLabel, inputTag,
         style=if (!is.null(width))
           paste0("width: ", validateCssUnit(width), ";",
                  "border-radius: 4px; padding-bottom:5px;"))
  )
}


# user interface ----------------
ui <- fluidPage(
  shinyFeedback::useShinyFeedback(),
  shinyjs::useShinyjs(),
  prompter::use_prompt(),
  shinybusy::add_busy_spinner(),

  titlePanel(tagList(
    span(HTML(paste(tags$span(style="font-size:14pt",
                              "Enrollment and Treatment Prediction"))),
         span(actionButton(
           "predict", "Predict",
           style="color: #fff; background-color: #337ab7;
           border-color: #2e6da4"),

           downloadButton("saveInputs", "Save inputs"),
           fileInputNoExtra("loadInputs", label=NULL, accept=".rds",
                            buttonLabel=list(icon("upload"), "Load inputs"),
                            width="116px"),
           tags$a(tags$span(icon(name = "question-circle")), target="_blank",
                  href="manual.pdf"),
           style="position:absolute;right:0.5em;",
           tags$style(type='text/css', "#saveInputs{margin-top: -5px;}")
         ))),
    windowTitle = "Enrollment and Treatment Prediction"),


  sidebarLayout(
    sidebarPanel(

      fluidRow(
        column(7,
               radioButtons(
                 "stage",
                 label = "Stage of the study",
                 choices = c("Design stage",
                             "Real-time before enrollment completion",
                             "Real-time after enrollment completion"),
                 selected = "Real-time after enrollment completion",
                 inline = FALSE)),

        column(5,
               conditionalPanel(
                 condition = "input.stage == 'Design stage' ||
                 input.stage == 'Real-time before enrollment completion'",

                 radioButtons(
                   "to_predict",
                   label = tags$span(
                     "What to predict?",
                     tags$span(icon(name = "question-circle")) %>%
                       add_prompt(
                         message = "Event refers to treatment discontinuation",
                         position = "right")),

                   choices = c("Enrollment only",
                               "Enrollment and event"),
                   selected = "Enrollment and event",
                   inline = FALSE)
               ),


               conditionalPanel(
                 condition =
                   "input.stage == 'Real-time after enrollment completion'",

                 radioButtons(
                   "to_predict2",
                   label = "What to predict?",
                   choices = c("Event only"),
                   selected = "Event only",
                   inline = FALSE)
               )
        )
      ),


      fluidRow(
        column(7,
               conditionalPanel(
                 condition = "input.stage == 'Design stage' ||
                 input.stage == 'Real-time before enrollment completion'",

                 numericInput(
                   "target_n",
                   label = "Target enrollment",
                   value = 300,
                   min = 1, max = 20000, step = 1)
               )
        ),

        column(5,
               conditionalPanel(
                 condition = "input.to_predict == 'Enrollment and event' ||
                 input.stage == 'Real-time after enrollment completion'",

                 numericInput(
                   "target_d",
                   label = "Target events",
                   value = 200,
                   min = 1, max = 10000, step = 1)
               )
        )
      ),


      conditionalPanel(
        condition = "input.stage != 'Design stage'",

        fileInput(
          "file1",
          label = "Upload enrollment and event data",
          accept = ".xlsx"
        )
      ),


      fluidRow(
        column(7, radioButtons(
          "pilevel",
          label = "Prediction interval",

          choices = c("95%" = "0.95", "90%" = "0.90", "80%" = "0.80"),
          selected = "0.95",
          inline = TRUE)
        ),

        column(5, numericInput(
          "nyears",
          label = "Years after cutoff",
          value = 1,
          min = 1, max = 10, step = 1)
        )
      ),


      conditionalPanel(
        condition = "input.to_predict == 'Enrollment and event' ||
        input.stage == 'Real-time after enrollment completion'",

        checkboxGroupInput(
          "to_show",
          label = "What to show on event prediction plot?",
          choices = c("Enrollment", "Event", "Ongoing"),
          selected = c("Enrollment", "Event"),
          inline = TRUE
        )
      ),


      fluidRow(
        column(7, checkboxInput(
          "by_treatment", label = "By treatment?", value = TRUE),

          conditionalPanel(
            condition = "input.by_treatment &&
            (input.to_predict == 'Enrollment and event' ||
            input.stage == 'Real-time after enrollment completion')",

            checkboxInput(
              "predict_dosing", label = "Predict dosing?", value = TRUE),
          )
        ),


        column(5, conditionalPanel(
          condition = "input.stage == 'Design stage' || input.by_treatment",

          selectInput(
            "k", label = "Treatments", choices = seq_len(6), selected = 2))
        )
      ),


      conditionalPanel(
        condition = "input.stage != 'Design stage' &
        !(input.stage == 'Real-time before enrollment completion' &
        input.to_predict == 'Enrollment only') &
        input.predict_dosing",

        fileInput(
          "file2",
          label = "Upload drug dispensing data",
          accept = ".xlsx"
        ),

        fluidRow(
          column(5, checkboxInput(
            "pred_pp_only",
            label = "Protocol based prediction only?",
            value = FALSE)),

          column(7,
                 conditionalPanel(
                   condition = "!input.pred_pp_only",
                   checkboxGroupInput(
                     "to_show_dosing",
                     label = "What to show on dosing prediction plot?",
                     choices = c("Model based", "Protocol based"),
                     selected = c("Model based"),
                     inline = TRUE
                   )
                 )
          )
        )
      ),


      conditionalPanel(
        condition = "input.stage == 'Design stage' ||
        (input.by_treatment &&
        input.stage != 'Real-time after enrollment completion')",

        lapply(2:6, f_treatment_allocation)
      ),


      fluidRow(
        column(7, numericInput(
          "nreps",
          label = "Simulation runs",
          value = 200,
          min = 100, max = 10000, step = 1)
        ),

        column(5, numericInput(
          "seed",
          label = "Seed",
          value = 2000,
          min = 0, max = 100000, step = 1
        ))
      )
    ),


    mainPanel(
      tabsetPanel(
        id = "results",
        observedPanel,
        enrollmentPanel,
        eventPanel,
        dosingPanel,
        eventPredictPanel,
        dosingPredictPanel
      )
    )
  )
)


# server function -------------
server <- function(input, output, session) {
  # session$onSessionEnded(function() {
  #   stopApp()
  # })


  # whether to show or hide the observed data panel
  observeEvent(input$stage, {
    if (input$stage != 'Design stage') {
      showTab(inputId = "results", target = "observed_data_panel")
    } else {
      hideTab(inputId = "results", target = "observed_data_panel")
    }
  })


  # whether to allow the user to specify the number of treatments
  observeEvent(input$stage, {
    shinyjs::toggleState("k", input$stage == "Design stage")
  })


  # what to predict at different stages
  to_predict <- reactive({
    if (input$stage != 'Real-time after enrollment completion') {
      input$to_predict
    } else {
      input$to_predict2
    }
  })


  # whether to show or hide enrollment and event panels
  observeEvent(to_predict(), {
    if (to_predict() == 'Enrollment only') {
      showTab(inputId = "results", target = "enroll_model_panel")
      hideTab(inputId = "results", target = "event_model_panel")
    } else if (to_predict() == 'Enrollment and event') {
      showTab(inputId = "results", target = "enroll_model_panel")
      showTab(inputId = "results", target = "event_model_panel")
    } else if (to_predict() == 'Event only') {
      hideTab(inputId = "results", target = "enroll_model_panel")
      showTab(inputId = "results", target = "event_model_panel")
    }
  })


  # input data set
  df <- reactive({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    inFile <- input$file1

    if (is.null(inFile))
      return(NULL)

    df <- readxl::read_excel(inFile$datapath)

    if (to_predict() == "Enrollment only") {
      req_cols <- c('trialsdt', 'usubjid', 'randdt', 'cutoffdt')
    } else {
      req_cols <- c('trialsdt', 'usubjid', 'randdt', 'cutoffdt',
                    'time', 'event', 'dropout')
    }

    if (input$by_treatment) {
      req_cols <- c(req_cols, 'treatment')
    }

    cols <- colnames(df)

    shiny::validate(
      need(all(req_cols %in% cols),
           paste("The following columns are missing from the input data:",
                 paste(req_cols[!(req_cols %in% cols)], collapse = ", "))))

    if (any(is.na(df[, req_cols]))) {
      stop(paste("The following columns have missing values:",
                 paste(req_cols[sapply(df, function(x) any(is.na(x)))],
                       collapse = ", ")))
    }

    if ('treatment' %in% cols && !('treatment_description' %in% cols)) {
      df <- df %>%
        mutate(treatment_description = paste("Treatment", treatment))
    }

    tibble(df) %>%
      mutate(trialsdt = as.Date(trialsdt),
             randdt = as.Date(randdt),
             cutoffdt = as.Date(cutoffdt))
  })


  # summarize observed data
  observed <- reactive({
    if (!is.null(df()))
      summarizeObserved(df(), to_predict(), showplot = FALSE,
                        input$by_treatment)
  })


  target_n <- reactive({
    req(input$target_n)
    valid = (input$target_n > 0 && input$target_n == round(input$target_n))
    shinyFeedback::feedbackWarning(
      "target_n", !valid,
      "Target enrollment must be a positive integer")
    req(valid)
    as.numeric(input$target_n)
  })


  target_d <- reactive({
    req(input$target_d)
    valid1 = (input$target_d > 0 && input$target_d == round(input$target_d))
    shinyFeedback::feedbackWarning(
      "target_d", !valid1,
      "Target events must be a positive integer")

    if (to_predict() == "Enrollment and event") {
      valid2 = (input$target_d <= input$target_n)
      shinyFeedback::feedbackWarning(
        "target_d", !valid2,
        "Target events must be less than or equal to target enrollment")
    } else {
      valid2 = (input$target_d <= observed()$n0)
      shinyFeedback::feedbackWarning(
        "target_d", !valid2,
        "Target events must be less than or equal to sample size")
    }

    req(valid1 && valid2)

    as.numeric(input$target_d)
  })


  pilevel <- reactive(as.numeric(input$pilevel))


  nyears <- reactive({
    req(input$nyears)
    valid = (input$nyears > 0)
    shinyFeedback::feedbackWarning(
      "nyears", !valid,
      "Years after cutoff must be a positive number")
    req(valid)
    as.numeric(input$nyears)
  })


  showEnrollment <- reactive({
    "Enrollment" %in% input$to_show
  })


  showEvent <- reactive({
    "Event" %in% input$to_show
  })


  showOngoing <- reactive({
    "Ongoing" %in% input$to_show
  })


  k <- reactive({
    if (!input$by_treatment && input$stage != "Design stage") {
      k = 1
    } else if (input$stage != "Design stage" && !is.null(df())) {
      k = length(table(df()$treatment))
      updateSelectInput(session, "k", selected=k)
    } else {
      k = as.numeric(input$k)
    }
    k
  })


  # whether to predict dosing
  predict_dosing <- reactive({
    if (input$by_treatment &&
        (to_predict() == 'Enrollment and event' ||
         input$stage == 'Real-time after enrollment completion')) {
      input$predict_dosing
    } else {
      FALSE
    }
  })


  # input data set
  visitview <- reactive({
    # input$file2 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    inFile <- input$file2

    if (is.null(inFile))
      return(NULL)

    df <- readxl::read_excel(inFile$datapath)

    req_cols <- c('usubjid', 'date', 'drug', 'drug_name', 'dose_unit',
                  'dispensed_quantity')

    cols <- colnames(df)

    shiny::validate(
      need(all(req_cols %in% cols),
           paste("The following columns are missing from the input data:",
                 paste(req_cols[!(req_cols %in% cols)], collapse = ", "))))

    if (any(is.na(df[, req_cols]))) {
      stop(paste("The following columns have missing values:",
                 paste(req_cols[sapply(df, function(x) any(is.na(x)))],
                       collapse = ", ")))
    }

    tibble(df) %>%
      mutate(date = as.Date(date))
  })


  dose_observed <- reactive({
    if (!is.null(df()) & !is.null(visitview())) {
      f_dose_observed(df(), visitview(), showplot = FALSE)
    } else {
      NULL
    }
  })


  # whether to make protocol based predictions only
  pred_pp_only <- reactive({
    if (predict_dosing()) {
      input$pred_pp_only
    } else {
      TRUE
    }
  })


  observeEvent(list(input$stage, to_predict(), predict_dosing()), {
    if (input$stage != 'Design stage' &&
        to_predict() != 'Enrollment only' &&
        predict_dosing()) {
      showTab(inputId = "observed_data_tabset", target = "dosing_data_panel")
    } else {
      hideTab(inputId = "observed_data_tabset", target = "dosing_data_panel")
    }
  })


  observeEvent(predict_dosing(), {
    if (predict_dosing()) {
      showTab(inputId = "results", target = "dosing_model_panel")
      showTab(inputId = "results", target = "dosing_prediction_panel")
    } else {
      hideTab(inputId = "results", target = "dosing_model_panel")
      hideTab(inputId = "results", target = "dosing_prediction_panel")
    }
  })


  showModelBased <- reactive({
    "Model based" %in% input$to_show_dosing
  })


  showProtocolBased <- reactive({
    "Protocol based" %in% input$to_show_dosing
  })


  treatment_allocation <- reactive({
    req(k())
    if (k() > 1) {
      d = input[[paste0("treatment_allocation_", k())]]
      d <- as.numeric(d)

      valid = all(d > 0 & d == round(d))
      if (!valid) {
        showNotification("Treatment allocation must be positive integers")
      }
      req(valid)
      d
    } else {
      1
    }
  })


  treatment_description <- reactive({
    req(k())
    if (k() > 1) {
      if (!input$by_treatment && input$stage != "Design stage") {
        a = "Overall"
      } else if (input$stage != "Design stage" && !is.null(df())) {
        treatment_mapping <- df() %>%
          group_by(treatment, treatment_description) %>%
          slice(n()) %>%
          select(treatment, treatment_description)
        a = treatment_mapping$treatment_description
      } else {
        a = rownames(input[[paste0("treatment_allocation_", k())]])
      }
    } else {
      a = "Overall"
    }
    a
  })


  observeEvent(treatment_description(), {
    if (input$by_treatment && !is.null(df())) {
      updateMatrixInput(
        session, paste0("treatment_allocation_", k()),
        value=matrix(treatment_allocation(), ncol = 1,
                     dimnames = list(treatment_description(), "Size")))
    }
  })


  nreps <- reactive({
    req(input$nreps)
    valid = (input$nreps > 0 && input$nreps == round(input$nreps))
    shinyFeedback::feedbackWarning(
      "nreps", !valid,
      "Number of simulations must be a positive integer")
    req(valid)
    as.numeric(input$nreps)
  })


  # observed data panel outputs
  output$dates <- renderText({
    if (!is.null(observed())) {
      str1 <- paste("Trial start date:", observed()$trialsdt)
      str2 <- paste("Data cutoff date:", observed()$cutoffdt)
      str3 <- paste("Days since trial start:", observed()$t0)
      paste(str1, str2, str3, sep='<br/>')
    }
  })


  output$statistics <- renderPrint({
    if (!is.null(df())) {

      if (input$by_treatment && k() > 1) {
        if (to_predict() == "Enrollment and event" ||
            to_predict() == "Event only") {
          sum_by_trt <- df() %>%
            bind_rows(df() %>% mutate(
              treatment = 9999, treatment_description = "Overall")) %>%
            group_by(treatment, treatment_description) %>%
            summarise(n0 = n(),
                      d0 = sum(event),
                      r0 = sum(!(event | dropout)),
                      rp = sum((time < as.numeric(
                        cutoffdt - randdt + 1)) & !(event | dropout)),
                      .groups = "drop")

          if (any(sum_by_trt$rp) > 0) {
            table <- t(sum_by_trt %>% select(n0, d0, r0, rp))
            colnames(table) <- sum_by_trt$treatment_description
            rownames(table) <- c("Current number of subjects",
                                 "Current number of events",
                                 "Number of ongoing subjects",
                                 "  With ongoing date before cutoff")
          } else {
            table <- t(sum_by_trt %>% select(n0, d0, r0))
            colnames(table) <- sum_by_trt$treatment_description
            rownames(table) <- c("Current number of subjects",
                                 "Current number of events",
                                 "Number of ongoing subjects")
          }
        } else {
          sum_by_trt <- df() %>%
            bind_rows(df() %>% mutate(
              treatment = 9999, treatment_description = "Overall")) %>%
            group_by(treatment, treatment_description) %>%
            summarise(n0 = n(), .groups = "drop")

          table <- t(sum_by_trt %>% select(n0))
          colnames(table) <- sum_by_trt$treatment_description
          rownames(table) <- c("Current number of subjects")
        }
      } else {
        if (to_predict() == "Enrollment and event" ||
            to_predict() == "Event only") {
          sum_overall <- tibble(n0 = observed()$n0,
                                d0 = observed()$d0,
                                r0 = observed()$r0,
                                rp = observed()$rp)

          if (sum_overall$rp > 0) {
            table <- t(sum_overall %>% select(n0, d0, r0, rp))
            colnames(table) <- "Overall"
            rownames(table) <- c("Current number of subjects",
                                 "Current number of events",
                                 "Number of ongoing subjects",
                                 "  With ongoing date before cutoff")
          } else {
            table <- t(sum_overall %>% select(n0, d0, r0))
            colnames(table) <- "Overall"
            rownames(table) <- c("Current number of subjects",
                                 "Current number of events",
                                 "Number of ongoing subjects")
          }
        } else {
          table <- t(tibble(n0 = observed()$n0))
          colnames(table) <- "Overall"
          rownames(table) <- c("Current number of subjects")
        }
      }

      print(table, quote = FALSE)
    }
  })


  output$cum_accrual_plot <- renderPlotly({
    cum_accrual_plot <- observed()$cum_accrual_plot
    if (!is.null(cum_accrual_plot)) cum_accrual_plot
  })


  output$daily_accrual_plot <- renderPlotly({
    daily_accrual_plot <- observed()$daily_accrual_plot
    if (!is.null(daily_accrual_plot)) daily_accrual_plot
  })


  output$event_km_plot <- renderPlotly({
    event_km_plot <- observed()$event_km_plot
    if (!is.null(event_km_plot)) event_km_plot
  })


  output$input_df <- renderDataTable(
    df(), options = list(pageLength = 10)
  )


  output$dates_copy <- renderText({
    if (!is.null(observed())) {
      str1 <- paste("Trial start date:", observed()$trialsdt)
      str2 <- paste("Data cutoff date:", observed()$cutoffdt)
      str3 <- paste("Days since trial start:", observed()$t0)
      paste(str1, str2, str3, sep='<br/>')
    }
  })


  output$cum_dose_t0 <- renderPrint({
    if (!is.null(dose_observed())) {
      table <- as.data.frame(dose_observed()$dosing_summary_t0)

      colnames(table) <- c("Drug", "Drug Name", "Dose Unit",
                           "Cumulative Dose")
      print(table, quote = FALSE, row.names = FALSE)
    }
  })


  output$cum_dispense_plot <- renderPlotly({
    if (!is.null(dose_observed())) {
      dose_observed()$cum_dispense_plot
    }
  })


  output$bar_t0_plot <- renderPlotly({
    if (!is.null(dose_observed())) {
      dose_observed()$bar_t0_plot
    }
  })


  output$bar_ti_plot <- renderPlotly({
    if (!is.null(dose_observed())) {
      dose_observed()$bar_ti_plot
    }
  })


  output$bar_di_plot <- renderPlotly({
    if (!is.null(dose_observed())) {
      dose_observed()$bar_di_plot
    }
  })


  output$input_visitview <- renderDataTable(
    visitview(), options = list(pageLength = 10)
  )


  # enrollment panel outputs
  poisson_rate <- reactive({
    req(input$poisson_rate)
    valid = (input$poisson_rate > 0)
    shinyFeedback::feedbackWarning(
      "poisson_rate", !valid,
      "Daily enrollment rate must be a positive number")
    req(valid)
    as.numeric(input$poisson_rate)
  })


  mu <- reactive({
    req(input$mu)
    valid = (input$mu > 0)
    shinyFeedback::feedbackWarning(
      "mu", !valid,
      "Base rate must be a positive number")
    req(valid)
    as.numeric(input$mu)
  })


  delta <- reactive({
    req(input$delta)
    valid = (input$delta > 0)
    shinyFeedback::feedbackWarning(
      "delta", !valid,
      "Decay rate must be a positive number")
    req(valid)
    as.numeric(input$delta)
  })


  piecewise_poisson_rate <- reactive({
    req(input$piecewise_poisson_rate)
    t = as.numeric(input$piecewise_poisson_rate[,1])
    lambda = as.numeric(input$piecewise_poisson_rate[,2])

    valid1 = all(diff(t) > 0) && (t[1] == 0)
    if (!valid1) {
      showNotification(
        "Starting time must be increasing and start at zero"
      )
    }

    valid2 = all(lambda >= 0)
    if (!valid2) {
      showNotification(
        "Enrollment rate must be nonnegative"
      )
    }

    valid3 = any(lambda > 0)
    if (!valid3) {
      showNotification(
        "At least one enrollment rate must be positive"
      )
    }

    req(valid1 && valid2 && valid3)

    matrix(c(t, lambda), ncol = 2,
           dimnames = list(paste("Interval", 1:length(t)),
                           c("Starting time", "Enrollment rate")))
  })


  observeEvent(input$add_piecewise_poisson_rate, {
    a = matrix(as.numeric(input$piecewise_poisson_rate),
               ncol=ncol(input$piecewise_poisson_rate))
    b = matrix(a[nrow(a),], nrow=1)
    b[1,1] = b[1,1] + 90
    c = rbind(a, b)
    rownames(c) = paste("Interval", seq(1,nrow(c)))
    colnames(c) = colnames(input$piecewise_poisson_rate)
    updateMatrixInput(session, "piecewise_poisson_rate", c)
  })


  observeEvent(input$del_piecewise_poisson_rate, {
    if (nrow(input$piecewise_poisson_rate) >= 2) {
      a = matrix(as.numeric(input$piecewise_poisson_rate),
                 ncol=ncol(input$piecewise_poisson_rate))
      b = matrix(a[-nrow(a),], ncol=ncol(a))
      rownames(b) = paste("Interval", seq(1,nrow(b)))
      colnames(b) = colnames(input$piecewise_poisson_rate)
      updateMatrixInput(session, "piecewise_poisson_rate", b)
    }
  })


  nknots <- reactive({
    req(input$nknots)
    valid = (input$nknots >= 0 && input$nknots == round(input$nknots))
    shinyFeedback::feedbackWarning(
      "nknots", !valid,
      "Number of inner knots must be a nonnegative integer")
    req(valid)
    as.numeric(input$nknots)
  })


  lags <- reactive({
    req(input$lags)
    valid = (input$lags >= 0 && input$lags == round(input$lags))
    shinyFeedback::feedbackWarning(
      "lags", !valid,
      "Number of day lags must be a nonnegative integer")
    req(valid)
    as.numeric(input$lags)
  })


  accrualTime <- reactive({
    t = as.numeric(input$accrualTime)
    valid = all(diff(t) > 0) && (t[1] == 0)
    if (!valid) {
      showNotification(
        "Starting time must be increasing and start at zero"
      )
    }
    req(valid)
    t
  })


  observeEvent(input$add_accrualTime, {
    a = matrix(as.numeric(input$accrualTime),
               ncol=ncol(input$accrualTime))
    b = matrix(a[nrow(a),] + 90, nrow=1)
    c = rbind(a, b)
    rownames(c) = paste("Interval", seq(1,nrow(c)))
    colnames(c) = colnames(input$accrualTime)
    updateMatrixInput(session, "accrualTime", c)
  })


  observeEvent(input$del_accrualTime, {
    if (nrow(input$accrualTime) >= 2) {
      a = matrix(as.numeric(input$accrualTime),
                 ncol=ncol(input$accrualTime))
      b = matrix(a[-nrow(a),], ncol=ncol(a))
      rownames(b) = paste("Interval", seq(1,nrow(b)))
      colnames(b) = colnames(input$accrualTime)
      updateMatrixInput(session, "accrualTime", b)
    }
  })


  enroll_fit <- reactive({
    if (!is.null(df()))
      fitEnrollment(df(), input$enroll_model, nknots(),
                    accrualTime(), showplot = FALSE)
  })


  output$enroll_fit <- renderPlotly({
    if (!is.null(enroll_fit())) enroll_fit()$fit_plot
  })


  # event panel outputs
  exponential_survival <- reactive({
    req(k())
    param = input[[paste0("exponential_survival_", k())]]
    lambda = as.numeric(param)
    valid = all(lambda > 0)
    if (!valid) {
      showNotification(
        "Hazard rate must be positive"
      )
    }
    req(valid)
    lambda
  })


  weibull_survival <- reactive({
    req(k())
    param = input[[paste0("weibull_survival_", k())]]
    shape = as.numeric(param[1,])
    scale = as.numeric(param[2,])

    valid1 = all(shape > 0)
    if (!valid1) {
      showNotification(
        "Weibull shape parameter must be positive"
      )
    }

    valid2 = all(scale > 0)
    if (!valid2) {
      showNotification(
        "Weibull scale parameter must be positive"
      )
    }

    req(valid1 && valid2)

    matrix(c(shape, scale), nrow = 2, byrow = TRUE)
  })


  llogis_survival <- reactive({
    req(k())
    param = input[[paste0("llogis_survival_", k())]]
    locationlog = as.numeric(param[1,])
    scalelog = as.numeric(param[2,])

    valid = all(scalelog > 0)
    if (!valid) {
      showNotification(
        "Scale on the log scale must be positive"
      )
    }

    req(valid)

    matrix(c(locationlog, scalelog), nrow = 2, byrow = TRUE)
  })


  lnorm_survival <- reactive({
    req(k())
    param = input[[paste0("lnorm_survival_", k())]]
    meanlog = as.numeric(param[1,])
    sdlog = as.numeric(param[2,])

    valid = all(sdlog > 0)
    if (!valid) {
      showNotification(
        "SD on the log scale must be positive"
      )
    }

    req(valid)

    matrix(c(meanlog, sdlog), nrow = 2, byrow = TRUE)
  })


  piecewise_exponential_survival <- reactive({
    req(k())
    param = input[[paste0("piecewise_exponential_survival_", k())]]
    t = as.numeric(param[,1])
    lambda = as.numeric(param[,-1])

    valid1 = all(diff(t) > 0) && (t[1] == 0)
    if (!valid1) {
      showNotification(
        "Starting time must be increasing and start at zero"
      )
    }

    valid2 = all(lambda > 0)
    if (!valid2) {
      showNotification(
        "Hazard rate must be positive"
      )
    }

    req(valid1 && valid2)

    matrix(c(t, lambda), nrow = length(t))
  })


  lapply(1:6, function(i) {
    pwexp <- paste0("piecewise_exponential_survival_", i)
    observeEvent(input[[paste0("add_piecewise_exponential_survival_", i)]], {
      a = matrix(as.numeric(input[[pwexp]]), ncol=ncol(input[[pwexp]]))
      b = matrix(a[nrow(a),], nrow=1)
      b[1,1] = b[1,1] + 90
      c = rbind(a, b)
      rownames(c) = paste("Interval", seq(1,nrow(c)))
      colnames(c) = colnames(input[[pwexp]])
      updateMatrixInput(session, pwexp, c)
    })
  })


  lapply(1:6, function(i) {
    pwexp <- paste0("piecewise_exponential_survival_", i)
    observeEvent(input[[paste0("del_piecewise_exponential_survival_", i)]], {
      if (nrow(input[[pwexp]]) >= 2) {
        a = matrix(as.numeric(input[[pwexp]]), ncol=ncol(input[[pwexp]]))
        b = matrix(a[-nrow(a),], ncol=ncol(a))
        rownames(b) = paste("Interval", seq(1,nrow(b)))
        colnames(b) = colnames(input[[pwexp]])
        updateMatrixInput(session, pwexp, b)
      }
    })
  })


  observeEvent(treatment_description(), {
    if (input$stage == "Design stage") {
      updateMatrixInput(
        session, paste0("exponential_survival_", k()),
        value=matrix(exponential_survival(), ncol=k(),
                     dimnames = list(NULL, treatment_description())))
      updateMatrixInput(
        session, paste0("weibull_survival_", k()),
        value=matrix(weibull_survival(), nrow=2, ncol=k(),
                     dimnames = list(c("Shape", "Scale"),
                                     treatment_description())))
      updateMatrixInput(
        session, paste0("llogis_survival_", k()),
        value=matrix(llogis_survival(), nrow=2, ncol=k(),
                     dimnames = list(c("Location on log scale",
                                       "Scale on log scale"),
                                     treatment_description())))
      updateMatrixInput(
        session, paste0("lnorm_survival_", k()),
        value=matrix(lnorm_survival(), nrow=2, ncol=k(),
                     dimnames = list(c("Mean on log scale",
                                       "SD on log scale"),
                                     treatment_description())))
      npieces = nrow(piecewise_exponential_survival())
      updateMatrixInput(
        session, paste0("piecewise_exponential_survival_", k()),
        value=matrix(piecewise_exponential_survival(),
                     nrow=npieces, ncol=k()+1,
                     dimnames = list(
                       paste("Interval", seq_len(npieces)),
                       c("Starting time", treatment_description()))))
    }
  })


  piecewiseSurvivalTime <- reactive({
    t = as.numeric(input$piecewiseSurvivalTime)
    valid = all(diff(t) > 0) && (t[1] == 0)
    if (!valid) {
      showNotification(
        "Starting time must be increasing and start at zero"
      )
    }
    req(valid)
    t
  })


  observeEvent(input$add_piecewiseSurvivalTime, {
    a = matrix(as.numeric(input$piecewiseSurvivalTime),
               ncol=ncol(input$piecewiseSurvivalTime))
    b = matrix(a[nrow(a),] + 90, nrow=1)
    c = rbind(a, b)
    rownames(c) = paste("Interval", seq(1,nrow(c)))
    colnames(c) = colnames(input$piecewiseSurvivalTime)
    updateMatrixInput(session, "piecewiseSurvivalTime", c)
  })


  observeEvent(input$del_piecewiseSurvivalTime, {
    if (nrow(input$piecewiseSurvivalTime) >= 2) {
      a = matrix(as.numeric(input$piecewiseSurvivalTime),
                 ncol=ncol(input$piecewiseSurvivalTime))
      b = matrix(a[-nrow(a),], ncol=ncol(a))
      rownames(b) = paste("Interval", seq(1,nrow(b)))
      colnames(b) = colnames(input$piecewiseSurvivalTime)
      updateMatrixInput(session, "piecewiseSurvivalTime", b)
    }
  })


  spline_k <- reactive({
    req(input$spline_k)
    valid = (input$spline_k >= 0 && input$spline_k == round(input$spline_k))
    shinyFeedback::feedbackWarning(
      "spline_k", !valid,
      "Number of inner knots must be a nonnegative integer")
    req(valid)
    as.numeric(input$spline_k)
  })


  event_fit <- reactive({
    if (!is.null(df()))
      fitEvent(df(), input$event_model, piecewiseSurvivalTime(),
               spline_k(), input$spline_scale, showplot = FALSE,
               input$by_treatment)
  })


  output$event_fit_ic <- renderText({
    if (input$by_treatment && k() > 1 && !is.null(event_fit())) {
      aic = sum(sapply(event_fit(), function(fit) fit$fit$aic))
      bic = sum(sapply(event_fit(), function(fit) fit$fit$bic))
      aictext = paste("Total AIC:", formatC(aic, format = "f", digits = 2))
      bictext = paste("Total BIC:", formatC(bic, format = "f", digits = 2))
      text1 = paste0("<i>", aictext, ", ", bictext, "</i>")
    } else {
      text1 = NULL
    }

    if (!is.null(text1)) text1
  })


  observe({
    walk(1:6, function(i) {
      output[[paste0("event_fit_output", i)]] <- renderPlotly({
        if (i <= k() && !is.null(event_fit())) {
          if (input$by_treatment && k() > 1) {
            event_fit()[[i]]$fit_plot
          } else {
            event_fit()$fit_plot
          }
        } else {
          NULL
        }
      })
    })
  })


  event_fit_outputs <- reactive({
    outputs <- map(1:k(), function(i) {
      plotlyOutput(paste0("event_fit_output", i))
    })

    tagList(outputs)
  })

  output$event_fit <- renderUI({
    event_fit_outputs()
  })


  # enrollment and event prediction panel outputs
  pred <- eventReactive(input$predict, {
    set.seed(as.numeric(input$seed))

    if (to_predict() != "Enrollment only") {
      shiny::validate(
        need(showEnrollment() || showEvent() || showOngoing(),
             "Need at least one parameter to show on prediction plot"))
    }

    if (input$stage == "Design stage") {
      w = treatment_allocation()/sum(treatment_allocation())

      # enroll model specifications
      if (input$enroll_prior == "Poisson") {
        theta = log(poisson_rate())
      } else if (input$enroll_prior == "Time-decay") {
        theta = c(log(mu()), log(delta()))
      } else if (input$enroll_prior == "Piecewise Poisson") {
        theta = log(piecewise_poisson_rate()[,2])
        accrualTime = piecewise_poisson_rate()[,1]
      }

      enroll_prior <- list(
        model = input$enroll_prior,
        theta = theta,
        vtheta = diag(length(theta))*1e-8)

      if (input$enroll_prior == "Piecewise Poisson") {
        enroll_prior$accrualTime = accrualTime
      }

      # event model specifications
      if (to_predict() == "Enrollment and event") {
        model = input$event_prior
        event_prior <- list()

        for (i in 1:k()) {
          if (model == "Exponential") {
            theta = log(exponential_survival()[i])
          } else if (model == "Weibull") {
            theta = c(log(weibull_survival()[2,i]),
                      -log(weibull_survival()[1,i]))
          } else if (model == "Log-logistic") {
            theta = c(llogis_survival()[1,i], log(llogis_survival()[2,i]))
          } else if (model == "Log-normal") {
            theta = c(lnorm_survival()[1,i], log(lnorm_survival()[2,i]))
          } else if (model == "Piecewise exponential") {
            theta = log(piecewise_exponential_survival()[,i+1])
            piecewiseSurvivalTime = piecewise_exponential_survival()[,1]
          }

          if (model != "Piecewise exponential") {
            event_prior[[i]] <- list(
              model = model,
              theta = theta,
              vtheta = diag(length(theta))*1e-8,
              w = w[i])
          } else {
            event_prior[[i]] <- list(
              model = model,
              theta = theta,
              vtheta = diag(length(theta))*1e-8,
              piecewiseSurvivalTime = piecewiseSurvivalTime,
              w = w[i])
          }
        }

        if (k() == 1) event_prior <- event_prior[[1]]
      }

      # get prediction results based on what to predict
      if (to_predict() == "Enrollment only") {
        getPrediction(
          to_predict = to_predict(),
          target_n = target_n(),
          enroll_prior = enroll_prior,
          pilevel = pilevel(),
          nyears = nyears(),
          nreps = nreps(),
          showsummary = FALSE,
          showplot = FALSE,
          by_treatment = input$by_treatment,
          ngroups = k(),
          alloc = treatment_allocation(),
          treatment_label = treatment_description())
      } else if (to_predict() == "Enrollment and event") {
        getPrediction(
          to_predict = to_predict(),
          target_n = target_n(),
          target_d = target_d(),
          enroll_prior = enroll_prior,
          event_prior = event_prior,
          pilevel = pilevel(),
          nyears = nyears(),
          nreps = nreps(),
          showEnrollment = showEnrollment(),
          showEvent = showEvent(),
          showOngoing = showOngoing(),
          showsummary = FALSE,
          showplot = FALSE,
          by_treatment = input$by_treatment,
          ngroups = k(),
          alloc = treatment_allocation(),
          treatment_label = treatment_description())
      }
    } else { # real-time prediction
      shiny::validate(
        need(!is.null(df()),
             paste("Please upload enrollment and event data for",
                   "real-time prediction.")))

      if (to_predict() == "Enrollment only") {
        shiny::validate(
          need(target_n() > observed()$n0,
               "Target enrollment has been reached."))

        getPrediction(
          df = df(),
          to_predict = to_predict(),
          target_n = target_n(),
          enroll_model = input$enroll_model,
          nknots = nknots(),
          lags = lags(),
          accrualTime = accrualTime(),
          pilevel = pilevel(),
          nyears = nyears(),
          nreps = nreps(),
          showsummary = FALSE,
          showplot = FALSE,
          by_treatment = input$by_treatment,
          alloc = treatment_allocation())
      } else if (to_predict() == "Enrollment and event") {
        shiny::validate(
          need(target_n() > observed()$n0,
               "Target enrollment has been reached."))

        shiny::validate(
          need(target_d() > observed()$d0,
               "Target number of events has been reached."))

        getPrediction(
          df = df(),
          to_predict = to_predict(),
          target_n = target_n(),
          target_d = target_d(),
          enroll_model = input$enroll_model,
          nknots = nknots(),
          lags = lags(),
          accrualTime = accrualTime(),
          event_model = input$event_model,
          piecewiseSurvivalTime = piecewiseSurvivalTime(),
          k = spline_k(),
          scale = input$spline_scale,
          dropout_model = "none",
          pilevel = pilevel(),
          nyears = nyears(),
          nreps = nreps(),
          showEnrollment = showEnrollment(),
          showEvent = showEvent(),
          showOngoing = showOngoing(),
          showsummary = FALSE,
          showplot = FALSE,
          by_treatment = input$by_treatment,
          alloc = treatment_allocation())
      } else if (to_predict() == "Event only") {
        shiny::validate(
          need(target_d() > observed()$d0,
               "Target number of events has been reached."))

        getPrediction(
          df = df(),
          to_predict = to_predict(),
          target_d = target_d(),
          event_model = input$event_model,
          piecewiseSurvivalTime = piecewiseSurvivalTime(),
          k = spline_k(),
          scale = input$spline_scale,
          dropout_model = "none",
          pilevel = pilevel(),
          nyears = nyears(),
          nreps = nreps(),
          showEnrollment = showEnrollment(),
          showEvent = showEvent(),
          showOngoing = showOngoing(),
          showsummary = FALSE,
          showplot = FALSE,
          by_treatment = input$by_treatment)
      }
    }
  })


  output$enroll_pred_date <- renderText({
    if (to_predict() == 'Enrollment only' ||
        to_predict() == 'Enrollment and event') {

      req(pred()$enroll_pred)
      req(pred()$stage == input$stage && pred()$to_predict == to_predict())


      if (input$stage != 'Design stage') {

        shiny::validate(
          need(!is.null(df()),
               "Please upload data for real-time prediction."))

        shiny::validate(
          need(target_n() > observed()$n0,
               "Target enrollment has been reached."))

        if (!is.null(pred()$enroll_pred$enroll_pred_date)) {
          str1 <- paste0("Time from cutoff until ",
                         pred()$enroll_pred$target_n, " subjects: ",
                         pred()$enroll_pred$enroll_pred_date[1] -
                           observed()$cutoffdt + 1, " days")
          str2 <- paste0("Median prediction date: ",
                         pred()$enroll_pred$enroll_pred_date[1])
          str3 <- paste0("Prediction interval: ",
                         pred()$enroll_pred$enroll_pred_date[2], ", ",
                         pred()$enroll_pred$enroll_pred_date[3])
          text1 <- paste(paste('<b>', str1, '</b>'), str2, str3, sep='<br/>')
        } else {
          text1 <- NULL
        }
      } else {
        if (!is.null(pred()$enroll_pred$enroll_pred_day)) {
          str1 <- paste0("Time from trial start until ",
                         pred()$enroll_pred$target_n, " subjects")
          str2 <- paste0("Median prediction day: ",
                         pred()$enroll_pred$enroll_pred_day[1])
          str3 <- paste0("Prediction interval: ",
                         pred()$enroll_pred$enroll_pred_day[2], ", ",
                         pred()$enroll_pred$enroll_pred_day[3])
          text1 <- paste(paste('<b>', str1, '</b>'), str2, str3, sep='<br/>')
        } else {
          text1 <- NULL
        }
      }
    } else {
      text1 <- NULL
    }

    if (!is.null(text1)) text1
  })


  output$event_pred_date <- renderText({
    if (to_predict() == 'Enrollment and event' ||
        to_predict() == 'Event only') {

      req(pred()$event_pred)
      req(pred()$stage == input$stage && pred()$to_predict == to_predict())

      if (input$stage != 'Design stage') {
        shiny::validate(
          need(!is.null(df()),
               "Please upload data for real-time prediction."))

        shiny::validate(
          need(target_d() > observed()$d0,
               "Target number of events has been reached."))


        if (!is.null(pred()$event_pred$event_pred_date)) {
          str1 <- paste0("Time from cutoff until ",
                         pred()$event_pred$target_d, " events: ",
                         pred()$event_pred$event_pred_date[1] -
                           observed()$cutoffdt + 1, " days")
          str2 <- paste0("Median prediction date: ",
                         pred()$event_pred$event_pred_date[1])
          str3 <- paste0("Prediction interval: ",
                         pred()$event_pred$event_pred_date[2], ", ",
                         pred()$event_pred$event_pred_date[3])
          text2 <- paste(paste('<b>', str1, '</b>'), str2, str3, sep='<br/>')
        } else {
          text2 <- NULL
        }
      } else {
        if (!is.null(pred()$event_pred$event_pred_day)) {
          str1 <- paste0("Time from trial start until ",
                         pred()$event_pred$target_d, " events")
          str2 <- paste0("Median prediction day: ",
                         pred()$event_pred$event_pred_day[1])
          str3 <- paste0("Prediction interval: ",
                         pred()$event_pred$event_pred_day[2], ", ",
                         pred()$event_pred$event_pred_day[3])
          text2 <- paste(paste('<b>', str1, '</b>'), str2, str3, sep='<br/>')
        } else {
          text2 <- NULL
        }
      }
    } else {
      text2 <- NULL
    }

    if (!is.null(text2)) text2
  })


  output$pred_date <- renderUI({
    if (to_predict() == 'Enrollment only') {
      htmlOutput("enroll_pred_date")
    } else if (to_predict() == 'Event only') {
      htmlOutput("event_pred_date")
    } else {
      fluidRow(column(6, htmlOutput("enroll_pred_date")),
               column(6, htmlOutput("event_pred_date")))
    }
  })


  # enrollment and event prediction plot
  pred_plot <- reactive({
    if (to_predict() == "Enrollment only") {
      req(pred()$enroll_pred)
      req(pred()$stage == input$stage && pred()$to_predict == to_predict())

      if (input$stage != 'Design stage') {
        shiny::validate(
          need(!is.null(df()),
               paste("Please upload enrollment and event data",
                     "for real-time prediction.")))

        shiny::validate(
          need(target_n() > observed()$n0,
               "Target enrollment has been reached."))
      }

      enroll_pred_plot <- pred()$enroll_pred$enroll_pred_plot
      enroll_pred_df <- pred()$enroll_pred$enroll_pred_df
      if ((!input$by_treatment || k() == 1) ||
          ((input$by_treatment || input$stage == 'Design stage') &&
           k() > 1 && "treatment" %in% names(enroll_pred_df) &&
           length(table(enroll_pred_df$treatment)) == k() + 1)) {
        g <- enroll_pred_plot
      } else {
        g <- NULL
      }
    } else { # predict event only or predict enrollment and event
      shiny::validate(
        need(showEnrollment() || showEvent() || showOngoing(),
             "Need at least one parameter to show on event prediction plot"))

      req(pred()$event_pred)
      req(pred()$stage == input$stage && pred()$to_predict == to_predict())


      if (input$stage != 'Design stage') {
        shiny::validate(
          need(!is.null(df()),
               paste("Please upload enrollment and event data",
                     "for real-time prediction.")))

        if (to_predict() == "Enrollment and event")
          shiny::validate(
            need(target_n() > observed()$n0,
                 "Target enrollment has been reached."))

        shiny::validate(
          need(target_d() > observed()$d0,
               "Target number of events has been reached."))
      }


      dfs <- tibble()

      if (showEnrollment())
        dfs <- dfs %>% bind_rows(pred()$event_pred$enroll_pred_df)

      if (showEvent())
        dfs <- dfs %>% bind_rows(pred()$event_pred$event_pred_df)

      if (showOngoing())
        dfs <- dfs %>% bind_rows(pred()$event_pred$ongoing_pred_df)


      if ((!input$by_treatment || k() == 1) &&
          !("treatment" %in% names(dfs))) { # overall
        if (input$stage != 'Design stage') {
          dfa <- dfs %>% filter(is.na(lower))
          dfb <- dfs %>% filter(!is.na(lower))

          dfa_enrollment <- dfa %>% filter(parameter == "Enrollment")
          dfb_enrollment <- dfb %>% filter(parameter == "Enrollment")
          dfa_event <- dfa %>% filter(parameter == "Event")
          dfb_event <- dfb %>% filter(parameter == "Event")
          dfa_ongoing <- dfa %>% filter(parameter == "Ongoing")
          dfb_ongoing <- dfb %>% filter(parameter == "Ongoing")

          g <- plotly::plot_ly() %>%
            plotly::add_lines(
              data = dfa_enrollment, x = ~date, y = ~n,
              line = list(shape="hv", width=2),
              name = "observed enrollment") %>%
            plotly::add_lines(
              data = dfb_enrollment, x = ~date, y = ~n,
              line = list(width=2),
              name = "median prediction enrollment") %>%
            plotly::add_ribbons(
              data = dfb_enrollment, x = ~date, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval enrollment") %>%
            plotly::add_lines(
              data = dfa_event, x = ~date, y = ~n,
              line = list(shape="hv", width=2),
              name = "observed_event") %>%
            plotly::add_lines(
              data = dfb_event, x = ~date, y = ~n,
              line = list(width=2),
              name = "median prediction event") %>%
            plotly::add_ribbons(
              data = dfb_event, x = ~date, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval event") %>%
            plotly::add_lines(
              data = dfa_ongoing, x = ~date, y = ~n,
              line = list(shape="hv", width=2),
              name = "observed ongoing") %>%
            plotly::add_lines(
              data = dfb_ongoing, x = ~date, y = ~n,
              line = list(width=2),
              name = "median prediction ongoing") %>%
            plotly::add_ribbons(
              data = dfb_ongoing, x = ~date, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval ongoing") %>%
            plotly::add_lines(
              x = rep(observed()$cutoffdt, 2),
              y = c(min(dfa$n), max(dfb$upper)),
              name = "cutoff", line = list(dash="dash"),
              showlegend = FALSE) %>%
            plotly::layout(
              annotations = list(
                x = observed()$cutoffdt, y = 0, text = 'cutoff',
                xanchor = "left", yanchor = "bottom", font = list(size = 12),
                showarrow = FALSE),
              xaxis = list(title = "", zeroline = FALSE),
              yaxis = list(zeroline = FALSE))

          if (observed()$tp < observed()$t0) {
            g <- g %>%
              plotly::add_lines(
                x = rep(observed()$cutofftpdt, 2),
                y = c(min(dfa$n), max(dfb$upper)),
                name = "prediction start",
                line = list(dash="dash", color="grey"),
                showlegend = FALSE) %>%
              plotly::layout(
                annotations = list(
                  x = observed()$cutofftpdt, y = 0,
                  text = 'prediction start',
                  xanchor = "left", yanchor = "bottom",
                  font = list(size=12), showarrow = FALSE))
          }

          if (showEvent()) {
            g <- g %>%
              plotly::add_lines(
                x = range(dfs$date), y = rep(target_d(), 2),
                name = 'target events', showlegend = FALSE,
                line = list(dash="dot", color="rgba(128, 128, 128, 0.5")) %>%
              plotly::layout(
                annotations = list(
                  x = 0.95, xref = "paper", y = target_d(),
                  text = 'target events', xanchor = "right",
                  yanchor = "bottom", font = list(size = 12),
                  showarrow = FALSE))
          }
        } else { # Design stage
          dfs_enrollment <- dfs %>% filter(parameter == "Enrollment")
          dfs_event <- dfs %>% filter(parameter == "Event")
          dfs_ongoing <- dfs %>% filter(parameter == "Ongoing")

          g <- plotly::plot_ly() %>%
            plotly::add_lines(
              data = dfs_enrollment, x = ~t, y = ~n,
              line = list(width=2),
              name = "median prediction enrollment") %>%
            plotly::add_ribbons(
              data = dfs_enrollment, x = ~t, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval enrollment") %>%
            plotly::add_lines(
              data = dfs_event, x = ~t, y = ~n,
              line = list(width=2),
              name = "median prediction event") %>%
            plotly::add_ribbons(
              data = dfs_event, x = ~t, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval event") %>%
            plotly::add_lines(
              data = dfs_ongoing, x = ~t, y = ~n,
              line = list(width=2),
              name = "median prediction ongoing") %>%
            plotly::add_ribbons(
              data = dfs_ongoing, x = ~t, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width=0),
              name = "prediction interval ongoing") %>%
            plotly::layout(
              xaxis = list(title = "Days since trial start",
                           zeroline = FALSE),
              yaxis = list(zeroline = FALSE))

          if (showEvent()) {
            g <- g %>%
              plotly::add_lines(
                x = range(dfs$t), y = rep(target_d(), 2),
                name = 'target events', showlegend = FALSE,
                line = list(dash="dot", color="rgba(128, 128, 128, 0.5")) %>%
              plotly::layout(
                annotations = list(
                  x = 0.95, xref = "paper", y = target_d(),
                  text = 'target events', xanchor = "right",
                  yanchor = "bottom", font = list(size = 12),
                  showarrow = FALSE))
          }
        }
      } else if (((input$by_treatment || input$stage == 'Design stage') &&
                  k() > 1) && ("treatment" %in% names(dfs)) &&
                 (length(table(dfs$treatment)) == k() + 1)) { # by treatment

        if (input$stage != 'Design stage') {
          dfa <- dfs %>% filter(is.na(lower))
          dfb <- dfs %>% filter(!is.na(lower))

          g <- list()
          for (i in c(9999, 1:k())) {
            dfsi <- dfs %>% filter(treatment == i)
            dfbi <- dfb %>% filter(treatment == i)
            dfai <- dfa %>% filter(treatment == i)

            dfai_enrollment <- dfai %>% filter(parameter == "Enrollment")
            dfbi_enrollment <- dfbi %>% filter(parameter == "Enrollment")
            dfai_event <- dfai %>% filter(parameter == "Event")
            dfbi_event <- dfbi %>% filter(parameter == "Event")
            dfai_ongoing <- dfai %>% filter(parameter == "Ongoing")
            dfbi_ongoing <- dfbi %>% filter(parameter == "Ongoing")

            g[[(i+1) %% 9999]] <- plotly::plot_ly() %>%
              plotly::add_lines(
                data = dfai_enrollment, x = ~date, y = ~n,
                line = list(shape="hv", width=2),
                name = "observed enrollment") %>%
              plotly::add_lines(
                data = dfbi_enrollment, x = ~date, y = ~n,
                line = list(width=2),
                name = "median prediction enrollment") %>%
              plotly::add_ribbons(
                data = dfbi_enrollment, x = ~date,
                ymin = ~lower, ymax = ~upper,
                fill = "tonexty", line = list(width=0),
                name = "prediction interval enrollment") %>%
              plotly::add_lines(
                data = dfai_event, x = ~date, y = ~n,
                line = list(shape="hv", width=2),
                name = "observed event") %>%
              plotly::add_lines(
                data = dfbi_event, x = ~date, y = ~n,
                line = list(width=2),
                name = "median prediction event") %>%
              plotly::add_ribbons(
                data = dfbi_event, x = ~date, ymin = ~lower, ymax = ~upper,
                fill = "tonexty", line = list(width=0),
                name = "prediction interval event") %>%
              plotly::add_lines(
                data = dfai_ongoing, x = ~date, y = ~n,
                line = list(shape="hv", width=2),
                name = "observed ongoing") %>%
              plotly::add_lines(
                data = dfbi_ongoing, x = ~date, y = ~n,
                line = list(width=2),
                name = "median prediction ongoing") %>%
              plotly::add_ribbons(
                data = dfbi_ongoing, x = ~date, ymin = ~lower, ymax = ~upper,
                fill = "tonexty", line = list(width=0),
                name = "prediction interval ongoing") %>%
              plotly::add_lines(
                x = rep(observed()$cutoffdt, 2),
                y = c(min(dfai$n), max(dfbi$upper)),
                name = "cutoff", line = list(dash="dash"),
                showlegend = FALSE) %>%
              plotly::layout(
                xaxis = list(title = "", zeroline = FALSE),
                yaxis = list(zeroline = FALSE)) %>%
              plotly::layout(
                annotations = list(
                  x = 0.5, y = 1,
                  text = paste0("<b>", dfsi$treatment_description[1], "</b>"),
                  xanchor = "center", yanchor = "bottom",
                  showarrow = FALSE, xref='paper', yref='paper'))


            if (observed()$tp < observed()$t0) {
              g[[(i+1) %% 9999]] <- g[[(i+1) %% 9999]] %>%
                plotly::add_lines(
                  x = rep(observed()$cutofftpdt, 2),
                  y = c(min(dfai$n), max(dfbi$upper)),
                  name = "prediction start",
                  line = list(dash="dash", color="grey"),
                  showlegend = FALSE)
            }


            if (i == 9999) {
              g[[1]] <- g[[1]] %>%
                plotly::layout(
                  annotations = list(
                    x = observed()$cutoffdt, y = 0, text = 'cutoff',
                    xanchor = "left", yanchor = "bottom",
                    font = list(size = 12), showarrow = FALSE))

              if (observed()$tp < observed()$t0) {
                g[[1]] <- g[[1]] %>%
                  plotly::layout(
                    annotations = list(
                      x = observed()$cutofftpdt, y = 0,
                      text = 'prediction start',
                      xanchor = "left", yanchor = "bottom",
                      font = list(size=12), showarrow = FALSE))
              }

              if (showEvent()) {
                g[[1]] <- g[[1]] %>%
                  plotly::add_lines(
                    x = range(dfsi$date), y = rep(target_d(), 2),
                    name = 'target events', showlegend = FALSE,
                    line = list(dash="dot",
                                color="rgba(128, 128, 128, 0.5")) %>%
                  plotly::layout(
                    annotations = list(
                      x = 0.95, xref = "paper", y = target_d(),
                      text = 'target events', xanchor = "right",
                      yanchor = "bottom", font = list(size = 12),
                      showarrow = FALSE))
              }
            }
          }

        } else {  # Design stage
          g <- list()

          for (i in c(9999, 1:k())) {
            dfsi <- dfs %>% filter(treatment == i)

            dfsi_enrollment <- dfsi %>% filter(parameter == "Enrollment")
            dfsi_event <- dfsi %>% filter(parameter == "Event")
            dfsi_ongoing <- dfsi %>% filter(parameter == "Ongoing")

            g[[(i+1) %% 9999]] <- plotly::plot_ly() %>%
              plotly::add_lines(
                data = dfsi_enrollment, x = ~t, y = ~n,
                line = list(width=2),
                name = "median prediction enrollment") %>%
              plotly::add_ribbons(
                data = dfsi_enrollment, x = ~t, ymin = ~lower, ymax = ~upper,
                fill = "tonexty", line = list(width=0),
                name = "prediction interval enrollment") %>%
              plotly::add_lines(
                data = dfsi_event, x = ~t, y = ~n,
                line = list(width=2),
                name = "median prediction event") %>%
              plotly::add_ribbons(
                data = dfsi_event, x = ~t, ymin = ~lower, ymax = ~upper,
                fill = "tonexty", line = list(width=0),
                name = "prediction interval event") %>%
              plotly::add_lines(
                data = dfsi_ongoing, x = ~t, y = ~n,
                line = list(width=2),
                name = "median prediction ongoing") %>%
              plotly::add_ribbons(
                data = dfsi_ongoing, x = ~t, ymin = ~lower, ymax = ~upper,
                fill = "tonexty", line = list(width=0),
                name = "prediction interval ongoing") %>%
              plotly::layout(
                xaxis = list(title = "Days since trial start",
                             zeroline = FALSE),
                yaxis = list(zeroline = FALSE)) %>%
              plotly::layout(
                annotations = list(
                  x = 0.5, y = 1,
                  text = paste0("<b>", dfsi$treatment_description[1], "</b>"),
                  xanchor = "center", yanchor = "bottom",
                  showarrow = FALSE, xref='paper', yref='paper'))


            if (i == 9999) {
              if (showEvent()) {
                g[[1]] <- g[[1]] %>%
                  plotly::add_lines(
                    x = range(dfsi$t), y = rep(target_d(), 2),
                    name = 'target events', showlegend = FALSE,
                    line = list(dash="dot",
                                color="rgba(128, 128, 128, 0.5")) %>%
                  plotly::layout(
                    annotations = list(
                      x = 0.95, xref = "paper", y = target_d(),
                      text = 'target events', xanchor = "right",
                      yanchor = "bottom", font = list(size = 12),
                      showarrow = FALSE))
              }
            }
          }
        }
      } else {
        g <- NULL
      }
    }

    g
  })


  mult_plot <- reactive({
    (to_predict() == "Enrollment only" &&
       (input$by_treatment || input$stage == 'Design stage') && k() > 1 &&
       "treatment" %in% names(pred()$enroll_pred$enroll_pred_df) &&
       length(table(pred()$enroll_pred$enroll_pred_df$treatment)) == k()+1) ||
      (to_predict() != "Enrollment only" &&
         (input$by_treatment || input$stage == 'Design stage') && k() > 1 &&
         "treatment" %in% names(pred()$event_pred$event_pred_df) &&
         length(table(pred()$event_pred$event_pred_df$treatment)) == k()+1)
  })


  observe({
    walk(1:6, function(i) {
      output[[paste0("pred_plot_output", i)]] <- renderPlotly({
        if (i <= k() + 1) {
          if (mult_plot()) {
            pred_plot()[[i]]
          } else {
            pred_plot()
          }
        } else {
          NULL
        }
      })
    })
  })


  pred_plot_outputs <- reactive({
    n = ifelse(mult_plot(), k() + 1, 1)
    outputs <- map(1:n, function(i) {
      plotlyOutput(paste0("pred_plot_output", i))
    })

    tagList(outputs)
  })


  output$pred_plot <- renderUI({
    pred_plot_outputs()
  })


  output$downloadEventSummaryData <- downloadHandler(
    filename = function() {
      paste0("event_summary_data_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      if (to_predict() == "Enrollment only") {
        eventsummarydata <- pred()$enroll_pred$enroll_pred_df
      } else {
        eventsummarydata <- pred()$event_pred$enroll_pred_df %>%
          bind_rows(pred()$event_pred$event_pred_df) %>%
          bind_rows(pred()$event_pred$ongoing_pred_df)
      }
      writexl::write_xlsx(eventsummarydata, file)
    }
  )


  output$downloadEventSubjectData <- downloadHandler(
    filename = function() {
      paste0("event_subject_data_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      eventsubjectdata <- pred()$subject_data
      writexl::write_xlsx(eventsubjectdata, file)
    }
  )


  # dosing panel outputs
  l <- reactive({
    if (input$stage != "Design stage" && !is.null(visitview())) {
      l = length(table(visitview()$drug))
      updateSelectInput(session, "l", selected=l)
    } else {
      l = as.numeric(input$l)
    }
    l
  })


  drug_description_df <- reactive({
    req(l())

    if (input$stage != "Design stage" && !is.null(visitview())) {
      a <- visitview() %>%
        group_by(drug, drug_name, dose_unit) %>%
        slice(n()) %>%
        select(drug, drug_name, dose_unit)
    } else {
      param = input[[paste0("drug_description_", l())]]
      a <- tibble(drug = 1:l(),
                  drug_name = as.character(param[,1]),
                  dose_unit = as.character(param[,2]))
    }

    a
  })


  drug_name <- reactive({
    drug_description_df()$drug_name
  })


  dose_unit <- reactive({
    drug_description_df()$dose_unit
  })


  treatment_by_drug <- reactive({
    req(k())
    req(l())
    if (input$stage != 'Design stage' && !is.null(df()) &&
        !is.null(visitview())) {

      vf <- dose_observed()$vf

      treatment_by_drug_df <- vf %>%
        group_by(treatment, drug, drug_name, dose_unit) %>%
        slice(n()) %>%
        select(treatment, drug, drug_name, dose_unit)

      treatment_by_drug <- matrix(
        0, nrow = k(), ncol = l(),
        dimnames = list(treatment_description(), drug_name()))
      for (i in 1:nrow(treatment_by_drug_df)) {
        treatment_by_drug[treatment_by_drug_df$treatment[i],
                          treatment_by_drug_df$drug[i]] = 1
      }

      updateMatrixInput(
        session, paste0("treatment_by_drug_", k(), "_", l()),
        value=treatment_by_drug)
    } else {
      param = input[[paste0("treatment_by_drug_", k(), "_", l())]]
      t = as.numeric(param)

      valid = all(t==1 | t==0)
      if (!valid) {
        showNotification(
          "Entries of treatment by drug matrix must be 1 or 0"
        )
      }

      req(valid)

      treatment_by_drug <- matrix(
        t, nrow = k(), ncol = l(),
        dimnames = list(treatment_description(), drug_name()))
    }

    treatment_by_drug
  })


  treatment_by_drug_df <- reactive({
    t = as.numeric(treatment_by_drug())

    treatment_by_drug_df <- tibble(
      treatment = rep(1:k(), l()),
      drug = rep(1:l(), each=k()),
      drug_name = rep(drug_name(), each=k()),
      dose_unit = rep(dose_unit(), each=k()),
      included = as.logical(t)) %>%
      filter(included) %>%
      select(treatment, drug, drug_name, dose_unit)

    treatment_by_drug_df
  })


  dosing_schedule_df <- reactive({
    req(l())

    x <- input[[paste0("dosing_schedule_", l())]]

    param <- tibble(
      drug = 1:l(),
      target_days = as.numeric(x[, "Days per Cycle"]),
      target_dose = as.numeric(x[, "Dose per Cycle"]),
      max_cycles = as.numeric(x[, "Number of Cycles"]))

    valid1 = all(param$target_days > 0 &
                   param$target_days == round(param$target_days))
    if (!valid1) {
      showNotification("Days per Cycle must be positive integers")
    }

    valid2 = all(param$target_dose >= 0)
    if (!valid2) {
      showNotification("Dose per Cycle must be nonnegative")
    }

    valid3 = all(param$max_cycles > 0 &
                   param$max_cycles == round(param$max_cycles))
    if (!valid3) {
      showNotification("Number of Cycles must be positive integers")
    }

    req(valid1 && valid2 && valid3)

    param
  })


  observeEvent(list(drug_name(), dose_unit(), l()), {
    updateMatrixInput(
      session, paste0("drug_description_", l()),
      value=matrix(c(drug_name(), dose_unit()),
                   nrow = l(), ncol = 2,
                   dimnames = list(NULL,
                                   c("Drug Name",
                                     "Dose Unit"))))
  })


  observeEvent(list(treatment_by_drug(), k(), l()), {
    updateMatrixInput(
      session, paste0("treatment_by_drug_", k(), "_", l()),
      value=matrix(treatment_by_drug(), nrow = k(), ncol = l(),
                   dimnames = dimnames(treatment_by_drug())))
  })


  observeEvent(list(drug_name(), l()), {
    updateMatrixInput(
      session, paste0("dosing_schedule_", l()),
      value=matrix(as.matrix(dosing_schedule_df()[,2:4]),
                   nrow = l(), ncol = 3,
                   dimnames = list(drug_name(),
                                   c("Days per Cycle",
                                     "Dose per Cycle",
                                     "Number of Cycles"))))
  })


  model_k0 <- reactive({
    if (input$stage != 'Design stage' && !is.null(df()) &&
        !is.null(visitview()) && !is.null(dosing_schedule_df())) {

      vf <- dose_observed()$vf %>%
        left_join(dosing_schedule_df(), by = "drug")

      df_k0 <- vf %>%
        filter(row_id == 1) %>%
        mutate(time = day,
               skipped = floor((time - target_days/2)/target_days) + 1)

      a = sum(df_k0$skipped)

      ifelse(a == 0, "Constant", input$model_k0)
    } else {
      input$model_k0
    }
  })


  # whether to show or hide k0, t0, t1, k0, ti, and di panels
  observeEvent(list(input$stage, model_k0(), pred_pp_only()), {
    if (input$stage == 'Design stage' || pred_pp_only()) {
      hideTab(inputId = "dosing_tabset", target = "k0_model_panel")
      hideTab(inputId = "dosing_tabset", target = "t0_model_panel")
      hideTab(inputId = "dosing_tabset", target = "t1_model_panel")
      hideTab(inputId = "dosing_tabset", target = "ki_model_panel")
      hideTab(inputId = "dosing_tabset", target = "ti_model_panel")
      hideTab(inputId = "dosing_tabset", target = "di_model_panel")
    } else {
      if (model_k0() != "Constant") {
        showTab(inputId = "dosing_tabset", target = "k0_model_panel")
        showTab(inputId = "dosing_tabset", target = "t1_model_panel")
      } else {
        hideTab(inputId = "dosing_tabset", target = "k0_model_panel")
        hideTab(inputId = "dosing_tabset", target = "t1_model_panel")
      }
      showTab(inputId = "dosing_tabset", target = "t0_model_panel")
      showTab(inputId = "dosing_tabset", target = "ki_model_panel")
      showTab(inputId = "dosing_tabset", target = "ti_model_panel")
      showTab(inputId = "dosing_tabset", target = "di_model_panel")
    }
  })


  model_t0 <- reactive({
    input$model_t0
  })


  model_t1 <- reactive({
    input$model_t1
  })


  model_ki <- reactive({
    input$model_ki
  })


  model_ti <- reactive({
    input$model_ti
  })


  model_di <- reactive({
    input$model_di
  })


  common_time_model <- reactive({
    if (!is.null(dosing_schedule_df())) {
      ifelse(length(unique(dosing_schedule_df()$target_days)) == 1,
             TRUE, FALSE)
    }
  })


  k0_fit <- reactive({
    if (!is.null(dose_observed()) & !is.null(dosing_schedule_df())) {
      vf = dose_observed()$vf %>%
        left_join(dosing_schedule_df(), by = "drug")

      if (common_time_model()) {
        # time from randomization to the first drug dispensing visit
        # use data from all drugs to inform the model
        df_k0 <- vf %>%
          filter(row_id == 1) %>%
          mutate(time = day,
                 skipped = floor((time - target_days/2)/target_days) + 1)

        k0_fit <- f_fit_ki(df_k0, model_k0(), nreps(), showplot = FALSE)
      } else {
        k0_fit <- list()
        for (h in 1:l()) {
          # observed dosing data for the drug under consideration
          vf1 <- vf %>% filter(drug == h)

          # time from randomization to the first drug dispensing visit
          df_k0 <- vf1 %>%
            filter(row_id == 1) %>%
            mutate(time = day,
                   skipped = floor((time - target_days/2)/target_days) + 1)

          k0_fit[[h]] <- f_fit_ki(df_k0, model_k0(), nreps(), showplot = FALSE)

          k0_fit[[h]]$fit_plot <- k0_fit[[h]]$fit_plot %>%
            plotly::layout(
              annotations = list(
                x = 0.5, y = 1,
                text = paste0("<b>", drug_name()[h], "</b>"),
                xanchor = "center", yanchor = "middle",
                showarrow = FALSE, xref='paper', yref='paper'))
        }
      }

      k0_fit
    }
  })


  t0_fit <- reactive({
    if (!is.null(dose_observed()) & !is.null(dosing_schedule_df())) {
      vf = dose_observed()$vf %>%
        left_join(dosing_schedule_df(), by = "drug")

      if (common_time_model()) {
        # time from randomization to the first drug dispensing visit
        # use data from all drugs to inform the model
        df_k0 <- vf %>%
          filter(row_id == 1) %>%
          mutate(time = day,
                 skipped = floor((time - target_days/2)/target_days) + 1)

        # no skipping
        df_t0 <- df_k0 %>%
          filter(skipped == 0) %>%
          mutate(left = time - 1, right = time)

        t0_fit <- f_fit_t0(df_t0, model_t0(), nreps(), showplot = FALSE)
      } else {
        t0_fit <- list()

        for (h in 1:l()) {
          # observed dosing data for the drug under consideration
          vf1 <- vf %>% filter(drug == h)

          # time from randomization to the first drug dispensing visit
          df_k0 <- vf1 %>%
            filter(row_id == 1) %>%
            mutate(time = day,
                   skipped = floor((time - target_days/2)/target_days) + 1)

          # no skipping
          df_t0 <- df_k0 %>%
            filter(skipped == 0) %>%
            mutate(left = time - 1, right = time)

          t0_fit[[h]] <- f_fit_t0(df_t0, model_t0(), nreps(), showplot = FALSE)

          t0_fit[[h]]$fit_plot <- t0_fit[[h]]$fit_plot %>%
            plotly::layout(
              annotations = list(
                x = 0.5, y = 1,
                text = paste0("<b>", drug_name()[h], "</b>"),
                xanchor = "center", yanchor = "middle",
                showarrow = FALSE, xref='paper', yref='paper'))
        }
      }

      t0_fit
    }
  })


  t1_fit <- reactive({
    if (!is.null(dose_observed()) & !is.null(dosing_schedule_df())) {
      vf = dose_observed()$vf %>%
        left_join(dosing_schedule_df(), by = "drug")

      if (common_time_model()) {
        # time from randomization to the first drug dispensing visit
        # use data from all drugs to inform the model
        df_k0 <- vf %>%
          filter(row_id == 1) %>%
          mutate(time = day,
                 skipped = floor((time - target_days/2)/target_days) + 1)

        # skipping
        df_t1 <- df_k0 %>%
          filter(skipped > 0) %>%
          mutate(k1 = skipped)

        if (nrow(df_t1) == 0) {
          t1_fit = NULL
        } else {
          t1_fit <- f_fit_ti(df_t1, model_t1(), nreps(), showplot = FALSE)
        }
      } else {
        t1_fit <- list()

        for (h in 1:l()) {
          # observed dosing data for the drug under consideration
          vf1 <- vf %>% filter(drug == h)

          # time from randomization to the first drug dispensing visit
          df_k0 <- vf1 %>%
            filter(row_id == 1) %>%
            mutate(time = day,
                   skipped = floor((time - target_days/2)/target_days) + 1)

          # skipping
          df_t1 <- df_k0 %>%
            filter(skipped > 0) %>%
            mutate(k1 = skipped)

          if (nrow(df_t1) == 0) {
            t1_fit[[h]] <- NULL
          } else {
            t1_fit[[h]] <- f_fit_ti(df_t1, model_t1(), nreps(),
                                    showplot = FALSE)

            t1_fit[[h]]$fit_plot <- t1_fit[[h]]$fit_plot %>%
              plotly::layout(
                annotations = list(
                  x = 0.5, y = 1,
                  text = paste0("<b>", drug_name()[h], "</b>"),
                  xanchor = "center", yanchor = "middle",
                  showarrow = FALSE, xref='paper', yref='paper'))
          }
        }
      }

      t1_fit
    }
  })


  ki_fit <- reactive({
    if (!is.null(dose_observed()) & !is.null(dosing_schedule_df())) {
      vf = dose_observed()$vf %>%
        left_join(dosing_schedule_df(), by = "drug")

      if (common_time_model()) {
        # gap time and number of skipped visits between drug dispensing visits
        df_ti <- vf %>%
          mutate(time = lead(day) - day,
                 skipped = pmax(floor((time - target_days/2)/target_days), 0),
                 k1 = skipped + 1) %>%
          filter(row_id < n())

        ki_fit <- f_fit_ki(df_ti, model_ki(), nreps(), showplot = FALSE)
      } else {
        ki_fit <- list()

        for (h in 1:l()) {
          # observed dosing data for the drug under consideration
          vf1 <- vf %>% filter(drug == h)

          # gap time & number of skipped visits between drug dispensing visits
          df_ti <- vf1 %>%
            mutate(time = lead(day) - day,
                   skipped = pmax(floor((time - target_days/2)/target_days),
                                  0),
                   k1 = skipped + 1) %>%
            filter(row_id < n())

          ki_fit[[h]] <- f_fit_ki(df_ti, model_ki(), nreps(), showplot = FALSE)

          ki_fit[[h]]$fit_plot <- ki_fit[[h]]$fit_plot %>%
            plotly::layout(
              annotations = list(
                x = 0.5, y = 1,
                text = paste0("<b>", drug_name()[h], "</b>"),
                xanchor = "center", yanchor = "middle",
                showarrow = FALSE, xref='paper', yref='paper'))
        }
      }

      ki_fit
    }
  })


  ti_fit <- reactive({
    if (!is.null(dose_observed()) & !is.null(dosing_schedule_df())) {
      vf = dose_observed()$vf %>%
        left_join(dosing_schedule_df(), by = "drug")

      if (common_time_model()) {
        # gap time and number of skipped visits between drug dispensing visits
        df_ti <- vf %>%
          mutate(time = lead(day) - day,
                 skipped = pmax(floor((time - target_days/2)/target_days), 0),
                 k1 = skipped + 1) %>%
          filter(row_id < n())

        ti_fit <- f_fit_ti(df_ti, model_ti(), nreps(), showplot = FALSE)
      } else {
        ti_fit <- list()

        for (h in 1:l()) {
          # observed dosing data for the drug under consideration
          vf1 <- vf %>% filter(drug == h)

          # gap time & number of skipped visits between drug dispensing visits
          df_ti <- vf1 %>%
            mutate(time = lead(day) - day,
                   skipped = pmax(floor((time - target_days/2)/target_days),
                                  0),
                   k1 = skipped + 1) %>%
            filter(row_id < n())

          ti_fit[[h]] <- f_fit_ti(df_ti, model_ti(), nreps(), showplot = FALSE)

          ti_fit[[h]]$fit_plot <- ti_fit[[h]]$fit_plot %>%
            plotly::layout(
              annotations = list(
                x = 0.5, y = 1,
                text = paste0("<b>", drug_name()[h], "</b>"),
                xanchor = "center", yanchor = "middle",
                showarrow = FALSE, xref='paper', yref='paper'))
        }
      }

      ti_fit
    }
  })


  di_fit <- reactive({
    if (!is.null(dose_observed()) & !is.null(dosing_schedule_df())) {
      vf = dose_observed()$vf

      di_fit <- list()
      for (h in 1:l()) {
        # observed dosing data for the drug under consideration
        vf1 <- vf %>% filter(drug == h)

        di_fit[[h]] <- f_fit_di(vf1, model_di(), nreps(), showplot = FALSE)

        di_fit[[h]]$fit_plot <- di_fit[[h]]$fit_plot %>%
          plotly::layout(
            annotations = list(
              x = 0.5, y = 1,
              text = paste0("<b>", drug_name()[h], "</b>"),
              xanchor = "center", yanchor = "middle",
              showarrow = FALSE, xref='paper', yref='paper'))
      }

      di_fit
    }
  })


  # loop over k0_fit, t0_fit, t1_fit, ki_fit, ti_fit plots
  fit_expr <- purrr::map(1:5, function(iter) {
    reactive({
      param = switch(iter, "k0", "t0", "t1", "ki", "ti")

      if (!pred_pp_only()) {
        fit = switch(iter, k0_fit(), t0_fit(), t1_fit(), ki_fit(), ti_fit())

        if (common_time_model()) {
          g <- fit$fit_plot
        } else {
          g <- purrr::map(1:l(), function(h) fit[[h]]$fit_plot)
        }
      } else {
        if (common_time_model()) {
          g <- NULL
        } else {
          g <- purrr::map(1:l(), function(h) NULL)
        }
      }

      list(param = param, fit_plot = g)
    })
  })


  fit_output_expr <- purrr::map(fit_expr, function(expr) {
    observe({
      walk(1:12, function(i) {
        output[[paste0(expr()$param, "_fit_output", i)]] <- renderPlotly({
          if (i <= l()) {
            if (!common_time_model()) {
              expr()$fit_plot[[i]]
            } else {
              expr()$fit_plot
            }
          } else {
            NULL
          }
        })
      })
    })

    reactive({
      n = ifelse(!common_time_model(), l(), 1)
      outputs <- purrr::map(1:n, function(i) {
        plotlyOutput(paste0(expr()$param, "_fit_output", i))
      })

      list(param = expr()$param, fit_outputs = tagList(outputs))
    })
  })


  observe({
    walk(fit_output_expr, function(expr) {
      output[[paste0(expr()$param, "_fit")]] <- renderUI({
        expr()$fit_outputs
      })
    })
  })


  # di fit information criteria
  output$di_fit_ic <- renderText({
    if (l() > 1 && !is.null(di_fit())) {
      aic = sum(sapply(di_fit(), function(fit) fit$fit$aic))
      bic = sum(sapply(di_fit(), function(fit) fit$fit$bic))
      aictext = paste("Total AIC:", formatC(aic, format = "f", digits = 2))
      bictext = paste("Total BIC:", formatC(bic, format = "f", digits = 2))
      text1 = paste0("<i>", aictext, ", ", bictext, "</i>")
    } else {
      text1 = NULL
    }

    if (!is.null(text1)) text1
  })


  # di_fit plot
  observe({
    walk(1:12, function(i) {
      output[[paste0("di_fit_output", i)]] <- renderPlotly({
        if (i <= l() && !is.null(di_fit())) {
          if (!pred_pp_only()) {
            di_fit()[[i]]$fit_plot
          } else {
            NULL
          }
        } else {
          NULL
        }
      })
    })
  })


  di_fit_outputs <- reactive({
    outputs <- map(1:l(), function(i) {
      plotlyOutput(paste0("di_fit_output", i))
    })

    tagList(outputs)
  })


  output$di_fit <- renderUI({
    di_fit_outputs()
  })


  # dosing prediction panel outputs
  dosing <- eventReactive(input$predict, {
    if (predict_dosing()) {
      req(pred()$event_pred)
      req(pred()$stage == input$stage && pred()$to_predict == to_predict())

      if (input$stage != 'Design stage') {
        shiny::validate(
          need(!is.null(df()),
               paste("Please upload enrollment and event data for",
                     "real-time prediction.")))

        shiny::validate(
          need(!is.null(visitview()),
               "Please upload dosing data for real-time prediction."))

        if (to_predict() == "Enrollment and event")
          shiny::validate(
            need(target_n() > observed()$n0,
                 "Target enrollment has been reached."))

        shiny::validate(
          need(target_d() > observed()$d0,
               "Target number of events has been reached."))
      }

      newEvents <- pred()$event_pred$newEvents
      if (!("treatment" %in% colnames(newEvents))) {
        newEvents$treatment = 1
        newEvents$treatment_description = "Treatment 1"
      }


      if (input$stage == 'Design stage') {
        drug_demand <- f_drug_demand(
          df = NULL,
          newEvents = newEvents,
          visitview = NULL,
          drug_description_df = drug_description_df(),
          treatment_by_drug = treatment_by_drug(),
          dosing_schedule_df = dosing_schedule_df(),
          pilevel = pilevel(),
          nyears = nyears(),
          pred_pp_only = pred_pp_only(),
          showplot = FALSE)
      } else {
        if (pred_pp_only()) {
          drug_demand <- f_drug_demand(
            df = df(),
            newEvents = newEvents,
            visitview = visitview(),
            dosing_schedule_df = dosing_schedule_df(),
            pilevel = pilevel(),
            nyears = nyears(),
            pred_pp_only = pred_pp_only(),
            showplot = FALSE)
        } else {
          drug_demand <- f_drug_demand(
            df = df(),
            newEvents = newEvents,
            visitview = visitview(),
            dosing_schedule_df = dosing_schedule_df(),
            model_k0 = model_k0(),
            model_t0 = model_t0(),
            model_t1 = model_t1(),
            model_ki = model_ki(),
            model_ti = model_ti(),
            model_di = model_di(),
            pilevel = pilevel(),
            nyears = nyears(),
            pred_pp_only = pred_pp_only(),
            showplot = FALSE)
        }
      }

      drug_demand
    }
  })


  dosing_plot <- reactive({
    if (!is.null(dosing())) {
      req(pred()$stage == input$stage && pred()$to_predict == to_predict())

      if (input$stage == 'Design stage') {
        g <- purrr::map(1:l(), function(h) {
          dosing()$dosing_pred_plot[[h]] %>%
            plotly::layout(
              annotations = list(
                x = 0.5, y = 1,
                text = paste0("<b>", drug_name()[h], "</b>"),
                xanchor = "center", yanchor = "bottom",
                showarrow = FALSE, xref='paper', yref='paper'))
        })
      } else {
        dfs <- dosing()$dosing_pred_df %>%
          filter(parameter == "observed data")

        if (showModelBased()) {
          dfs <- dfs %>% bind_rows(
            dosing()$dosing_pred_df %>%
              filter(parameter == "model based prediction"))
        }

        if (pred_pp_only() || showProtocolBased()) {
          dfs <- dfs %>% bind_rows(
            dosing()$dosing_pred_df %>%
              filter(parameter == "protocol based prediction"))
        }

        cutoffdt = dosing()$cutoffdt

        g <- purrr::map(1:l(), function(h) {
          dfa <- dfs %>% filter(drug == h & parameter == "observed data")
          dfb <- dfs %>% filter(drug == h &
                                  parameter == "model based prediction")
          dfb_pp <- dfs %>% filter(drug == h &
                                     parameter == "protocol based prediction")

          fig <- plotly::plot_ly() %>%
            plotly::add_lines(
              data = dfa, x = ~date, y = ~n,
              line = list(shape ="hv", width = 2),
              name = "observed") %>%
            plotly::add_lines(
              data = dfb, x = ~date, y = ~n, line = list(width = 2),
              name = "median prediction model") %>%
            plotly::add_ribbons(
              data = dfb, x = ~date, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width = 0),
              name = "prediction interval model") %>%
            plotly::add_lines(
              data = dfb_pp, x = ~date, y = ~n, line = list(width = 2),
              name = "median prediction protocol") %>%
            plotly::add_ribbons(
              data = dfb_pp, x = ~date, ymin = ~lower, ymax = ~upper,
              fill = "tonexty", line = list(width = 0),
              name = "prediction interval protocol") %>%
            plotly::layout(
              xaxis = list(title = "", zeroline = FALSE),
              yaxis = list(title = paste0("Doses to dispense ",
                                          "(", dfa$dose_unit[1], ")"),
                           zeroline = FALSE),
              annotations = list(
                x = 0.5, y = 1,
                text = paste0("<b>", dfa$drug_name[1], "</b>"),
                xanchor = "center", yanchor = "bottom",
                showarrow = FALSE, xref = 'paper', yref = 'paper'))

          if (nrow(dfb) > 0 || nrow(dfb_pp) > 0) {
            fig <- fig %>% plotly::add_lines(
              x = rep(cutoffdt, 2),
              y = c(min(dfa$n), max(dfb$upper, dfb_pp$upper, na.rm = TRUE)),
              line = list(dash = "dash"), showlegend = FALSE,
              name = "cutoff")
          }

          if (h==1) {
            fig <- fig %>%
              plotly::layout(
                annotations = list(
                  x = cutoffdt, y = 0, text = 'cutoff',
                  xanchor = "left", yanchor = "bottom",
                  font = list(size = 12), showarrow = FALSE))
          }

          fig
        })
      }
    } else {
      g <- NULL
    }

    g
  })


  observe({
    walk(1:12, function(i) {
      output[[paste0("dosing_plot_output", i)]] <- renderPlotly({
        if (i <= l()) {
          dosing_plot()[[i]]
        } else {
          NULL
        }
      })
    })
  })


  dosing_plot_outputs <- reactive({
    outputs <- map(1:l(), function(i) {
      plotlyOutput(paste0("dosing_plot_output", i))
    })

    tagList(outputs)
  })


  output$dosing_plot <- renderUI({
    dosing_plot_outputs()
  })


  output$downloadDosingSummaryData <- downloadHandler(
    filename = function() {
      paste0("dosing_summary_data_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      dosingsummarydata <- dosing()$dosing_pred_df
      writexl::write_xlsx(dosingsummarydata, file)
    }
  )


  output$downloadDosingSubjectData <- downloadHandler(
    filename = function() {
      paste0("dosing_subject_data_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      dosingsubjectdata <- dosing()$dosing_subject
      writexl::write_xlsx(dosingsubjectdata, file)
    }
  )


  # save inputs
  output$saveInputs <- downloadHandler(
    filename = function() {
      paste0("inputs_", Sys.Date(), "_drug_demand.rds")
    },

    content = function(file) {
      x <- list(
        stage = input$stage,
        to_predict = input$to_predict,
        to_predict2 = input$to_predict2,
        target_n = target_n(),
        target_d = input$target_d,
        pilevel = pilevel(),
        nyears = nyears(),
        to_show = input$to_show,
        to_show_dosing = input$to_show_dosing,
        by_treatment = input$by_treatment,
        predict_dosing = predict_dosing(),
        pred_pp_only = pred_pp_only(),
        k = k(),
        treatment_allocation = matrix(
          treatment_allocation(), ncol=1,
          dimnames = list(treatment_description(), "Size")),
        nreps = nreps(),
        seed = input$seed,

        enroll_prior = input$enroll_prior,
        poisson_rate = poisson_rate(),
        mu = mu(),
        delta = delta(),
        piecewise_poisson_rate = piecewise_poisson_rate(),
        enroll_model = input$enroll_model,
        nknots = nknots(),
        lags = lags(),
        accrualTime = matrix(
          accrualTime(), ncol = 1,
          dimnames = list(paste("Interval", 1:length(accrualTime())),
                          "Starting time")),

        event_prior = input$event_prior,
        exponential_survival = matrix(
          exponential_survival(), nrow = 1,
          dimnames = list(NULL, treatment_description())),
        weibull_survival = matrix(
          weibull_survival(), nrow = 2,
          dimnames = list(c("Shape", "Scale"), treatment_description())),
        llogis_survival = matrix(
          llogis_survival(), nrow = 2,
          dimnames = list(c("Location on log scale", "Scale on log scale"),
                          treatment_description())),
        lnorm_survival = matrix(
          lnorm_survival(), nrow = 2,
          dimnames = list(c("Mean on log scale", "SD on log scale"),
                          treatment_description())),
        piecewise_exponential_survival = matrix(
          piecewise_exponential_survival(), ncol = k()+1,
          dimnames = list(paste("Interval",
                                1:nrow(piecewise_exponential_survival())),
                          c("Starting time", treatment_description()))),
        event_model = input$event_model,
        piecewiseSurvivalTime = matrix(
          piecewiseSurvivalTime(), ncol = 1,
          dimnames = list(paste("Interval",
                                1:length(piecewiseSurvivalTime())),
                          "Starting time")),
        spline_k = spline_k(),
        spline_scale = input$spline_scale,

        l = l(),
        drug_description = matrix(
          as.character(input[[paste0("drug_description_", l())]]), ncol = 2,
          dimnames = list(NULL, c("Drug Name", "Dose Unit"))),
        treatment_by_drug = treatment_by_drug(),
        dosing_schedule = matrix(
          input[[paste0("dosing_schedule_", l())]], ncol = 3,
          dimnames = list(drug_name(),
                          c("Days per Cycle", "Dose per Cycle",
                            "Number of Cycles"))),
        model_k0 = model_k0(),
        model_t0 = model_t0(),
        model_t1 = model_t1(),
        model_ki = model_ki(),
        model_ti = model_ti(),
        model_di = model_di()
      )

      save(x, file = file)
    }
  )


  # load inputs
  observeEvent(input$loadInputs, {
    file <- input$loadInputs
    ext <- tools::file_ext(file$datapath)

    req(file)

    valid <- (ext == "rds")
    if (!valid) showNotification("Please upload an rds file")
    req(valid)

    load(file=file$datapath)

    updateRadioButtons(session, "stage", selected=x$stage)

    if (x$stage == 'Design stage' ||
        x$stage == 'Real-time before enrollment completion') {
      updateRadioButtons(session, "to_predict", selected=x$to_predict)
      updateNumericInput(session, "target_n", value=x$target_n)
    } else {
      updateRadioButtons(session, "to_predict2", selected=x$to_predict2)
    }

    if (x$to_predict == 'Enrollment and event' ||
        x$stage == 'Real-time after enrollment completion') {
      updateNumericInput(session, "target_d", value=x$target_d)
      updateCheckboxGroupInput(session, "to_show", selected=x$to_show)
    }

    if (x$stage != 'Design stage' &&
        !(x$stage == 'Real-time before enrollment completion' &&
          x$to_predict == 'Enrollment only') &&
        x$predict_dosing) {
      updateCheckboxGroupInput(session, "to_show_dosing",
                               selected=x$to_show_dosing)
    }

    updateNumericInput(session, "pilevel", value=x$pilevel)
    updateNumericInput(session, "nyears", value=x$nyears)

    updateCheckboxInput(session, "by_treatment", value=x$by_treatment)

    if (x$by_treatment &&
        (x$to_predict == 'Enrollment and event' ||
         x$stage == 'Real-time after enrollment completion')) {
      updateCheckboxInput(session, "predict_dosing", value=x$predict_dosing)

      if (x$predict_dosing) {
        updateCheckboxInput(session, "pred_pp_only", value=x$pred_pp_only)
      }
    }

    if (x$stage == 'Design stage' || x$by_treatment) {
      updateSelectInput(session, "k", selected=x$k)
    }

    if ((x$stage == 'Design stage' ||
         (x$by_treatment &&
          x$stage != 'Real-time after enrollment completion')) && x$k > 1) {
      updateMatrixInput(
        session, paste0("treatment_allocation_", x$k),
        value=x$treatment_allocation)
    }

    updateNumericInput(session, "nreps", value=x$nreps)
    updateNumericInput(session, "seed", value=x$seed)


    if (x$stage == 'Design stage') {
      updateRadioButtons(session, "enroll_prior", selected=x$enroll_prior)

      if (x$enroll_prior == "Poisson") {
        updateNumericInput(session, "poisson_rate", value=x$poisson_rate)
      } else if (x$enroll_prior == "Time-decay") {
        updateNumericInput(session, "mu", value=x$mu)
        updateNumericInput(session, "delta", value=x$delta)
      } else if (x$enroll_prior == "Piecewise Poisson") {
        updateMatrixInput(
          session, "piecewise_poisson_rate", value=x$piecewise_poisson_rate)
      }
    } else {
      if (x$stage == 'Real-time before enrollment completion') {
        updateRadioButtons(session, "enroll_model", selected=x$enroll_model)

        if (x$enroll_model == "B-spline") {
          updateNumericInput(session, "nknots", value=x$nknots)
          updateNumericInput(session, "lags", value=x$lags)
        } else if (x$enroll_model == "Piecewise Poisson") {
          updateMatrixInput(
            session, "accrualTime", value=x$accrualTime)
        }
      }
    }


    if (x$stage == 'Design stage') {
      if (x$to_predict == 'Enrollment and event') {
        updateRadioButtons(session, "event_prior", selected=x$event_prior)
      }

      if (x$event_prior == 'Exponential') {
        updateMatrixInput(
          session, paste0("exponential_survival_", x$k),
          value=x$exponential_survival)
      }

      if (x$event_prior == 'Weibull') {
        updateMatrixInput(
          session, paste0("weibull_survival_", x$k),
          value=x$weibull_survival)
      }

      if (x$event_prior == 'Log-logistic') {
        updateMatrixInput(
          session, paste0("llogis_survival_", x$k),
          value=x$llogis_survival)
      }

      if (x$event_prior == 'Log-normal') {
        updateMatrixInput(
          session, paste0("lnorm_survival_", x$k),
          value=x$lnorm_survival)
      }

      if (x$event_prior == 'Piecewise exponential') {
        updateMatrixInput(
          session, paste0("piecewise_exponential_survival_", x$k),
          value=x$piecewise_exponential_survival)
      }
    } else {
      if ((x$stage == 'Real-time before enrollment completion' &&
           x$to_predict == 'Enrollment and event') ||
          x$stage == 'Real-time after enrollment completion') {

        updateRadioButtons(session, "event_model", selected=x$event_model)

        if (x$event_model == "Piecewise exponential") {
          updateMatrixInput(
            session, "piecewiseSurvivalTime", value=x$piecewiseSurvivalTime)
        } else if (x$event_model == "Spline") {
          updateNumericInput(session, "spline_k", value=x$spline_k)
          updateRadioButtons(session, "spline_scale", selected=x$spline_scale)
        }
      }
    }


    if (x$predict_dosing) {
      updateSelectInput(session, "l", selected=x$l)

      updateMatrixInput(
        session, paste0("drug_description_", x$l),
        value=x$drug_description)

      updateMatrixInput(
        session, paste0("treatment_by_drug_", x$k, "_", x$l),
        value=x$treatment_by_drug)

      updateMatrixInput(
        session, paste0("dosing_schedule_", x$l),
        value=x$dosing_schedule)

      if (x$model_k0 != "Constant") {
        updateRadioButtons(session, "model_k0", selected=x$model_k0)
        updateRadioButtons(session, "model_t1", selected=x$model_t1)
      }

      updateRadioButtons(session, "model_t0", selected=x$model_t0)
      updateRadioButtons(session, "model_ki", selected=x$model_ki)
      updateRadioButtons(session, "model_ti", selected=x$model_ti)
      updateRadioButtons(session, "model_di", selected=x$model_di)
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
