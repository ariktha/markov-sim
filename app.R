library(shiny)
library(shinyMatrix)
library(tidyverse)
library(msm)
library(here)

set.seed(1)

get_summ <- function(data_tab){
  summ_dat <- data_tab %>%
    group_by(subject, state) %>%
    summarise(n = n()) %>% ungroup() %>% 
    group_by(subject) %>%
    do ({                        
      if (max(.$state) == 4) {    
        summ_op1 <- .
      } else {                
        new <- tibble(subject = .$subject[1],
                      n = 1,
                      state = 5)
        summ_op1 <- rbind(., new)
      }
      summ_op1
    }) %>% ungroup()
  
  dat_char <- summ_dat %>% group_by(subject) %>%
    summarise(max_state = max(state),
              all_states = list(as.character(state))) %>% 
    ungroup() %>% rowwise() %>%
    mutate(ever_col = ifelse("2" %in% all_states, TRUE, FALSE),
           ever_inf = ifelse("3" %in% all_states, TRUE, FALSE)) %>%
    dplyr::select(-all_states) %>% ungroup() %>%
    mutate(across(max_state, as.factor))
  
  dat_char %>% summarise(percent_discharged = nrow(subset(dat_char, max_state == 4))/nrow(dat_char),
                         percent_obs_colonized = sum(ever_col)/nrow(dat_char),
                         percent_obs_infected = sum(ever_inf)/nrow(dat_char))
  
}

ui <- fluidPage(
  titlePanel("Generic Markov Simulation Model"),
  sidebarLayout(
    sidebarPanel(
      numericInput(inputId = "subj_count", "Number of subjects", value = 50),
      numericInput(inputId = "bytime", "Observe every ___ days", value = 2),
      numericInput(inputId = "endtime", "Total observation time", value = 30),
      matrixInput(inputId = "qmat", "Transition Intensity Matrix", class = "numeric",
                  value = matrix(c(0,    0.1, 0,    0.02, 
                                   0.02, 0,   0.05, 0.04,
                                   0,    0.2, 0,    0.1,
                                   0,    0,   0,    0), 
                                  nrow = 4, ncol = 4, byrow = T, 
                                  dimnames = list(c("U", "C", "I", "D"),
                                                  c("U", "C", "I", "D")))),
      matrixInput(inputId = "init_ct", "Initial States", class = "numeric",
                  value = matrix(c(30, 10, 10),
                                 nrow = 3, ncol = 1, byrow = T,
                                 dimnames = list(c("U", "C", "I"),
                                                 c("Num of subjects")))),
      actionButton("goButton", "Go!", class = "btn-success"),
      verbatimTextOutput("Text_op1"),
      imageOutput("Image")),
    mainPanel(
      tabsetPanel(
        tabPanel("Simulated and observed data characteristics",
                 fluidRow(
                   column(width = 8,
                          tableOutput("char_simData")),
                   column(width = 8,
                          tableOutput("char_obsData"))
                 )),
        tabPanel("Markov model output", 
                 fluidRow(
                   column(width = 8,
                          plotOutput("Plot1", width = "600px", height = "600px")),
                   column(width = 8,
                          plotOutput("Plot2", width = "500px", height = "500px")),
                   column(width = 8,
                          tableOutput("sojournTime")),
                   column(width = 8,
                          tableOutput("qMat")),
                   column(width = 8,
                          tableOutput("pNext")))
                 ),
        tabPanel("Potential Additions",
                 verbatimTextOutput("AddnDir"))
        )))
)

server <- function(input, output) {
  
  output$Image <- renderImage({
    list(src = here("model_structure.png"),
         contentType = 'image/png',
         alt = "Model Structure", height="400px", width="400px", align="center")
  }, deleteFile = F)
  
  output$Text_op1 <- renderText("Assumptions:\n- Constant force of infection\n- History of antibiotics doesn't contribute to C/I risk")
  output$Text_op2 <- renderText("Sojourn Times") %>%
    bindEvent(input$goButton)
  
  simData <- reactive({
    subj_ids <- c(1:input$subj_count)
    data_cf <- cross_join(tibble(time = seq(1, input$endtime, by = 1)), tibble(subject = subj_ids)) %>%
      arrange(subject, time)
    init_states <- rep(1:3, c(input$init_ct[1,1], input$init_ct[2,1], input$init_ct[3, 1])) 
    op_cf_raw <- simmulti.msm(data_cf, input$qmat, start = init_states, death = TRUE)
    op_cf_raw %>% mutate(time = floor(time)) %>%
      group_by(subject) %>% mutate(max_time = max(time), is_max = ifelse(time == max_time, TRUE, FALSE)) %>%
      mutate(sum_is_max = sum(is_max)) %>% ungroup() %>% rowwise() %>%
      mutate(keep = ifelse(is_max == TRUE & sum_is_max > 1 & state != 4, FALSE, TRUE)) %>%
      dplyr::filter(keep) %>% dplyr::select(subject, time, state)
  }) %>%
  bindEvent(input$goButton)
  
  obsData <- reactive({
    obstimes <- seq(1, input$endtime, by = input$bytime)
    subj_ids <- c(1:input$subj_count)
    data1 <- cross_join(tibble(time = obstimes), tibble(subject = subj_ids)) %>%
      arrange(subject, time)
    op_cf_disc <- simData() %>% dplyr::filter(state == 4)
    data1 %>% left_join(simData(), by = join_by(time, subject)) %>% rbind(op_cf_disc) %>% distinct() %>% drop_na(state) %>%
      arrange(subject, time)
  }) %>%
  bindEvent(input$goButton)
  
  model_obsData <- reactive({
    msm(state ~ time, subject = subject, qmatrix = input$qmat, data = obsData(), deathexact = TRUE)
  }) %>%
  bindEvent(input$goButton)
  
  output$char_simData <- renderTable({
    data.frame(get_summ(simData()))
  }, caption = "Simulated Data",
  caption.placement = getOption("xtable.caption.placement", "top")) %>%
  bindEvent(input$goButton)
  
  output$char_obsData <- renderTable({
    data.frame(get_summ(obsData()))
  }, caption = "Observed Data",
  caption.placement = getOption("xtable.caption.placement", "top")) %>%
  bindEvent(input$goButton)
  
  output$Plot1 <- renderPlot(({
    plot.prevalence.msm(model_obsData())
    }), height = 600, width = 600)
  
  output$Plot2 <- renderPlot(({
    plot.msm(model_obsData())
    }), height = 500, width = 500)
  
  output$sojournTime <- renderTable(
    sojourn.msm(model_obsData()), rownames = TRUE, caption = "Sojourn Times",
    caption.placement = getOption("xtable.caption.placement", "top")
  )
  
  output$qMat <- renderTable(
    qmatrix.msm(model_obsData())$estimate, rownames = TRUE, caption = "Estimated Transition Intensity Matrix",
    caption.placement = getOption("xtable.caption.placement", "top")
  )
  
  output$pNext <- renderTable(
    pnext.msm(model_obsData())$estimate, rownames = TRUE, caption = "Probability that the next move of a process in (row) is to (column)",
    caption.placement = getOption("xtable.caption.placement", "top")
  )
  
  
  output$AddnDir <- renderText("To add:\n
          - Currently replacing the state on date of discharge with discharge (losing 1 \"observation\")\n
          - Currently not using simmulti function's in-built method to specify observation times\n
          - Look at colonization status at discharge
          - Add reasonable values as default q matrix")
}

shinyApp(ui = ui, server = server)
