# A shiny web app to explore the discharge record of the Gribble Gap watershed in Cullowhee, NC
# J.P. Gannon
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(EcoHydRology)
library(shiny)
library(shinythemes)
library(tidyverse)
library(lubridate)
library(patchwork)
library(readr)
library(grid)



# Define UI
ui <- fluidPage(
              #google analytics
              tags$head(includeHTML(("google-analytics.html"))),
                
                #set theme for app
                theme = shinytheme("lumen"),
                
                #title
                titlePanel("Gribble Gap Discharge (v2)"),
                
                sidebarLayout(
                  sidebarPanel(
                    h4("In this app you can explore Precipitation, Streamflow, and Stream Temperature for the Gribble Gap Watershed"),
                    h5(strong("Location"), "Cullowhee, NC, USA"),
                    h5(strong("Coordinates:"), "35.303667, -83.205814"),
                    h5(strong("Land Cover: "), "Forested"),
                    h5(strong("Size: "), "43 ha"),
                    h5(strong("Climate: "), "Temperate Oceanic"),
                    h5(strong("TO ZOOM: On the streamflow plot, click and drag and then double click. 
                      Double click again to zoom to full extent.")),
                    # Select type data to plot
                    selectInput(inputId = "type", label = strong("Streamflow plot units:"),
                                choices = c("Discharge (L/s)","Discharge (mm/day)"),
                                selected = "Discharge (L/s)"),
                    
                    # Select date range to be plotted
                    #dateRangeInput("date", strong("Date range"), start = "2015-08-01", end = "2019-02-12",
                     #              min = "2015-08-01", max = "2019-2-12"),
                    
                    # Select whether to log primary axis
                    checkboxInput(inputId = "log", label = strong("Log streamflow y axis."), value = FALSE),
                    
                    # Select to show calculated baseflow, and choose filter parameter
                    checkboxInput(inputId = "bf", label = strong("Show calculated baseflow."), value = FALSE),
                    numericInput(inputId = "filter", label = "Filter parameter for baseflow calculation", value = 0.925, min = 0, max = 1),
                   
                    #select to show water temperature
                    checkboxInput(inputId = "temp", label = strong("Show water temperature."), value = FALSE),
                    
                    #add a horizontal line to the plot to facilitate data comparisons
                    numericInput(inputId = "horizLine", label = strong("Add a horizontal line to discharge plot at:"), 
                                 value = 0, min = 0, max = 100),
                    
                    #Identiy range of precip to highlight
                    sliderInput("Prange", "Range of Precip Values to Highlight", value = c(0, 150), min = 0, max = 100),
                    
                    #change the size of the text in the plot using a slider
                    sliderInput(inputId = "mag", label = strong("Plot text size."), min = 11, max = 24, value = 12),
                    
                    #info
                    tags$div("Built by: JP Gannon", tags$br(), "Collegiate Assistant Professor", tags$br(), "Virginia Tech"),
                    a(href = "https://github.com/jpgannon/GG_Discharge_Shiny_App/blob/master/app.R", "Source Code"),
                    a(href = "https://github.com/jpgannon/GG_Discharge_Shiny_App/blob/master/app.R", "Source Code")
                    
                  ),
                  
                  # Output: Set up plot panel
                  mainPanel(
                    plotOutput(outputId = "precip", height = "150px"),
                    plotOutput(outputId = "lineplot", height = "250px",
                               dblclick = "plot1_dblclick",
                               brush = brushOpts(
                                 id = "plot1_brush",
                                 resetOnNew = TRUE)),
                    plotOutput(outputId = "therest")
                  )
                  
                ))


# Define server function
server <- function(input, output) {
  #load data when running locally
  setwd("/Volumes/GoogleDrive/My Drive/Gribble Gap Discharge/GribbleGap_Discharge_v2")
  
  ranges <- reactiveValues(x = as.POSIXct(c(start = "2015-08-01", end = "2019-03-01")))
  maxrange <- reactiveValues(x = as.POSIXct(c(start = "2015-08-01", end = "2019-03-01")))
  
  
  #Read hourly discharge data and set up date time column
  Q <- read_csv("GGQ_to_FEB19_hourly.csv") %>%
    mutate(time = mdy_hm(time))
  
  
  #Q_narm must not have NA's
  Q_narm <- Q
  
  #Read precip data. Currently from: 
  #https://www.ncdc.noaa.gov/cdo-web/datasets/GHCND/stations/GHCND:USC00312200/detail
  P <- read_csv("Cullowhee Precip 2015-2018.csv") %>%
    mutate(Date = mdy(Date)) %>%
    mutate(Highlight = NA)
 
  #make static ECDF of entire data record
  ecdf_all <- ecdf(Q$disch)
  Qrange_all <- seq(0,max(Q$disch, na.rm = TRUE),.1)
  Probs_all <- ecdf_all(Qrange_all)
  Qall_E <- tibble(Probs_all, Qrange_all)
  colnames(Qall_E) <- c("Prob","Q")
  
  #run baseflow separation and create ecdf for baseflow
  QBF <- reactive({
    req(input$filter)
    
    Q_narm$BFQ <- BaseflowSeparation(Q_narm$disch[is.na(Q_narm$disch)==FALSE], input$filter)$bt
    Q_narm$BFQmmd <- BaseflowSeparation(Q_narm$GGwsd[is.na(Q_narm$GGwsd)==FALSE], input$filter)$bt
    
    left_join(Q, select(Q_narm, time, BFQ, BFQmmd), by = "time")
  })
  
  bfecdf <- reactive({
    req(input$filter)

    Q_narm$BFQ <- BaseflowSeparation(Q_narm$disch[is.na(Q_narm$disch)==FALSE], input$filter)$bt
    Q <- left_join(Q, select(Q_narm, time, BFQ),by = "time")

    #make static ECDF of entire baseflow data record
    bfecdf_all <- ecdf(Q$BFQ)
    bfQrange_all <- seq(0,max(Q$BFQ),.1)
    bfProbs_all <- bfecdf_all(bfQrange_all)

    bfecdf <- tibble(bfProbs_all,bfQrange_all)
    colnames(bfecdf) <- c("Prob","Q")
    bfecdf
  })
  
  output$precip <- renderPlot({
    PrecipTrim <- P %>% filter(between(as.POSIXct(Date), ranges$x[1], ranges$x[2]))
    Ptotal <- round(sum(PrecipTrim$Precmm, na.rm = TRUE),2)
    
    # Create text
    grob <- grobTree(textGrob(paste("Total Precip:", Ptotal, "mm"), x=0.1,  y=0.1, hjust=0,
                              gp=gpar(col="black", fontsize=13, fontface="italic")))
    
    P <- P %>% mutate(Highlight = ifelse(Precmm >= input$Prange[1] & Precmm <= input$Prange[2], Precmm, NA))
    
    P %>% filter(between(as.POSIXct(Date), ranges$x[1], ranges$x[2])) %>%
      ggplot(aes(Date, Precmm))+
      geom_bar(stat = "identity")+
      scale_y_reverse()+
      theme_classic()+
      geom_bar(aes(Date, Highlight), fill = "blue", stat = "identity")+
      ylab("Precip (mm)")+
      xlab(element_blank())+
      theme(text = element_text(size=input$mag))+
      theme(plot.margin = unit(c(5.5,5.5,5.5,20), "pt"))+
      annotation_custom(grob)
    })
  
  # Create scatterplot object the plotOutput function is expecting
  output$lineplot <- renderPlot({

    Qsub <- QBF() %>% filter(between(time, ranges$x[1], ranges$x[2]))
    
    #calculate total discharge in mm for time period
    #ggwsd is in mmd
    perhour <- Qsub$GGwsd/24 #days to hours
    
    Qtotal <- sum(perhour, na.rm = TRUE)  
    Qtotal <- round(Qtotal, 2)
    
    # Create text
    grob <- grobTree(textGrob(paste("Total Discharge:", Qtotal, "mm"), x=0.1,  y=0.95, hjust=0,
                              gp=gpar(col="black", fontsize=13, fontface="italic")))
   
    
    
    if(input$type == "Discharge (L/s)") {
      plot1 <- Qsub %>%
        ggplot(aes(time, disch))+
        geom_line()+
        ylab(input$type)+
        xlab(element_blank())+
        theme_classic()+
        theme(text = element_text(size=input$mag))+
        geom_hline(yintercept = input$horizLine, color = "red")+
        theme(plot.margin = unit(c(5.5,5.5,5.5,15), "pt"))+
        annotation_custom(grob)
        
    }
    
    if(input$type == "Discharge (mm/day)") {
      plot1 <- Qsub %>% 
        ggplot(aes(time, GGwsd))+
        geom_line()+
        ylab(input$type)+
        xlab(element_blank())+
        theme_classic()+
        theme(text = element_text(size=input$mag))+
        geom_hline(yintercept = input$horizLine, color = "red")+
        theme(plot.margin = unit(c(5.5,5.5,5.5,15), "pt"))+
        annotation_custom(grob)
    }
    
  
    if(input$log) plot1 <- plot1 + scale_y_log10()  
    
    if(input$bf & input$type == "Discharge (L/s)")
      plot1 <- plot1 + geom_line(aes(time, BFQ), color = "blue")
    
    if(input$bf & input$type == "Discharge (mm/day)")
      plot1 <- plot1 + geom_line(aes(time, BFQmmd), color = "blue")
    
   plot1
     })
  
  output$therest <- renderPlot({
    Qsub <- QBF() %>% filter(between(time, ranges$x[1], ranges$x[2]))
    
    plottemp <- Qsub %>% 
      ggplot(aes(time, TEMPERATURE))+
      geom_line()+
      ylab("Water Temperature (C)")+
      xlab(element_blank())+
      theme_classic()+
      theme(text = element_text(size=input$mag))
    
    #discharge ecdf calcs for selected data
    ecdf_user <- ecdf(Qsub$disch)
    Qrange_sel <- seq(0,max(Qsub$disch, na.rm = TRUE),.1)
    probs_user <- ecdf_user(Qrange_sel)
    Q_Esel <- tibble(probs_user, Qrange_sel)
    colnames(Q_Esel) <- c("Prob","Q")

    Q_ecdfs <- bind_rows("Whole Record" = Qall_E,"Selected Record" =  Q_Esel, .id = "Period")

    #baseflow ecdf calcs for selected data
    bfecdf_user <- ecdf(Qsub$BFQ)
    bfQrange_sel <- seq(0,max(Qsub$BFQ),.1)
    bfprobs_user <- bfecdf_user(bfQrange_sel)
    bf_ecdf_sel <- tibble(bfprobs_user, bfQrange_sel)
    colnames(bf_ecdf_sel) <- c("Prob","Q")
    
    BFecdfs <- bind_rows("Selected Record" = bf_ecdf_sel,"Whole Record" = bfecdf(), .id = "Period")
    
    QECDFplot <- Q_ecdfs %>%
      ggplot(aes(Q, Prob, color = Period))+
      geom_line()+
      theme_classic()+
      ylab("Probability of Non-Exceedance")+
      xlab("Discharge (L/s)")+
      theme(legend.position = "bottom")+
      theme(text = element_text(size=input$mag))

    BF_ecfd_plot <- BFecdfs %>%
      ggplot(aes(Q, Prob, color = Period))+
      geom_line()+
      theme_classic()+
      ylab("Probability of Non-Exceedance")+
      xlab("Baseflow (L/s)")+
      theme(legend.position = "none")+
      theme(text = element_text(size=input$mag))


    #QECDFplot #+ BF_ecfd_plot
    if(input$temp){
      plottemp / (QECDFplot + BF_ecfd_plot)
    }else{
      plottemp <- Qsub %>%
        ggplot(aes(time, TEMPERATURE)) + geom_blank() + theme_void()
      (QECDFplot + BF_ecfd_plot) / plottemp
    }

  })
  
   observeEvent(input$plot1_dblclick,
               {
                 brush <- input$plot1_brush
                 if (!is.null(brush)) 
                 {
                   ranges$x <- c(brush$xmin, brush$xmax)
                   
                 } else {
                   ranges$x <- maxrange$x
                   }
               })
}

# Create Shiny object
shinyApp(ui = ui, server = server)

