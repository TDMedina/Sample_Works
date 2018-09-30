
library(shiny)

#Define UI for application with two tab pages.
ui <- navbarPage('Tyler\'s App',
  
  #Makes tab 1.               
  tabPanel('Random Histogram Simulator',
           
           #Page title.
           titlePanel(h3('Histogram from Random Normal Distribution')),
           
           #Splits the page into a side panel and a main area.
           sidebarLayout(
             
             #Creates the side panel, which has a numeric input and two sliders.
             sidebarPanel('',
                          numericInput('sample_size','Number of Samples',1000,0,step=1),
                          sliderInput('mean_slider','Mean',0,100,50,1,sep=','),
                          sliderInput('sd_slider','Standard Deviation',0,25,10,1,sep=',')),
             
             #Creates the main panel, which contains the histogram plot.
             mainPanel(h4('Generated Histogram:'),
                       plotOutput('histoplot'))
             
                        )
           ),
  
  #Makes tab 2.
  tabPanel('CSV to Boxplot Display',
           
           #Page title on the second tab.
           titlePanel(h3('CSV Displayer and Boxplot Generator')),
           
           #Splits the page into a side panel and a main area.
           sidebarLayout(
             
             #Creates the side panel on page two, which has a file input, column choice, and
             #color select radio button.
             sidebarPanel(fileInput('file_upload','CSV File',
                                    accept = c(
                                      "text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                          
                          selectInput('col_choice','Column to Display',
                                      choices=c('Age'='age', 'Height'='height',
                                                'Weight'='weight', 'BMI'='bmi')),
                          
                          radioButtons('color_button', 'Boxplot Color',
                                       c('Blue'='deepskyblue3','Green'='seagreen2',
                                         'Purple'='mediumorchid1','Red'='firebrick'))
                           ),
             
             #Creates the main panel, which has a table display and a boxplot.
             mainPanel(tableOutput('csv_display'),
                       plotOutput('col_boxplot'))
             
                        )
           
           )
)

#Define server logic.
server <- function(input, output) {
  
  #Tab 1: Logic to take the values from the inputs and apply them to the histogram.
  output$histoplot <- renderPlot({
    hist(rnorm(input$sample_size, input$mean_slider, input$sd_slider),
         main='Histogram from Random Normal Distribution',
         xlab='Bins',
         ylab='Number in Bin',
         border='turquoise4',
         col='turquoise3')})
  
  #Tab 2: Logic to display the contents of the uploaded file in a table.
  output$csv_display <- renderTable({
    req(input$file_upload)
    read.csv(input$file_upload$datapath)})
  
  #Tab 2: Logic to take the file data and make a boxplot from it.
  output$col_boxplot <- renderPlot({
    req(input$file_upload)
    the_file <- read.csv(input$file_upload$datapath)
    column_choice <- input$col_choice
    label_names <- c('age'='Age','weight'='Weight','height'='Height','bmi'='BMI')
    boxplot(the_file[column_choice], 
            main=paste('Boxplot of ',label_names[column_choice],sep=''),
            ylab='Value',
            col = input$color_button)})
  
}

#Run the application.
shinyApp(ui = ui, server = server)