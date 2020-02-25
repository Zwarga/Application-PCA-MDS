library(shiny)
library(stringr)
library(ggplot2)
library(MASS)
library(FactoMineR)
library(factoextra)

# Define UI for data upload app ----
ui <- fluidPage(
  #theme = "bootstrap.css",
  # App title ----
  titlePanel("PCA MDS Generator"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      
      
      # Input: Select a file ----
      fileInput("file1", "Choose CSV File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Checkbox if file has header ----
      checkboxInput("header", "Header", TRUE),
      
      # Input: Select separator ----
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),
      
      # Input: Select quotes ----
      radioButtons("quote", "Quote",
                   choices = c(None = "",
                               "Double Quote" = '"',
                               "Single Quote" = "'"),
                   selected = '"'),
      
      # Horizontal line ----
      tags$hr(),
      
      
      # Remove first column ---
      radioButtons("removecol", "Remove First Column",
                   choices = c("Yes" = '1',
                               "No" = "0"),
                   selected = "0"),
      
      
      #Use Label or Not
      radioButtons("label", "Use Label",
                   choices = c("Yes" = '1',
                               "No" = "0"),
                   selected = "0"),
      
      #Input label groupe 1
      textInput("label1","Label Group 1",value="Defaut"),
      
      #Input label groupe 2
      textInput("label2","Label Group 2",value="Defaut "),
      
      #Input label groupe 3
      textInput("label3","Label Group 3",value="Defaut  "),
      
      #Input label groupe 4
      textInput("label4","Label Group 4",value="Defaut   ")
      
     
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      tabsetPanel(
        tabPanel("Data",
                 tableOutput("contents")
        ),
        tabPanel("Scree-Plot", plotOutput("scree")),
        tabPanel("PCA", plotOutput("pca"),
                 #Input Title of the PCA
                 textInput("title_pca","PCA Title",value=""),
                 #Use Ellipse or Not
                 radioButtons("ellipse", "Add Ellipse",
                              choices = c("Yes" = "TRUE",
                                          "No" = "FALSE"),
                              selected = "FALSE"),
                 #Save the PCA plot
                 downloadLink("downloadPCA", "Download")),
        tabPanel("MDS", plotOutput("mds"),
                 #Input Title of the MDS
                 textInput("title_mds","MDS Title",value=""),
                 #Save the MDS plot
                 downloadLink("downloadMDS", "Download"))
      )
    )
  )
    
      
   
  
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  dataInput <- reactive({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    # when reading semicolon separated files,
    # having a comma separator causes `read.csv` to error
    tryCatch(
      {
        df <- read.csv(input$file1$datapath,
                       header = input$header,
                       sep = input$sep,
                       quote = input$quote)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    if(input$removecol==1){df<-df[,-1]}
    return(df)
    
  })
 
  output$contents <- renderTable({head(dataInput())})
  
 
  
  output$scree<-renderPlot({
    Aging_pca<-prcomp(t(dataInput()),scale=FALSE)
    plot(fviz_eig(Aging_pca))
  })
  
  plotInput <- reactive({
    if(input$label==1){
      Aging_pca<-prcomp(t(dataInput()),scale=FALSE)
      grouping<-rep("Defaut",ncol(dataInput()))
      grouping[str_detect(colnames(dataInput()),input$label1)]=input$label1
      grouping[str_detect(colnames(dataInput()),input$label2)]=input$label2
      grouping[str_detect(colnames(dataInput()),input$label3)]=input$label3
      grouping[str_detect(colnames(dataInput()),input$label4)]=input$label4
      p<-fviz_pca_ind(Aging_pca, geom=c("point", "text"), pointsize = 3, habillage=grouping, addEllipses =input$ellipse) 
      p + scale_color_brewer(palette="Set1") +ggtitle(input$title_pca)+
        theme(legend.text=element_text(size=10),axis.title=element_text(size=15),title=element_text(size=20))
    }else{
      Aging_pca<-prcomp(t(dataInput()),scale=FALSE)
      p<-fviz_pca_ind(Aging_pca, geom=c("point", "text"), pointsize = 3) 
      p + scale_color_brewer(palette="Set1") +ggtitle(input$title_pca)+
        theme(legend.text=element_text(size=10),axis.title=element_text(size=15),title=element_text(size=20))
    }
    
 
  })
  
  output$pca<-renderPlot({
    plot(plotInput())
  })
  
  output$downloadPCA <- downloadHandler(
    filename = "plotpca.png",
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
      ggsave(file, plot = plotInput(), device = device)
    }
    )
  
  mdsInput<-reactive({
    d <- dist(t(dataInput())) # euclidean distances between the rows
    fit <- isoMDS(d, k=2) # k is the number of dim
    
    # plot solution
    x <- fit$points[,1]
    y <- fit$points[,2]
    if(input$label==1){
      label1<-mean(x[str_detect(names(x),input$label1)])
      label2<-mean(x[str_detect(names(x),input$label2)])
      label3<-mean(x[str_detect(names(x),input$label3)])
      label4<-mean(x[str_detect(names(x),input$label4)])
      x<-c(label1,label2,label3,label4,x)
      names(x)[1:4]<-c(input$label1,input$label2,input$label3,input$label4)
      
      label1<-mean(y[str_detect(names(y),input$label1)])
      label2<-mean(y[str_detect(names(y),input$label2)])
      label3<-mean(y[str_detect(names(y),input$label3)])
      label4<-mean(y[str_detect(names(y),input$label4)])
      y<-c(label1,label2,label3,label4,y)
      names(y)[1:4]<-c(input$label1,input$label2,input$label3,input$label4)
    }
    plot_Aging<-data.frame(x,y)
    
    
    #rownames(plot_Aging)<-names(y)
    colourgrp<-rep("black",nrow(plot_Aging))
    colourgrp[str_detect(rownames(plot_Aging),input$label1)]="red"
    colourgrp[str_detect(rownames(plot_Aging),input$label2)]="blue"
    colourgrp[str_detect(rownames(plot_Aging),input$label3)]="green"
    colourgrp[str_detect(rownames(plot_Aging),input$label4)]="yellow"

    gg <- ggplot(plot_Aging,aes(x=x, y=y,colour=colourgrp,label = rownames(plot_Aging))) + geom_point()+
      geom_text(aes(label=rownames(plot_Aging)),hjust=0, vjust=0)+ guides(color = FALSE, size = FALSE)+ggtitle(input$title_mds)+
      theme(
        plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text(color="black", size=10, face="bold"),
        axis.title.y = element_text(color="black", size=10, face="bold"))
    gg
  })
  
  output$mds<-renderPlot({
   plot(mdsInput())
  })
  
  output$downloadMDS <- downloadHandler(
    filename = "plotmds.png",
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
      ggsave(file, plot = mdsInput(), device = device)
    }
  )

}

shinyApp(ui = ui, server = server)

