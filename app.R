library(shiny)
library(shinydashboard)
library(leaflet)
library(data.table)
library(DT)
library(plyr)
library(ggplot2)
library(dplyr)
library(gplots)


krsa <-
  function (file ,iterations, numSigSubs, sigSubs) {
    subKinase <- read.table(file, header = TRUE, sep = "\t", as.is=T, quote = "\"", check.names=FALSE, stringsAsFactors=FALSE)
    kinaseSep=" "
    subKinaseOnly=toupper(subKinase[,2])

    splitKinase=NULL
    for(i in 1:length(subKinaseOnly)){
      splitKinase=paste0(splitKinase, subKinaseOnly[i], sep=kinaseSep )
    }
    uniqueKin=unique(unlist(strsplit(splitKinase, split=kinaseSep)))
    scores=data.frame(matrix(ncol = iterations, nrow = length(uniqueKin)), row.names=uniqueKin)
    for(i in 1:iterations){
      sampleList=sample(subKinaseOnly, numSigSubs) 
      scores[,i]=iterationScore(sampleList, kinaseSep, uniqueKin)
    }
    
    scores=cbind(scores, rowMeans(scores)) #get the average kinase amount
    yourSubs=subset(subKinase, subKinase$Substrates %in% sigSubs)
    yourScores=data.frame(matrix(ncol = 1, nrow = length(uniqueKin)), row.names=uniqueKin)
    yourScores[,1]=iterationScore(yourSubs[,2], kinaseSep, uniqueKin)
    
    
    #section to make overall table
    overallTable=as.data.frame(cbind(Kinase=uniqueKin, DxAverage=yourScores[,1], ReSamplingAverage=scores[,iterations+1]))
    overallTable[,2]=as.numeric(as.character(overallTable[,2]))
    overallTable[,3]=as.numeric(as.character(overallTable[,3]))
    overallTable=cbind(overallTable, AbsoluteDiference=abs(overallTable[,2]-overallTable[,3]))
    overallTable=cbind(overallTable, StandardDeviation=unlist(apply(scores[,1:iterations], 1, sd)))
    overallTable=cbind(overallTable, Minus2Sd=overallTable[,3]-2*(overallTable[,5]) , Plus2SD=overallTable[,3]+2*(overallTable[,5]))
    overallTable=cbind(overallTable, Down=(ifelse(overallTable[,2]<overallTable[,6], "YES", "x")), Up=(ifelse(overallTable[,2]>overallTable[,7], "YES", "x")))
    overallTable=cbind(overallTable, ZScore=abs((overallTable[,2]-overallTable[,3])/overallTable[,5]))
    overallTable2=overallTable[order(overallTable$ZScore, decreasing=TRUE),]
    
    #kittyYES<<-overallTable2
    print("Done")
    return(list(overallTable2, scores))
  }

#Function to get a score for each kinase for each iteration
iterationScore <- function(sampleList, kinaseSep, uniqueKin){
  
  splitKinase=NULL
  for(i in 1:length(sampleList)){
    splitKinase=paste0(splitKinase, sampleList[i], sep=kinaseSep)
  }
  splitKinase2=unlist(strsplit(splitKinase, split=kinaseSep))
  returnVector=c()
  for(i in 1:length(uniqueKin)){
    temp=paste("^",uniqueKin[i],"$", sep="")
    returnVector[i]=length(grep(temp, splitKinase2))
  }
  return(returnVector)
  
}

dropdownButton <- function(label = "", status = c("default", "primary", "success", "info", "warning", "danger"), ..., width = NULL) {
  
  status <- match.arg(status)
  # dropdown button content
  html_ul <- list(
    class = "dropdown-menu",
    style = if (!is.null(width)) 
      paste0("width: ", validateCssUnit(width), ";"),
    lapply(X = list(...), FUN = tags$li, style = "margin-left: 10px; margin-right: 10px;")
  )
  # dropdown button apparence
  html_button <- list(
    class = paste0("btn btn-", status," dropdown-toggle"),
    type = "button", 
    `data-toggle` = "dropdown"
  )
  html_button <- c(html_button, list(label))
  html_button <- c(html_button, list(tags$span(class = "caret")))
  # final result
  tags$div(
    class = "dropdown",
    do.call(tags$button, html_button),
    do.call(tags$ul, html_ul),
    tags$script(
      "$('.dropdown-menu').click(function(e) {
      e.stopPropagation();
});")
  )
  }


make_graphs <- function(controlSamplesTest, experimentalSamplesTest, fullDataSet, opt1, opt2, lower, upper){
  #Get the correct rows and columns and substrate names
  sampleNameCol=grep("Sample name", fullDataSet)
  sampleNameRow=which(fullDataSet[,sampleNameCol]=="Sample name")
  cycleCol=grep("Cycle", fullDataSet)
  cycleRow=which(fullDataSet[,cycleCol]=="Cycle")
  exposureTimeCol=grep("Exposure time", fullDataSet)
  exposureTimeRow=which(fullDataSet[,exposureTimeCol]=="Exposure time")
  sequenceCol=grep("Sequence", fullDataSet)
  sequenceRow=which(fullDataSet[,sequenceCol]=="Sequence")
  substrateCol=unlist(as.numeric(which(apply(fullDataSet, 2, grep, pattern="^ID$", value=FALSE, perl=TRUE)!=0)))
  substrateNames=fullDataSet[(sequenceRow+1):nrow(fullDataSet), substrateCol]
  print(substrateNames)
  substrateNames2=substrateNames[which(substrateNames!="#REF")]
  
  #Get list of unique control samples (if the user selected more than one sample type)
  controlSamples_u=unique(as.character(controlSamplesTest[sampleNameRow, ]))
  experimentalSamples_u=unique(as.character(experimentalSamplesTest[sampleNameRow, ]))
  allSamples_u=union(controlSamples_u, experimentalSamples_u)
  combinedDF=cbind(controlSamplesTest, experimentalSamplesTest)
  par(mfrow=c(1,1))
  
  #Clean the data
  goodRows=which(fullDataSet[,substrateCol]!="#REF" & fullDataSet[,substrateCol]!="" & fullDataSet[,substrateCol]!="ID" & fullDataSet[,substrateCol]!="UniprotAccession")
  goodRowsTable=fullDataSet[goodRows,(substrateCol):ncol(fullDataSet)]
  goodRowsNames<<-goodRowsTable[,1]
  
  
  #for ggplot
  graphs=list()
  for(sub in goodRows){
    graphs[[sub-sequenceRow]]=ggplot(as.data.frame(matrix(NA, ncol=2, nrow=5)))
  }
  
  barcodeColumns5=as.data.frame(matrix(nrow=nrow(controlSamplesTest), ncol=0))
  
  #For each unique sample (control)
  results_table_c=NULL
  for(samp in 1:length(controlSamples_u)){
    
    #Pull out the correct samples and reduce this to the post wash phase
    barcodeColumns=sapply(controlSamples_u[samp], grep, controlSamplesTest)
    barcodeColumns2=controlSamplesTest[,barcodeColumns]
    postWashCol=min(which(sapply(1:(length(as.numeric(barcodeColumns2[cycleRow,]))-1), temp_fun, as.numeric(barcodeColumns2[cycleRow,]))==0))
    barcodeColumns3=barcodeColumns2[,postWashCol:ncol(barcodeColumns2)]
    barcodeColumns5=cbind(barcodeColumns5,barcodeColumns3)
    
    #Some empty vectors to store the results of the linear regression
    r2s=NULL
    slopes=NULL
    
    #For each substrate
    for(sub in goodRows){
      
      #Do and save linear regression results (maybe plot)
      #plot(as.numeric(barcodeColumns3[exposureTimeRow,]), as.numeric(barcodeColumns3[sub,]), xlab="Exposure time", ylab="Values", main=substrateNames[sub-sequenceRow])
      
      #temp_graph=ggplot(as.data.frame(cbind(as.numeric(barcodeColumns3[exposureTimeRow,]), as.numeric(barcodeColumns3[sub,]))), aes(x=V1, y=V2)) +geom_point()
      temp_graph=graphs[[sub-sequenceRow]] + geom_line(aes(x=V1, y=V2), as.data.frame(cbind(as.numeric(barcodeColumns3[exposureTimeRow,]), as.numeric(barcodeColumns3[sub,])))) + labs( x= "Time", y= "Intensity")
      graphs[[sub-sequenceRow]]=temp_graph
      df=as.data.frame(cbind(as.numeric(barcodeColumns3[exposureTimeRow,]), as.numeric(barcodeColumns3[sub,])))
      #linear=lm(df)
      linear=lm(V2~V1, df)
      #linear=lm(y~0+.,df)
      r2s=c(r2s,summary(linear)$r.squared)
      slopes=c(slopes, coef(linear)[["V1"]])
    }
    
    #If there is more than one sample then average the r2s and slopes
    # if(samp == 1){
      results_table_c=cbind(results_table_c,goodRowsTable[,1], r2s, slopes)
    # }else{
    #   results_table_c[,2]=apply(rbind(r2s, as.numeric(results_table_c[,2])),2, mean)
    #   results_table_c[,3]=apply(rbind(slopes, as.numeric(results_table_c[,3])),2, mean)
    # }
  }
  results_table_c=cbind(results_table_c[,1], apply(apply(results_table_c[,seq(2,ncol(results_table_c), 3)],2,as.numeric), 1, mean),apply(apply(results_table_c[,seq(3,ncol(results_table_c), 3)],2,as.numeric), 1, mean))
  
  
  #For each unique sample (experimental)
  results_table_e=NULL
  for(samp in 1:length(experimentalSamples_u)){
    
    #Pull out the correct samples and reduce this to the post wash phase
    barcodeColumns=sapply(experimentalSamples_u[samp], grep, experimentalSamplesTest)
    barcodeColumns2=experimentalSamplesTest[,barcodeColumns]
    postWashCol=min(which(sapply(1:(length(as.numeric(barcodeColumns2[cycleRow,]))-1), temp_fun, as.numeric(barcodeColumns2[cycleRow,]))==0))
    barcodeColumns4=barcodeColumns2[,postWashCol:ncol(barcodeColumns2)]
    barcodeColumns5=cbind(barcodeColumns5,barcodeColumns4)
    
    #Some empty vectors to store the results of the linear regression
    r2s=NULL
    slopes=NULL
    
    #For each substrate
    for(sub in goodRows){
      
      #Do and save linear regression results (maybe plot)
      #plot(as.numeric(barcodeColumns3[exposureTimeRow,]), as.numeric(barcodeColumns3[sub,]), xlab="Exposure time", ylab="Values", main=substrateNames[sub-sequenceRow])
      temp_graph=graphs[[sub-sequenceRow]] + geom_line(color="red", aes(x=V1, y=V2), as.data.frame(cbind(as.numeric(barcodeColumns4[exposureTimeRow,]), as.numeric(barcodeColumns4[sub,])))) + labs( x= "Time", y= "Intensity")
      graphs[[sub-sequenceRow]]=temp_graph
      df=as.data.frame(cbind(as.numeric(barcodeColumns4[exposureTimeRow,]), as.numeric(barcodeColumns4[sub,])))
      linear=lm(V2~V1, df)
      r2s=c(r2s,summary(linear)$r.squared)
      slopes=c(slopes, coef(linear)[["V1"]])
    }
    
    #If there is more than one sample then average the r2s and slopes
    # if(samp == 1){
      results_table_e=cbind(results_table_e,goodRowsTable[,1], r2s, slopes)
    # }else{
    #   results_table_e[,2]=apply(rbind(r2s, as.numeric(results_table_e[,2])),2, mean)
    #   results_table_e[,3]=apply(rbind(slopes, as.numeric(results_table_e[,3])),2, mean)
    # }
  }
  results_table_e=cbind(results_table_e[,1], apply(apply(results_table_e[,seq(2,ncol(results_table_e), 3)],2,as.numeric), 1, mean),apply(apply(results_table_e[,seq(3,ncol(results_table_e), 3)],2,as.numeric), 1, mean))
  
  results_table=inner_join(as.data.frame(results_table_c), as.data.frame(results_table_e), by="V1")
  if(opt1==T){ #Remove substrates without a max value of above 2 for each sample
    print("in opt 1")
    temp=matrix(ncol=length(allSamples_u), nrow=length(goodRows))
    for(i in 1:(ncol(barcodeColumns5)/5)){
      temp[,i]=as.numeric(as.character(apply(barcodeColumns5[goodRows,((1:5)+(i-1)*5)],1,max)))
    }
    results_table=results_table[which(apply(temp, 1, min)>=2),]
    x<<-length(which(apply(temp, 1, min)<2))
    #results_table=results_table[which(barcodeColumns3[sequenceRow+1:nrow(barcodeColumns3),5]>=2 & barcodeColumns4[sequenceRow+1:nrow(barcodeColumns4),5]>=2, goodRows),]
    # results_table_c=results_table_c[intersect(which(barcodeColumns3[,5]>=2), goodRows)-goodRows[1],]
    # results_table_e=results_table_e[intersect(which(barcodeColumns4[,5]>=2), goodRows)-goodRows[1],]
  }
  if(opt2==T){ #Remove substrates without an average R2 of greater than 0.9
    print("in opt 2")
    #results_table_c=results_table_c[which(results_table_c[,2]>=0.9),]   #TODO: variablize this.
    #results_table_e=results_table_e[which(results_table_e[,2]>=0.9),]
    y<<-length(which(as.numeric(as.character(results_table[,2]))<0.9 | as.numeric(as.character(results_table[,4]))<0.9))
    results_table=results_table[which(as.numeric(as.character(results_table[,2]))>=0.9 & as.numeric(as.character(results_table[,4]))>=0.9),]
    
  }
  # reduced_peptides=intersect(results_table_c[,1], results_table_e[,1])
  # temp_results_table_c=log2(as.numeric(results_table_c[results_table_c[,1] %in% reduced_peptides,3])*100)
  # temp_results_table_e=log2(as.numeric(results_table_e[results_table_e[,1] %in% reduced_peptides,3])*100)
  # temp_results_table=temp_results_table_c/temp_results_table_e
  # results_table_c[,3]=sign(as.numeric(results_table_c[,3]))*log2(abs(as.numeric(results_table_c[,3]))*100)
  # results_table_e[,3]=sign(as.numeric(results_table_e[,3]))*log2(abs(as.numeric(results_table_e[,3]))*100)
  # goodResults_c=results_table_c[which(sign(as.numeric(results_table_c[,3]))==1),]
  # goodResults_e=results_table_e[which(sign(as.numeric(results_table_e[,3]))==1),]
  # goodResults=inner_join(as.data.frame(goodResults_c), as.data.frame(goodResults_e), by="V1") #This is what I'm going to return from my function
  #goodResults[,6]=as.numeric(as.character(goodResults[,5]))-as.numeric(as.character(goodResults[,3])) #fold change differnce between log slopes
  # goodResults[,6]=as.numeric(as.character(goodResults[,5]))/as.numeric(as.character(goodResults[,3])) #fold change differnce between slopes
  # goodResults[,7]=2^goodResults[,6]
  # goodResults=goodResults[which(goodResults[,7]<=lower | goodResults[,7]>=upper),]
  results_table[,6]=as.numeric(as.character(results_table[,3]))/as.numeric(as.character(results_table[,5]))
  results_table_prime=results_table
  #TODO: save previous results table for heatmap
  results_table=results_table[which(results_table[,6]<=lower | results_table[,6]>=upper),]
  
  return(list(results_table, graphs, results_table_prime))
}

#This function is used in the post wash phase identification to find when the cycle is no longer increasing
temp_fun <- function(ind,list){
  dist(c(list[ind], list[ind+1]))
}

###############################################################################
ui <- shinyUI(fluidPage(
  
  #Application title
  titlePanel("KRSA"),
  
  #Layout type
  sidebarLayout(
    
  #Sidebar Panel
  sidebarPanel(
    
    #TODO: data upload
    fileInput('file0', label = h3("Upload raw data file"),
              accept=c('text/csv', 'text/comma-separated-values','text/plain', '.csv')),
    #Options for file type including header, seperator, and quotes.
    h4("File options:"),
    checkboxInput('header', 'Header', TRUE),
    radioButtons('sep', 'Separator', inline=F,
                 c(Tab='\t',
                   Comma=',',
                   Semicolon=';'),
                 '\t'),
    radioButtons('quote', 'Quote', inline=F,
                 c(None='',
                   'Double Quote'='"',
                   'Single Quote'="'"),
                 '"'),
    tags$hr(),
    radioButtons("radio", label = h3("Chose your kinase-peptide association option"),
                 choices = list("PiNET" = 1, "Upload" = 2), 
                 selected = 2),
    
    #This takes a tab or comma seperated file as input, with the first column "Substrates" and
    #the second column "Kinases" being the only two columns and both required.
    conditionalPanel(
      condition = "input.radio==2",
      fileInput('file1', label = h5("Kinase-peptide association file"),
                accept=c('text/csv', 'text/comma-separated-values','text/plain', '.csv')),
      h4("File options:"),
      checkboxInput('header', 'Header', TRUE),
      radioButtons('sep', 'Separator', inline=F,
                   c(Tab='\t',
                     Comma=',',
                     Semicolon=';'),
                   '\t'),
      radioButtons('quote', 'Quote', inline=F,
                   c(None='',
                     'Double Quote'='"',
                     'Single Quote'="'"),
                   '"')
    ),
  
    
    
    #Options for file type including header, seperator, and quotes.
    
    tags$hr()
    

    
    #Testing the chooser
    #chooserInput("mychooser", "Available frobs", "Selected frobs",
    #             row.names(USArrests), c(), size = 10, multiple = TRUE)
  ),
  
  #Main Panel
  mainPanel(
    tabsetPanel(
      tabPanel("Step 1: Options",
               dataTableOutput('contents0'),
               h4("Data input options:"),
               #celltype
               fluidRow(
                 column(4,
                        dropdownButton(
                          label = "Step 1: Experiment Barcode", status = "default", width = 80,
                          actionButton(inputId = "celltypea2z", label = "Sort A to Z", icon = icon("sort-alpha-asc")),
                          actionButton(inputId = "celltypez2a", label = "Sort Z to A", icon = icon("sort-alpha-desc")),
                          br(),
                          actionButton(inputId = "celltypeall", label = "(Un)select all"),
                          checkboxGroupInput(inputId = "celltype", label = "Choose", choices = c()))),
                 column(4,
                        dropdownButton(
                          label = "Step 2: Control Samples", status = "default", width = 80,
                          actionButton(inputId = "controlSamplesa2z", label = "Sort A to Z", icon = icon("sort-alpha-asc")),
                          actionButton(inputId = "controlSamplesz2a", label = "Sort Z to A", icon = icon("sort-alpha-desc")),
                          br(),
                          actionButton(inputId = "controlSamplesall", label = "(Un)select all"),
                          checkboxGroupInput(inputId = "controlSamples", label = "Choose", choices = c()))),
                 column(4,
                        dropdownButton(
                          label = "Step 3: Experimental Samples", status = "default", width = 80,
                          actionButton(inputId = "experimentalSamplesa2z", label = "Sort A to Z", icon = icon("sort-alpha-asc")),
                          actionButton(inputId = "experimentalSamplesz2a", label = "Sort Z to A", icon = icon("sort-alpha-desc")),
                          br(),
                          actionButton(inputId = "experimentalSamplesall", label = "(Un)select all"),
                          checkboxGroupInput(inputId = "experimentalSamples", label = "Choose", choices = c())))
                 ),
               h4("Control and Experimental Tables:"),
               fluidRow(
                 column(6,
                        DT::dataTableOutput('contents1')
                        ),
                 column(6,
                        DT::dataTableOutput('contents2')
                        )
               ),
               tags$hr()
             ),
      tabPanel("Step 2: Stringency",
               sliderInput("FCslider", label = h3("Select Fold Change"), min = 0,  max = 2, value = c(0.85,1.15), step=0.05),
               h3("Select Stringency Level"),
               fluidRow(
                 column(6,
                        checkboxGroupInput("StrinRadio",label=NULL, choices = list("Max exposure >= 2" = 1, "R^2 >=0.9" = 2), selected = c(1,2)),
                        actionButton("calculateData", "Find Significant Substrates")
                        ),
                 column(6,
                        textOutput("remain1"),
                        textOutput("remain2")
                        )
               ),
               hr(),
               selectInput("xyplot", "Choose a Substrate:", choices=c()), 
               plotOutput("xyplot2", width = "100%", height="600px"),
               hr()
               #dataTableOutput("subTable")

      ),
      tabPanel("Step 3: Iteration", 
               #This is the slider that allows you to choose the number of iterations of sampling you
               #want to do, from 0 to 5000 (1000 default) and jumping in increments of 100.
               sliderInput("slider1", label = h4("Select the number of iterations"), min = 0, 
                           max = 5000, value = 1000, step=100, width='100%'),
               tags$hr(),
               h4("Confirm and Run"),
               textOutput("subnum"),
               textOutput("subnum4"),
               textOutput("subnum3"),
               htmlOutput("subnum2"),
               actionButton("compute", "Run KRSA")

      ),
      tabPanel("Step 4: Results",
           dataTableOutput("overallTable")
      ),
      tabPanel("Step 5: Histograms",
          selectInput("kinase", "Choose a Kinase:", choices=c()), 
          plotOutput("histo", width = "100%", height="600px"),
          downloadButton('downloadData', 'Download Histogram')
      )
    )

  )
  )
))


###############################################################################
server <- function(input, output, session) {
  
  overallTable2 <<- reactiveValues()
  scores <- reactiveValues()
  barcodeOptions2 <<-reactiveValues()
  controlOptions2 <<-reactiveValues()
  experimentalOptions2 <<-reactiveValues()
  graphs<<-reactiveValues()
  
  #Creates an empty object to store the Substrate names in.
  dsnames <- c()
  
  #Creates an empty object to store the selected substrates in.
  #Not sure the below is correct...
  selectedSubnames<-reactive({input$inCheckboxGroup})
  #print(selectedSubnames())
  
  
  #Reactive dataset for the input file 0 (raw data).
  data_set0 <- reactive({
    req(input$file0)
    inFile <- input$file0
    data_set0<-read.csv(inFile$datapath, header=input$header, 
                       sep=input$sep, quote=input$quote, stringsAsFactors = F)
    fullDataSet<<-data_set0
  })
  
  #Output for rendering the preview table in the main panel.
  output$contents0 <- DT::renderDataTable({
    if(is.null(input$celltype)){
      return(data_set0())
    }else{
      barcodeColumns=sapply(input$celltype, grep, data_set0())
      barcodeColumns2=data_set0()[,barcodeColumns]
      return(barcodeColumns2)
    }
  }, caption = "Preview of raw data file")
  
  #Output for rendering the preview table in the main panel.
  output$contents1 <- DT::renderDataTable({
    if(is.null(input$controlSamples)){
      return(data_set0())
    }else{
      a=paste0("^", input$controlSamples, "$")
      b=NULL
      for(x in 1:length(a)){
        temp=as.numeric(which(apply(data_set0(), 2, grep, pattern=a[x], value=FALSE, perl=TRUE)!=0))
        print(temp)
        b=c(b, temp)
      }
      barcodeColumns=unlist(b)
      print(barcodeColumns)
      #barcodeColumns=sapply(input$controlSamples, grep, data_set0())
      barcodeColumns2=data_set0()[,unlist(barcodeColumns)]
      controlSamplesTest<<-barcodeColumns2
      return(barcodeColumns2)
    }
  }, caption = "Preview of control data", options=list(autoWidth=TRUE))
  
  #Output for rendering the preview table in the main panel.
  output$contents2 <- DT::renderDataTable({
    if(is.null(input$experimentalSamples)){
      return(data_set0())
    }else{
      a=paste0("^", input$experimentalSamples, "$")
      print(a)
      b=NULL
      for(x in 1:length(a)){
        b=c(b, as.numeric(which(apply(data_set0(), 2, grep, pattern=a[x], value=FALSE, perl=TRUE)!=0)))
      }
      barcodeColumns=unlist(b)
      print(barcodeColumns)
      #barcodeColumns=sapply(input$experimentalSamples, grep, data_set0())
      barcodeColumns2=data_set0()[,unlist(barcodeColumns)]
      experimentalSamplesTest<<-barcodeColumns2
      return(barcodeColumns2)
    }
  }, caption = "Preview of experimental data", options=list(autoWidth=TRUE))

  #barcode select
  observe({
    barcodeOptions=grep("Barcode", data_set0())
    barcodeOptions1=which(data_set0()[,barcodeOptions]=="Barcode")
    barcodeOptions1.5=as.numeric(data_set0()[barcodeOptions1,(barcodeOptions+1):ncol(data_set0())])
    barcodeOptions2<<-unique(sort(barcodeOptions1.5))
    updateCheckboxGroupInput(session, "celltype",
                             choices = barcodeOptions2,
                             selected = NA #This change made it so all of the boxes were selected by default
    )
  })
  # Sorting asc
  observeEvent(input$celltypea2z, {
    updateCheckboxGroupInput(
      session = session, inputId = "celltype", 
      choices = barcodeOptions2,
      selected = input$celltype
    )
  })
  # Sorting desc
  observeEvent(input$celltypez2a, {
    updateCheckboxGroupInput(
      session = session, inputId = "celltype", 
      choices = sort(barcodeOptions2,decreasing=T),
      selected = input$celltype
    )
  })
  # Select all / Unselect all
  observeEvent(input$celltypeall, {
    if (is.null(input$celltype)) {
      updateCheckboxGroupInput(
        session = session, inputId = "celltype", 
        selected = barcodeOptions2
      )
    } else {
      updateCheckboxGroupInput(
        session = session, 
        inputId = "celltype", 
        selected = ""
      )
    }
  })
  
  #control group select
  observeEvent(input$celltype,{
    controlOptions=grep("Sample name", data_set0())
    controlOptions1=which(data_set0()[,controlOptions]=="Sample name")
    controlOptions1.5=as.character(data_set0()[controlOptions1,sapply(input$celltype, grep, data_set0())])
    controlOptions2<<-unique(sort(controlOptions1.5))
    #print(controlOptions)
    #print(controlOptions1)
    #print(controlOptions1.5)
    #print(controlOptions2)
    updateCheckboxGroupInput(session, "controlSamples",
                             choices = controlOptions2,
                             selected = NA #This change made it so all of the boxes were selected by default
    )
  })
  # Sorting asc
  observeEvent(input$controlSamplesa2z, {
    updateCheckboxGroupInput(
      session = session, inputId = "controlSamples", 
      choices = controlOptions2,
      selected = input$controlSamples
    )
  })
  # Sorting desc
  observeEvent(input$controlSamplesz2a, {
    updateCheckboxGroupInput(
      session = session, inputId = "controlSamples", 
      choices = sort(controlOptions2,decreasing=T),
      selected = input$controlSamples
    )
  })
  # Select all / Unselect all
  observeEvent(input$controlSamplesall, {
    if (is.null(input$controlSamples)) {
      updateCheckboxGroupInput(
        session = session, inputId = "controlSamples", 
        selected = controlOptions2
      )
    } else {
      updateCheckboxGroupInput(
        session = session, 
        inputId = "controlSamples", 
        selected = ""
      )
    }
  })
  
  #experimental group select
  observeEvent(input$controlSamples,{
    experimentalOptions=grep("Sample name", data_set0())
    experimentalOptions1=which(data_set0()[,experimentalOptions]=="Sample name")
    experimentalOptions1.5=as.character(data_set0()[experimentalOptions1,sapply(input$celltype, grep, data_set0())])
    experimentalOptions2<<-unique(sort(experimentalOptions1.5))
    #print(experimentalOptions)
    #print(experimentalOptions1)
    #print(experimentalOptions1.5)
    #print(experimentalOptions2)
    updateCheckboxGroupInput(session, "experimentalSamples",
                             choices = experimentalOptions2,
                             selected = NA #This change made it so all of the boxes were selected by default
    )
  })
  # Sorting asc
  observeEvent(input$experimentalSamplesa2z, {
    updateCheckboxGroupInput(
      session = session, inputId = "experimentalSamples", 
      choices = experimentalOptions2,
      selected = input$experimentalSamples
    )
  })
  # Sorting desc
  observeEvent(input$experimentalSamplesz2a, {
    updateCheckboxGroupInput(
      session = session, inputId = "experimentalSamples", 
      choices = sort(experimentalOptions2,decreasing=T),
      selected = input$experimentalSamples
    )
  })
  # Select all / Unselect all
  observeEvent(input$experimentalSamplesall, {
    if (is.null(input$experimentalSamples)) {
      updateCheckboxGroupInput(
        session = session, inputId = "experimentalSamples", 
        selected = experimentalOptions2
      )
    } else {
      updateCheckboxGroupInput(
        session = session, 
        inputId = "experimentalSamples", 
        selected = ""
      )
    }
  })
  #Reactive dataset for the input file 1  (kinase peptide associations).
  data_set <- reactive({
    req(input$file1)
    inFile <- input$file1
    data_set<-read.csv(inFile$datapath, header=input$header, 
                       sep=input$sep, quote=input$quote)
  })
  
  #Output for rendering the preview table in the main panel.
  output$contents <- renderTable({
    data_set()
  })
  
  #For the stringency tab
  output$FCslider <- renderPrint({ input$FCslider })
  output$StrinRadio <- renderPrint({ input$StrinRadio })
  
  
  #Display graph
  observeEvent(input$calculateData, {
    #print(input$xyplot)
    output$xyplot2 <- renderPlot({
      plotIndex=which(goodRowsNames==input$xyplot)
      print(plotIndex)
      graphs[plotIndex]
    })
  })
  
  #Run make_graphs
  observeEvent(input$calculateData, {
    options1=input$StrinRadio
    print(options1)
    if(1 %in% options1){
      opt1=T
    }else{
      opt1=F
    }
    if(2 %in% options1){
      opt2=T
    }else{
      opt2=F
    }
    results<<-make_graphs(controlSamplesTest, experimentalSamplesTest, fullDataSet, opt1, opt2, input$FCslider[1], input$FCslider[2])
    #print(class(results))
    resultsTable<<-results[[1]]
    graphs <<- results[[2]]
    resultsTablePrime<<-results[[3]]
    print(class(graphs))
    updateSelectInput(
      session = session, inputId = "xyplot", 
      choices = resultsTable[,1],
      selected = resultsTable[1,1]
    )
    output$overallTable <- DT::renderDataTable({
      datatable(subTable, rownames = TRUE, caption = "Substrates Table")
    })
  })
  
  #Guts for making the dynamic checkboxes work.
  observe({
    req(input$file1)
    dsnames <- data_set()$Substrates
    updateCheckboxGroupInput(session, "inCheckboxGroup",
                             label = "",
                             choices = dsnames,
                             selected = "")
  })
  
  #You can access the value of the widget with input$slider1, e.g.
  #output$value <- renderPrint({ input$slider1 })
  
  #State how many of each element there are.
  observeEvent(input$calculateData,{
    
    #Update checkboxes to reflect the number of subs removed by each option
    if(1 %in% input$StrinRadio & 2 %in% input$StrinRadio){
      kittyx=paste0("Max exposure >= 2 (", x, " substrates removed with Option 1)")
      kittyy=paste0("R^2 >=0.9 (", y, " substrates removed with Option 2)")
      updateCheckboxGroupInput(session, inputId="StrinRadio", choiceNames = list(kittyx, kittyy), choiceValues=c(1,2), selected = c(1,2))
    }else if(1 %in% input$StrinRadio){
      kittyx=paste0("Max exposure >= 2 (", x, " substrates removed with Option 1)")
      updateCheckboxGroupInput(session, inputId="StrinRadio", choiceNames = list(kittyx, "R^2 >=0.9"), choiceValues=c(1,2), selected = 1)
    }else if(2 %in% input$StrinRadio){
      kittyy=paste0("R^2 >=0.9 (", y, " substrates removed with Option 2)")
      updateCheckboxGroupInput(session, inputId="StrinRadio", choiceNames = list("Max exposure >= 2", kittyy), choiceValues=c(1,2), selected = 2)
    }
    
    output$remain1 <- renderText({
      paste0("You have selected ", nrow(resultsTable), " substrates based on inclusion and exclusion criteria.")
    })
    
    output$remain2 <- renderText({
      paste0( nrow(resultsTable), "/", nrow(resultsTablePrime), " (", round((nrow(resultsTable)/nrow(resultsTablePrime))*100, 2), "%) substrates were chosen for KRSA analysis.")
    })
    
    output$subnum <- renderText({
      paste0("There are ", nrow(data_set()), " substrates in your data set.")
    })
    
    output$subnum3 <- renderText({
      paste0("You have selected ", nrow(resultsTable) , " significant substrates from your data set (listed below):")
    })
    
    output$subnum2 <- renderUI({
      str1 <- paste0("")
      str2 <- as.character(resultsTable[,1])
      HTML(paste(str2, str1, sep = '<br/>'))
    })
    
    output$subnum4 <- renderText({
      paste0("You have selected ", input$slider1 , " iterations.")
    })
  })
  
  observeEvent(input$compute, {
    functionresults=krsa(input$file1$datapath,iterations=input$slider1, numSigSubs=nrow(resultsTable), sigSubs=as.character(resultsTable[,1]))
    showNotification("KRSA Finished")
    overallTable2<<-functionresults[[1]]
    meow<<-as.data.frame(overallTable2)
    scores<<-functionresults[[2]]
    output$overallTable <- DT::renderDataTable({
      datatable(overallTable2, rownames = TRUE, caption = "Results Table")
    })
    #Guts for making the histogram drop down work.
    updateSelectInput(session, "kinase", choices=rownames(overallTable2))
  })
  
  plotInput <-function(){
    OTindex<<-which(rownames(overallTable2)==input$kinase)
    Sindex<<-which(rownames(scores)==input$kinase)
    kitty<<-as.matrix(scores[Sindex,1:input$slider1])
    kitty2=matrix(ncol=2, nrow=(nrow(resultsTable)+1))
    kitty2[,1]=0:nrow(resultsTable)
    for(i in 0:nrow(resultsTable)){
      kitty2[i+1,2]=length(which(kitty==i))
    }
    
    theData=as.data.frame(kitty2)
    colnames(theData)=c("Hits", "Frequency")
    SDrect <- data.frame (xmin=overallTable2[OTindex,6], xmax=overallTable2[OTindex,7], ymin=0, ymax=Inf)
    ggplot(theData, aes(Hits,Frequency)) +
      geom_bar(stat="identity") +
      scale_x_continuous(breaks=seq(0,nrow(resultsTable),ceiling(nrow(resultsTable)/20))) +
      geom_rect(data=SDrect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray", alpha=0.5, inherit.aes = FALSE) +
      geom_segment(x=overallTable2[OTindex,6], y=0, xend=overallTable2[OTindex,6], yend=Inf, colour="cyan") +
      geom_segment(x=overallTable2[OTindex,7], y=0, xend=overallTable2[OTindex,7], yend=Inf, colour="cyan") +
      geom_segment(x=overallTable2[OTindex,2], y=0, xend=overallTable2[OTindex,2], yend=Inf, colour="red") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }
  
  output$histo <- renderPlot({
    plotInput()
  })
  
  #Download Histogram
  output$downloadData <- downloadHandler(
    filename = 'test.pdf', #This is necessary even though it doesn't save it as 'test.pdf', you have to name it 'something.pdf' instead.
    content = function(file) {
      pdf(file,width=7,height=10)
      plotInput()
      dev.off()
    })
  
}

###############################################################################
#Make the app!
shinyApp(ui, server)
