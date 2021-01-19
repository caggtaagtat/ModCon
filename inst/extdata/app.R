#Load ggplot package to render the plots
library("ggplot2")
library("shinycssloaders")
library("shiny")
library("shinyFiles")
library("seqinr")
library("Biostrings")
library("BSgenome")
library("VarCon")

readRDS("exampleTransCoord")

###########################################
## Server and UI code

## Server in and output Skript
server <- function(input, output, session) {
  
  ## close the R session when app closes
  session$onSessionEnded(function() {
    stopApp()
  })
  
  uploadReferenceDNA <- eventReactive(path$pth,{
    
    testFASTA <- strsplit(path$pth,"\\.")[[1]]
    if(testFASTA[length(testFASTA)] %in% c("fa","fasta")){
      referenceDnaStringSet2 <- readDNAStringSet(path$pth, format="fasta",use.names=TRUE)
      ref_names <- as.character(lapply(names(referenceDnaStringSet2),
                                       function(x){ strsplit(x, " ")[[1]][[1]]}))
      names(referenceDnaStringSet2) <- ref_names
      referenceDnaStringSet <- referenceDnaStringSet2
    }else{
      load(path$pth)
    }   
    
    referenceDnaStringSet
  }) 
  
  
  
  uploadTranscriptTable <- eventReactive(path3$pth3,{
    
    ## Get human transcript tables e.g. from https://github.com/caggtaagtat/VarConTables
    testCSV <- strsplit(path3$pth3,"\\.")[[1]]
    if(testCSV[length(testCSV)] == "csv"){
      transCoord <- read.csv(path3$pth3, sep=";")
    }else{ transCoord <- readRDS(path3$pth3)}
    
    transCoord
  }) 
  
  ## Report when upload of reference genome is complete
  output$sum_text2 <- renderUI({
    
    ## Genome fasta to download e.g. from 
    ## ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/
    ## dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
    
    test <- uploadReferenceDNA()
    HTML("Upload of reference genome completed...")
    
  })
  
  ## Report when upload of reference genome is complete
  output$sum_text22 <- renderUI({
    
    test2 <- uploadTranscriptTable()
    HTML("Upload of transcript table completed...")
    
  })
  
  
  
  ## Generate the text for describing the difference between the Hexplorer Scores of both sequences
  output$sum_text <- renderUI({
    
    gene2transcript <- read.csv(path2$pth2, sep=";", stringsAsFactors=FALSE)
    referenceDnaStringSet <- uploadReferenceDNA()
    transCoord <- uploadTranscriptTable()
    
    #Get information about the SNV
    res <-  getSeqInfoFromVariation(referenceDnaStringSet, input$transcriptID, input$variation,
                                    ntWindow= input$ntWindow, transCoord, gene2transcript=gene2transcript)
    
    #Sum up the info
    HTML(paste0("For the given annotation  ",res$funcAnnotation,
                      " within transcript ", res$transcript,
                      " following sequence was found around the chromosomal coordinate ",
                      res$genomicCoordinate," on chromosome ",res$chromosome, " :"),
               "",paste0("Ref Seq: ",res$sequence),"", 
               paste0("Ref+vari:",res$altSeq) , sep="<br/>")
    
  })
  
  ## Generate the plot where you can mark the area to zoom in
  output$plot <- renderPlot({
    
    gene2transcript <- read.csv(path2$pth2, sep=";", stringsAsFactors=FALSE)
    referenceDnaStringSet <- uploadReferenceDNA()
    transCoord <- uploadTranscriptTable()
    
    ## Retrieve information form genome
    res <-  getSeqInfoFromVariation(referenceDnaStringSet, input$transcriptID,
                                    input$variation, ntWindow=input$ntWindow, transCoord,
                                    gene2transcript=gene2transcript)
    
    ## Calculate HZEI values
    durchzahl <-  calculateHZEIperNT(res$sequence)
    
    durchzahl$Sequence <- "sequence of interest"
    
    plot <- ggplot(durchzahl, aes(x = durchzahl, y = endhex, fill=Sequence )) +scale_y_continuous(name="Hexplorer score",breaks=c(seq(-75,0,5),seq(2,34,2)),limits=c(min(durchzahl$endhex)-5,max(c(durchzahl$hbs,durchzahl$endhex))+1) )+
      scale_fill_manual(values=c("#56B4E9", "#000000"))+
      geom_bar(stat='identity', position = "dodge")+ xlab("Sequence")+ylab("Hexplorer score")+
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
      annotate("text", label =substr(durchzahl$seq9[durchzahl$Sequence=="sequence of interest"],6,6 ) , x= 1:((nrow(durchzahl))), y = min(durchzahl$endhex-8), size = 3, colour = "black")
    
    plot
    
  })
  
  ## Generate the plot where you can mark the area to zoom in
  output$plot2 <- renderPlot({
    
    gene2transcript <- read.csv(path2$pth2, sep=";", stringsAsFactors=FALSE)
    referenceDnaStringSet <- uploadReferenceDNA()
    transCoord <- uploadTranscriptTable()
    
    ## Retrieve information form genome
    res <-  getSeqInfoFromVariation(referenceDnaStringSet, input$transcriptID, 
                                    input$variation, ntWindow=input$ntWindow, transCoord,
                                    gene2transcript=gene2transcript)
    
    generateHEXplorerPlot(res,input$ntWindow)
    
  })
  
  ## Create reactive value ranges, for the zooming plot
  ranges2 <- reactiveValues(x = NULL)
  
  ## Genereate the plot, where you can see a zoomed in version of the plot above
  output$plot_zoom <- renderPlot({
    
    gene2transcript <- read.csv(path2$pth2, sep=";", stringsAsFactors=FALSE)
    referenceDnaStringSet <- uploadReferenceDNA()
    transCoord <- uploadTranscriptTable()
    
    
    ## Retrieve information form genome
    res <-  getSeqInfoFromVariation(referenceDnaStringSet, input$transcriptID,
                                    input$variation, ntWindow=input$ntWindow,
                                    transCoord,gene2transcript=gene2transcript)
    
    
    ## Calculate HZEI values
    durchzahl <-  calculateHZEIperNT(res$sequence)
    
    durchzahl$Sequence <- "sequence of interest"
    
    plot <- ggplot(durchzahl, aes(x = durchzahl, y = endhex, fill=Sequence )) +scale_y_continuous(name="Hexplorer score",breaks=c(seq(-75,0,5),seq(2,34,2)),limits=c(min(durchzahl$endhex)-6,max(c(durchzahl$hbs,durchzahl$endhex))+1) )+
      scale_fill_manual(values=c("#56B4E9", "#000000"))+
      geom_bar(stat='identity', position = "dodge")+ xlab("Sequence")+ylab("Hexplorer score")+
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
      annotate("text", label =substr(durchzahl$seq9[durchzahl$Sequence=="sequence of interest"],6,6 ) , x= 1:((nrow(durchzahl))), y = min(durchzahl$endhex-4), size = 3, colour = "black")+
      coord_cartesian(xlim = ranges2$x, expand = FALSE)
    plot
    
    
    
  })
  
  ## Genereate the plot, where you can see a zoomed in version of the plot above
  output$plot2_zoom <- renderPlot({
    
    gene2transcript <- read.csv(path2$pth2, sep=";", stringsAsFactors=FALSE)
    referenceDnaStringSet <- uploadReferenceDNA()
    transCoord <- uploadTranscriptTable()
    
    
    ## Retrieve information form genome
    res <-  getSeqInfoFromVariation(referenceDnaStringSet, input$transcriptID, 
                                    input$variation, ntWindow=input$ntWindow, transCoord,
                                    gene2transcript=gene2transcript)
    
    results_plot <- generateHEXplorerPlot(res,input$ntWindow)
    results_plot+coord_cartesian(xlim = ranges2$x, expand = FALSE)
    
  })
  
  ## Create a function which keeps checking on the input with the mouse
  observe({
    brush <- input$plot2_brush
    if (!is.null(brush)) {
      ranges2$x <- c(brush$xmin, brush$xmax)
      
    } else {
      ranges2$x <- NULL
    }
  })
  
  ## Generate the text for describing the difference between the Hexplorer Scores of both sequences
  output$plot2_text <- renderUI({
    
    gene2transcript <- read.csv(path2$pth2, sep=";", stringsAsFactors=FALSE)
    referenceDnaStringSet <- uploadReferenceDNA()
    transCoord <- uploadTranscriptTable()
    
    ## Retrieve information form genome
    res <-  getSeqInfoFromVariation(referenceDnaStringSet, input$transcriptID,
                                    input$variation, ntWindow=input$ntWindow,
                                    transCoord,gene2transcript=gene2transcript)
    
    ## calculte HZEI values
    durchzahl <-  calculateHZEIperNT(res$sequence)
    
    durchzahl$Sequence <- "reference"
    
    durchzahl2 <-  calculateHZEIperNT(res$altSeq)
    durchzahl2$Sequence <- "alternative"
    
    durchzahl$mut_hex  <- durchzahl2$endhex
    durchzahl$diff_hex <- durchzahl$mut_hex- durchzahl$endhex
    
    pre <- ""
    if(sum(durchzahl$diff_hex) > 0) pre <- "+"
    
    ## Return the Difference in the score
    HTML(paste0("The difference in the HEXplorer Score integral amounts to ", paste0(pre, sum(durchzahl$diff_hex)), " in total."))
    
  })
  
  
  ## Define reactive paths
  path <- reactiveValues(
    pth= system.file("extdata", "referenceDnaStringSet.fa", package="VarCon")
  )
  
  
  path2 <- reactiveValues(
    pth2= system.file("extdata", "fastaEx.fa", package="Biostrings")
  )
  
  path3 <- reactiveValues(
    pth3= system.file("extdata", "exampleTransCoord", package="VarCon")
  )
  
  
  
  observeEvent(input$filechoose,{
    path$pth <- file.choose()
  })
  
  observeEvent(input$filechoose2,{
    path2$pth2 <- file.choose()
  })
  
  observeEvent(input$filechoose3,{
    path3$pth3 <- file.choose()
  })
  
  
  
}

#User Interface Script

ui <- fluidPage(
  
  ## Type Headline
  titlePanel("VarCon: Retrieve genomic sequence around sequence variation"),
  
  "VarCon retrieves the surrounding genomic sequence of a stated sequence variation and visualizes potential changes in sequence elements important for splicing. Please first upload the fasta file of the respective reference genome sequence. Loading and processing of the data will take up to 2 minutes.",
  
  br(),
  br(),
  
  ## Have different tabs in your programm
  tabsetPanel(type = "tabs",
              
              tabPanel("Upload reference data",
                       
                       
                       fluidRow(
                         
                         column(4,  h4("Fasta reference genome"),
                                actionButton("filechoose",label = "Select FASTA file")
                                
                                
                         ),
                         
                         
                         column(4, h4("Transcript table"),
                                actionButton("filechoose3",label = "Select transcript table")),
                         
                         column(3,  h4("Optional: gene/transcript table"),
                                actionButton("filechoose2",label = "Select gene2transcript table")
                                
                                
                         )
                         
                       ),
                       
                       
                       fluidRow(
                         
                         column(4,   withSpinner(htmlOutput("sum_text2"), type=6)),
                         column(4,   withSpinner(htmlOutput("sum_text22"), type=6))
                         
                         
                         
                       )
                       
                       
                       
              ),
              
              tabPanel("Retrieve sequence around SNV",
                       
                       
                       fluidRow(
                         
                         column(3,
                                h4("Please enter the required information"),
                                helpText("Please enter the transcript of interrest",
                                         "an the annotation of the functional variation."
                                )),
                         
                         
                         column(3, textInput("transcriptID", label = h4("Transcript ID (ENSEMBL)"),value= "ENST00000544455")),
                         
                         column(3, textInput("variation", label = h4("Functional variation"),value= "c.516+21A>T")),
                         
                         column(3, numericInput("ntWindow", label = h4("Seq x nt up/downstream"), value= 20, min=5, max=150))
                         
                         
                         
                       ),
                       br(),
                       
                       
                       fluidRow(
                         
                         column(12, withSpinner(htmlOutput("sum_text"), type=6))
                         
                         
                       ),
                       
                       br()
                       
                       
                       
                       
                       
                       
              ),
              
              tabPanel("Impact splice site strength and SREs",
                       
                       br(),
                       
                       br(),
                       
                       withSpinner(plotOutput("plot2", height = 200,
                                              brush = brushOpts(
                                                id = "plot2_brush",
                                                resetOnNew = TRUE
                                              )), type =6),
                       
                       h4("Zoomed in plot:"),
                       plotOutput("plot2_zoom", height = 200),
                       
                       br(),
                       
                       
                       htmlOutput("plot2_text")
                       
                       
              ),
              
              tabPanel("Manual",
                       h3("1. Upload reference genome fasta file"),
                       "First, please upload the fasta file (or zipped fasta.gz) of the reference genome sequence.",
                       br(),
                       "Potentially required data for the reference genome GRCh37 and GRCh38 is availible in the directory of this application.",
                       br(),
                       "If needed, the required file can be downloaded from the Ensembl ftp server ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
                       br(),
                       br(),
                       
                       h3("2. Select genome assembly"), 
                       "Next, select whether the uploaded genome reference file originated from assembly GRCh37 or GRCh38.",
                       "The respective transcript table, holding the genomic exon coordinates will be selected.",
                       br(),
                       br(),
                       
                       h3("3. Select gene to transcript conversion table (optional)"), 
                       "Select a csv-table holding gene names and gene transcripts which shall be used synonymously during the querries.",
                       br(),
                       br(),
                       
                       h3("4. Data entry"), 
                       "On the next panel, now enter the respective transcript name (or gene name) and the single nucleotide variation of interest.",
                       "The sequence variations can either refer to the nucleotide positions within the coding sequence or genomic coordiantes.",
                       "Example variations: c.142+2A>T  or g.12746124G>A",
                       br(),
                       br()
                       
              )
              
  )
  
)


shinyApp(ui = ui, server = server)
