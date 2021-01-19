#Load ggplot package to render the plots
library("shinycssloaders")
library("shiny")
library("shinyFiles")
library("ModCon")
library("shinydashboard")
library("shinyjs")
library("dplyr")


###########################################
#+++Server and UI code

#User Interface Script

header <- dashboardHeader(title = "ModCon")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Set general parameters", tabName = "model", icon = icon("list-alt"),
             sliderInput("optiRate",
                         "Extent of optimization [%]", min=0, max=100,
                         step=1, value=100),
             checkboxInput("increaseHZEI", "Optimize immediate context",
                           value = TRUE,
                           width = NULL),
             numericInput("upchangecodons",
                          "Exchange X codons up- and downstream",
                          16),
             numericInput("sdMaximalHBS",
                          "SD degradation treshold [HBond score]",
                          14),
             numericInput("acMaximalMaxent",
                          "SA degradation treshold [MaxEntScan score]",
                          4)))
)

body <- dashboardBody(
  
  fluidRow(shinyjs::useShinyjs(),
           tabBox(width = 12, height = NULL,
                  
                  tabPanel("Adjust SD context",
                           
                           br(),
                           fluidRow(
                             column(6,
                                    textInput("enteredCDS",
                                              "Enter coding sequence",
                                              "ATGGAAGACGCCAAAAACATAAAGAAAGGCCCGGCGCCATTCTATCCGCTGGAAGATGGAACCGCTGGAGAGCAACTGCATAAGGCTATGAAGAGATACGCCCTGGTTCCTGGAACAATTGCTTTTACAGATGCACATATCGAGGTGGACATCACTTACGCTGAGTACTTCGAAATGTCCGTTCGGTTGGCAGAAGCTATGAAACGATATGGGCTGAATACAAATCACAGAATCGTCGTATGCAGTGAAAACTCTCTTCAATTCTTTATGCCGGTGTTGGGCGCGTTATTTATCGGAGTTGCAGTTGCGCCCGCGAACGACATTTATAATGAACGTGAATTGCTCAACAGTATGGGCATTTCGCAGCCTACCGTGGTGTTCGTTTCCAAAAAGGGGTTGCAAAAAATTTTGAACGTGCAAAAAAAGCTCCCAATCATCCAAAAAATTATTATCATGGATTCTAAAACGGATTACCAGGGATTTCAGTCGATGTACACGTTCGTCACATCTCATCTACCTCCCGGTTTTAATGAATACGATTTTGTGCCAGAGTCCTTCGATAGGGACAAGACAATTGCACTGATCATGAACTCCTCTGGATCTACTGGTCTGCCTAAAGGTGTCGCTCTGCCTCATAGAACTGCCTGCGTGAGATTCTCGCATGCCAGAGATCCTATTTTTGGCAATCAAATCATTCCGGATACTGCGATTTTAAGTGTTGTTCCATTCCATCACGGTTTTGGAATGTTTACTACACTCGGATATTTGATATGTGGATTTCGAGTCGTCTTAATGTATAGATTTGAAGAAGAGCTGTTTCTGAGGAGCCTTCAGGATTACAAGATTCAAAGTGCGCTGCTGGTGCCAACCCTATTCTCCTTCTTCGCCAAAAGCACTCTGATTGACAAATACGATTTATCTAATTTACACGAAATTGCTTCTGGTGGCGCTCCCCTCTCTAAGGAAGTCGGGGAAGCGGTTGCCAAGAGGTTCCATCTGCCAGGTATCAGGCAAGGATATGGGCTCACTGAGACTACATCAGCTATTCTGATTACACCCGAGGGGGATGATAAACCGGGCGCGGTCGGTAAAGTTGTTCCATTTTTTGAAGCGAAGGTTGTGGATCTGGATACCGGGAAAACGCTGGGCGTTAATCAAAGAGGCGAACTGTGTGTGAGAGGTCCTATGATTATGTCCGGTTATGTAAACAATCCGGAAGCGACCAACGCCTTGATTGACAAGGATGGATGGCTACATTCTGGAGACATAGCTTACTGGGACGAAGACGAACACTTCTTCATCGTTGACCGCCTGAAGTCTCTGATTAAGTACAAAGGCTATCAGGTGGCTCCCGCTGAATTGGAATCCATCTTGCTCCAACACCCCAACATCTTCGACGCAGGTGTCGCAGGTCTTCCCGACGATGACGCCGGTGAACTTCCCGCCGCCGTTGTTGTTTTGGAGCACGGAAAGACGATGACGGAAAAAGAGATCGTGGATTACGTCGCCAGTCAAGTAACAACCGCGAAAAAGTTGCGCGGAGGAGTTGTGTTTGTGGACGAAGTACCGAAAGGTCTTACCGGAAAACTCGACGCAAGAAAAATCAGAGAGATCCTCATAAAGGCCAAGAAGGGCGGAAAGATCGCCGTG")
                             ),
                             column(6,
                                    numericInput("sdPosition", "Position of first SD nucleotide",
                                                 1001))
                           ),br(),
                           
                           ## Execute HEXCO
                           actionButton("goButton", strong("Optimize context!")), br(), br(), br(),
                           
                           
                           br(),
                           strong("Results"),
                           htmlOutput("quick_info")%>% withSpinner(color="#0dc5c1")
                           
                           , br(),
                           downloadButton("DownloadData", label = strong("Download optimized sequence"))
                           
                           
                  ),
                  
                  tabPanel(title="Set GA parameter (optional)",
                           
                           h4("Genetic algorithm for HZEI manipulation"),
                           
                           fluidRow(
                             column(6,
                                    numericInput("nGenerations",
                                                 "Max # generations",
                                                 20),
                                    numericInput("parentSize",
                                                 "# parental seqs",
                                                 200),
                                    numericInput("startParentSize",
                                                 "# parental seqs",
                                                 800),
                                    numericInput("bestRate",
                                                 "X% of fittest reproduce",
                                                 50)
                             ),
                             column(6,
                                    
                                    numericInput("semiLuckyRate",
                                                 "X% fitness-based chance",
                                                 20),
                                    numericInput("luckyRate",
                                                 "X% randomly reproduce",
                                                 5),
                                    numericInput("mutationRate",
                                                 "Rate of codon mutation (0-1)",
                                                 0.01, min=0, max=1)
                             )
                           )
                  )
           )
  )
)

ui <- dashboardPage(header,
                    sidebar,
                    body)



# Server in and output Skript
server <- function(input, output, session) {
  
  
  # close the R session when Chrome closes
  session$onSessionEnded(function() {
    stopApp()
    #q("no")
  })
  
  observe({
    if (input$goButton) {
      # enable the download button
      shinyjs::enable("DownloadData")
      # change the html of the download button
      shinyjs::html("DownloadData",
                    sprintf("<i class='fa fa-download'></i>
                             <b>Download optimized sequence<b>",
                            round(runif(1, 1, 10000))
                    )
      )
    }
  })
  
  hexFunc <- eventReactive(input$goButton, {
    
    
    res <- ModCon(input$enteredCDS, input$sdPosition, upChangeCodonsIn=input$upchangecodons,
                  downChangeCodonsIn=input$upchangecodons, optimizeContext=input$increaseHZEI,
                  sdMaximalHBS=input$sdMaximalHBS, acMaximalMaxent=input$acMaximalMaxent,
                  optiRate=input$optiRate, nGenerations=input$nGenerations, parentSize=input$parentSize,
                  startParentSize=input$startParentSize,  bestRate=input$bestRate,
                  semiLuckyRate=input$semiLuckyRate, luckyRate=input$luckyRate, 
                  mutationRate=input$mutationRate, nCores=-1)
    
    
    
    return(res)
    
  })
  
  ## Generate sequence
  output$quick_info <- renderText({
    
    ## Access population generation result
    res <- hexFunc()
    
    ## Saving wildtype and optimized sequence
    before <- input$enteredCDS
    after <-  res
    
    ## Define region of interest
    upChangeCodons <- (input$upchangecodons * 3) - 1
    downChangeCodons <- (input$upchangecodons * 3)
    inputSeq <- substr(after, (input$sdPosition - 56), (input$sdPosition+4))
    upSeqAfter <- inputSeq
    inputSeq <- substr(after, (input$sdPosition+5), (input$sdPosition+5+60))
    downSeqAfter <- inputSeq
    
    ## Define region of interest
    upChangeCodons <- (input$upchangecodons * 3) - 1
    downChangeCodons <- (input$upchangecodons * 3)
    inputSeq <- substr(before, (input$sdPosition - 56), (input$sdPosition+4))
    upSeqBefore <- inputSeq
    inputSeq <- substr(before, (input$sdPosition+5), (input$sdPosition+5+60))
    downSeqBefore <- inputSeq
    
    ## Calculate HZEI per subsequence
    upSeqAfterHZEI <- calculateHZEIint(upSeqAfter)
    downSeqAfter   <- calculateHZEIint(downSeqAfter)
    upSeqBefore    <- calculateHZEIint(upSeqBefore)
    downSeqBefore  <- calculateHZEIint(downSeqBefore)
    
    beforeSSHW <- upSeqBefore - downSeqBefore
    afterSSHW  <- upSeqAfterHZEI - downSeqAfter
    
    
    #Return the Difference in the scores
    HTML(paste0("SSHW was changed from ", round(beforeSSHW/nchar(upSeqBefore),1),
                " HZEI/nt to ", round(afterSSHW/nchar(upSeqAfter),1)," HZEI/nt."))
    
    
  })
  
  output$DownloadData <- downloadHandler(
    filename = function() {
      paste('HEXCO-', Sys.Date(), '.txt', sep='')
    },
    content = function(file) {
      
      ## Access population generation result
      res <- hexFunc()
      
      ## Saving wildtype and optimized sequence
      before <- input$enteredCDS
      after <-  res
      
      ## Define region of interest
      upChangeCodons <- (input$upchangecodons * 3) - 1
      downChangeCodons <- (input$upchangecodons * 3)
      inputSeq <- substr(after, (input$sdPosition - 56), (input$sdPosition+4))
      upSeqAfter <- inputSeq
      inputSeq <- substr(after, (input$sdPosition+5), (input$sdPosition+5+60))
      downSeqAfter <- inputSeq
      
      ## Define region of interest
      upChangeCodons <- (input$upchangecodons * 3) - 1
      downChangeCodons <- (input$upchangecodons * 3)
      inputSeq <- substr(before, (input$sdPosition - 56), (input$sdPosition+4))
      upSeqBefore <- inputSeq
      inputSeq <- substr(before, (input$sdPosition+5), (input$sdPosition+5+60))
      downSeqBefore <- inputSeq
      
      ## Calculate HZEI per subsequence
      upSeqAfterHZEI <- calculateHZEIint(upSeqAfter)
      downSeqAfter   <- calculateHZEIint(downSeqAfter)
      upSeqBefore    <- calculateHZEIint(upSeqBefore)
      downSeqBefore  <- calculateHZEIint(downSeqBefore)
      
      beforeSSHW <- upSeqBefore - downSeqBefore
      afterSSHW  <- upSeqAfterHZEI - downSeqAfter
      
      
      
      outputTable <- data.frame("Results:"=
                                  c(paste0("Before HZEI/nt: ",round(beforeSSHW/nchar(upSeqBefore),1)),
                                    paste0("After HZEI/nt: ",round(afterSSHW/nchar(upSeqAfter),1)),
                                    paste0("Output CDS: ", after)))
      row.names(outputTable) <- NULL
      names(outputTable) <- NULL
      write.table(outputTable, file, row.names = F, quote = F)
    }
  )
  
  # disable the downdload button on page load
  shinyjs::disable("DownloadData")
  
}
shinyApp(ui = ui, server = server)
