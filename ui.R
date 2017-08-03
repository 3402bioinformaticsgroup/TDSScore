# ui.R

library(shiny)
#library(edgeR)

shinyUI(fluidPage(
 titlePanel("Transcriptional Disease Signature analysis Pt3"),
 sidebarLayout(
   sidebarPanel(

      fileInput("TDS_reference", label = h3("TDS.reference(.csv)")),
      fileInput("Drug_data", label = h3("Drug.expression(.csv)")),
      textInput("LogFC.cutoff", "LogFC cutoff value:", "1"),
      textInput("outputfile", "Output file name","Drug.score.result"),
      downloadButton('downloadData', 'Download')

   ),

 mainPanel(
#h2(textOutput("caption", container = span)),
   tableOutput("tds.reference.view"),
   tableOutput("drug.data.view"),

fluidRow(
column(6, plotOutput("plotZScore")),
column(6, plotOutput("plotDrugScore"))
)

#tableOutput("view")

              )
           )
        )

     )


