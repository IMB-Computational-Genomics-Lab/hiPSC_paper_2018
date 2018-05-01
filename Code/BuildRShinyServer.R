# This script is used to build the Shiny server, available at: http://computationalgenomics.com.au/shiny/hipsc/

#global
library('shiny')
library('shinyjs')
library('plotly')
library('RColorBrewer')
library('data.table')

path='/Users/quan.nguyen/Documents/Powell_group_MacQuan/HiPSC/FullRun/Shiny_Server/Live_Shiny_HiPSC/'
summary_data_big <- function(id, cluster, dat_gene){
  cellCount <- length(cluster)
  pos_idx <- which(dat_gene[cluster]>0)
  posCount <- length(pos_idx)
  posPercent <- posCount/cellCount*100
  meanExprs <- mean(dat_gene[cluster])
  cluster_exprs <- dat_gene[cluster]
  MeanPositive <- mean(cluster_exprs[pos_idx])
  return(list(cellCount,posCount,round(posPercent,3),round(meanExprs,3),round(MeanPositive,3)))
  cluster_exprs=NULL
}

Names <-readRDS(paste0(path,"GeneNames.RDS"))
#server
server <- function(input, output, session) {
  #A function to run everything after the data is loaded
  LetsRock <- function(expression,dat3d_sh){
    cluster1 <- which(dat3d_sh$cluster==1)
    cluster2 <- which(dat3d_sh$cluster==2)
    cluster3 <- which(dat3d_sh$cluster==3)
    cluster4 <- which(dat3d_sh$cluster==4)

    expressionC1 <-expression[cluster1,]
    expressionC2 <-expression[cluster2,]
    expressionC3 <-expression[cluster3,]
    expressionC4 <-expression[cluster4,]
    listExpression <-list(expressionC1,expressionC2,expressionC3,expressionC4)  #list order matters

    dat3dC1 <-dat3d_sh[cluster1,]
    dat3dC2<-dat3d_sh[cluster2,]
    dat3dC3<-dat3d_sh[cluster3,]
    dat3dC4<-dat3d_sh[cluster4,]
    list_dat3d <-list(dat3dC1,dat3dC2,dat3dC3,dat3dC4) #list order matters
    #-----------------
    output$scatter3d <- renderPlotly({

      if (input$name == 'ColourCluster'){
        p <- plot_ly(dat3d_sh, x = ~tSNE1, y = ~tSNE2,  z= ~tSNE3, type="scatter3d", mode = 'markers',
                     marker = list(opacity = 0.7, size=5),
                     color = ~factor(cluster),
                     colors= c('#8A2022', '#CCA47C', '#788E2B', '#2DCBF2'),
                     hoverinfo = 'text',
                     text = ~paste('Cluster:', cluster, '<br> Batch:', batches)) %>%
          layout(title = paste0('Displaying Subpopulation'),
                 xaxis = list(showgrid = FALSE),
                 yaxis = list(showgrid = FALSE),
                 showlegend = T,
                 legend=list(bgcolor='#E2E2E2',bordercolor='#FFFFFF',borderwidth=2))
      } else if (input$name != 'ColourCluster'){
        if (input$cluster == 'All'){
          gene_idx <- which(Names==input$name)
          Log2GeneExpression <- log2(expression[,gene_idx]+1) #in expression matrix, genes are in columns
          sizeExpression <- as.numeric(as.logical(Log2GeneExpression))*10 + 3
          p <- plot_ly(dat3d_sh, x = ~tSNE1, y = ~tSNE2,  z= ~tSNE3, type="scatter3d", mode = 'markers',
                       color = ~Log2GeneExpression,
                       marker = list(opacity = 0.75, size=sizeExpression),
                       hoverinfo = 'text',
                       text = ~paste('Cluster:',
                                     cluster, '<br> Batch:', batches,'<br> Exprs:',
                                     round(expression[,gene_idx], digits=3))) %>%
            layout(title = paste0('Displaying gene ', input$name),
                   showlegend = T, legend=list(bgcolor='#E2E2E2',
                                               bordercolor='#FFFFFF',borderwidth=2))
        }
        else if (input$cluster != 'All'){
          cluster<- as.numeric(input$cluster)
          expression_cluster <- listExpression[[cluster]] #select cells in expression matrix (cells as rows)
          dat3d_cluster <-list_dat3d[[cluster]]
          gene_idx <- which(Names==input$name)
          Log2GeneExpression <- log2(expression_cluster[,gene_idx] + 1)
          sizeExpression <- (as.numeric(as.logical(Log2GeneExpression))*8) + 3
          p <- plot_ly(dat3d_cluster, x = ~tSNE1, y = ~tSNE2,  z= ~tSNE3, type="scatter3d", mode = 'markers',
                       color = ~Log2GeneExpression,
                       marker = list(opacity = 0.75, size=sizeExpression),
                       hoverinfo = 'text',
                       text = ~paste(' Cluster:', cluster, '<br> ', batches,'<br> Exprs:', round(expression_cluster[,gene_idx],digits=3))) %>%
            layout(title = paste0('Displaying gene ', input$name, ', Cluster ', input$cluster),
                   showlegend = T, legend=list(bgcolor='#E2E2E2',bordercolor='#FFFFFF',borderwidth=2))
          p}
      }
    })

    #what it input$name ="ColorClusters"

    density <- reactive({
      if (input$name != 'ColourCluster'){
        gene_idx  <- which(Names==input$name)
        Exprs_gene <-as.data.frame(expression[,gene_idx]) # in the expression matrix genes are columns
        Exprs_gene$Cluster <- dat3d_sh$cluster
        colnames(Exprs_gene) <-c('Gene', 'Cluster')
        p <- ggplot(data=(Exprs_gene[Exprs_gene$Gene>0,]) , aes(log2(Gene+1))) +
          geom_density(aes( y=..scaled.., fill=as.factor(Cluster)), alpha=0.5) +
          theme_bw() + theme(text = element_text(size=14)) +
          ylab('Scaled Density') + xlab('log2(Expression+1)') +
          scale_fill_manual(name="Cluster", values=c('#da5932','#8e5ad5','#6cac30','#d83b9a'),limits=c('1','2','3','4')) +
          ggtitle(paste0("Expression of Positive Cells for Gene: ", input$name))
        m <- list(l = 100, r = 0, b = 100, t = 100, pad = 4) #setting l=100, r=100 is better to make room for the y-axis label
        p <- ggplotly(p)  %>% layout(autosize = F, width = 800, height = 400, margin = m) #setting width =800, not 500 spans well the whole area, the height needs to be 400 not 500
        p
      }
    })
    output$density <- renderPlotly(density())

    ClusterSummary <- reactive({
      if (input$name != 'ColourCluster'){
        gene_idx <- which(Names==input$name)
        DayAll <- 1:nrow(dat3d_sh)
        dat_gene <- expression[,gene_idx]

        all_output <- summary_data_big(gene_idx, DayAll, dat_gene)
        cluster1_output <- summary_data_big(gene_idx, cluster1, dat_gene)
        cluster2_output <- summary_data_big(gene_idx, cluster2, dat_gene)
        cluster3_output <- summary_data_big(gene_idx, cluster3, dat_gene)
        cluster4_output <- summary_data_big(gene_idx, cluster4, dat_gene)

        data.frame(
          ClusterData = c("Cell Count","Positive Cells","Percent Positive Cells","Mean Expression All","Mean Expression Positive"),
          AllClusters = as.character(c(all_output[[1]],all_output[[2]],all_output[[3]],all_output[[4]],all_output[[5]])),
          Cluster1 = as.character(c(cluster1_output[[1]],cluster1_output[[2]],cluster1_output[[3]],cluster1_output[[4]],cluster1_output[[5]])),
          Cluster2 = as.character(c(cluster2_output[[1]],cluster2_output[[2]],cluster2_output[[3]],cluster2_output[[4]],cluster2_output[[5]])),
          Cluster3 = as.character(c(cluster3_output[[1]],cluster3_output[[2]],cluster3_output[[3]],cluster3_output[[4]],cluster3_output[[5]])),
          Cluster4 = as.character(c(cluster4_output[[1]],cluster4_output[[2]],cluster4_output[[3]],cluster4_output[[4]],cluster4_output[[5]]))
        )
      }
    })

    output$SummaryData <- renderTable({
      ClusterSummary()
    })


    output$AboutText <- renderText({ "This site contains user-friendly tools for displaying and summarising data from the single cell sequencing of 18,787 human induced pluripotent stem cells (hiPSC). The detailed manuscript is being under-review in Genome Research and the preprint is accessible from bioRxiv at: \n

      Nguyen, Q., Lukowski, S., Chiu, H., Senabouth, A., Bruxner, T., Christ, A., Palpant, N., and Powell, J. (2017). Single-Cell Transcriptome Sequencing Of 18,787 Human Induced Pluripotent Stem Cells Identifies Differentially Primed Subpopulations. bioRxiv, doi:https://doi.org/10.1101/119255. \n

      Usage instruction: \n

      1. The default display, when the plot is first loaded is a tSNE plot showing all  18,787 cells in four clusters. Each point represents one cell. The mouse-over action will interactively show information for the cluster ID and batch ID of each cell. There are 5 batches (each batch is an independent biological sample of hiPSC cells). An example for exploring this dataset is that: no batch effect is observed when mousing-over the cells in the tSNE display. \n
      2. Users can start exploring the data by selecting a gene in the gene selection panel. The selection action will trigger the calculation and display of: expression value in each cell for the selected gene in the tSNE plot, a summary table, and a density plot . \n
      3. The tSNE plot can be used to display all cells by selecting the option 'All', or to display cells in each cluster by selecting a cluster IDs (1, 2, 3, or 4). \n
      4. The summary table underneath the tSNE diplays the total cells, counts of the cells with positive expression for the selected gene, mean expression values of positive cells for all cells and specific fot each cluster. \n
      5. The density plot displays the distribution of the expression of a gene in each cluster, allowing the visual comparison of gene expression between clusters. \n
      6. Plots and Tables are downloadable, by clicking the 'Download' (for tables), or 'save image' buttons.

      For questions and comments, please contact Dr. Joseph E. Powell at j.powell@imb.uq.edu.au." })
    }


  observeEvent(input$radio,{
    if(input$radio == "checkbox"){
      withProgress(message = 'Reading Large Data! (PLEASE WAIT)', value = 0, {

        for(i in 1:2){
          incProgress(1/2, detail = paste("File #", i, " of 2")) #e.g. 1 of 3
          expression <- readRDS(paste0(path,'Expression_CPM_unLog_minus1_positive_Tp.RDS'))
          dat3d_sh  <- readRDS(paste0(path,'tSNE_3D.RDS'))
          LetsRock(expression=expression,dat3d_sh=dat3d_sh)
        } } ) } else {
          withProgress(message = 'Reading Data!', value = 0, {

            for(i in 1:2){
              incProgress(1/2, detail = paste("File #", i, " of 2"))
              expression <- readRDS(paste0(path,'Expression_CPM_unLog_minus1_positive_Tp_20pc.RDS'))
              dat3d_sh <- readRDS(paste0(path,'tSNE_3D_20pc.RDS'))
              LetsRock(expression=expression,dat3d_sh=dat3d_sh)
            } } ) }
  })

    }

shinyApp(ui, server)

#UI
ui <- navbarPage("HiPSC Single Cell Expression Data",fluid = TRUE,
                 column(8,
                        img(src=paste0(path,"/www/animated.gif"), height = 300, width = 350, align='left')
                 ),
                 tabPanel("DataMining",
                          column(12,
                                 radioButtons("radio","Data Size:",
                                              choices = list("Small 20% Data Set" = "checkboxsmall",
                                                             "Full Data Set" = "checkbox"),
                                              selected = "checkboxsmall"),
                                 #tSNE
                                 wellPanel(
                                   plotlyOutput("scatter3d", height = 600, width = 800),
                                   selectizeInput("cluster","Select the Cluster to view:",
                                                  choices = c('All',1,2,3,4),
                                                  multiple = T,
                                                  options = list(maxItems = 1, placeholder = 'Select a cluster'),
                                                  selected = 'All')
                                   #summary table
                                 ),
                                 conditionalPanel(condition = "input.name != 'ColourCluster'",
                                                  wellPanel(
                                                    tableOutput("SummaryData"),
                                                    downloadButton('downloadData', 'Download')
                                                  )
                                 ),
                                 #density
                                 fluidRow(
                                   conditionalPanel(condition = "input.name != 'ColourCluster'",
                                                    column(12,
                                                           wellPanel(
                                                             plotlyOutput("density")
                                                           )
                                                    )
                                   )
                                 ),
                                 #gene selection panelsidebarPanel(
                                 sidebarPanel(
                                   column(8,
                                          h3(""),
                                          img(src="animated.png", height = 300, width = 350, align='left')
                                   ),

                                   absolutePanel(
                                     top = 71, right = 20, width = 300,
                                     draggable = TRUE,
                                     fixed = TRUE,
                                     wellPanel(
                                       selectizeInput("name",
                                                      label="Select a gene to view:",
                                                      choices = c('ColourCluster',sort(unique(Names))),
                                                      multiple = T,
                                                      options = list(maxItems = 1, placeholder = 'Select a name'),
                                                      selected = "ColourCluster"),
                                       helpText("Select a gene to explore expression in each cell and supopulation (cluster) by tSNE and density plots and a summary table. To view tSNE for cells in each cluster, use the second selection panel underneath the tSNE plot.")
                                     ),
                                     style = "opacity: 0.92"
                                   )
                                 )
                          )
                 ) ,


                 #About tab ======================================
                 tabPanel("About and Instruction",
                          column(12,
                                 wellPanel(
                                   verbatimTextOutput("AboutText") #Need to use verbatimTextOutput so that is recognises the \n
                                 )
                          )
                 )
)



