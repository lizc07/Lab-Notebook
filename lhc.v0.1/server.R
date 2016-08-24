# server.R
library(shiny)
library(ggplot2)
library(reshape2)

load(file = "data.Rdata")

shinyServer(function(input, output) {
  sel_gene <- reactive({input$sel_gene})
  sel_mirna <- reactive({input$sel_mirna})
  analysis <- reactive({ input$analysis})
  category <- reactive({ input$category})
  heat.height <- reactive({ input$heat.height})
  
  p <- ggplot() + theme_grey() + theme(panel.background=element_rect(fill='transparent', color='gray'),
                                       legend.key=element_rect(fill='transparent', color='transparent'),
                                       axis.text=element_text(color='black'), legend.position = "none")
  
  

#   if(is.null(sel_gene()) && is.null(sel_mirna())){
#     return()
#   }
  
  
  output$plot <- renderPlot({
    category <- category()
    if(category == "N_C_TYPES")
      subData <- subset(brcaData,subset = !is.na(get(category)) , c(category, sel_gene(), sel_mirna()))
    else
      subData <- subset(brcaData,subset = !is.na(get(category)) & N_C_TYPES =="Tumor" , c(category, sel_gene(), sel_mirna()))
    if(analysis() == "Differential expression"){
      ggData <- melt(subData)
      
      q <- p + geom_boxplot(mapping = aes_string(x = category,y = "value", fill = category),data = ggData) + facet_grid(facets = . ~ variable) +xlab(category) + ylab("log2 expression")
      #t.test(value ~ Sample,ggData) #p-value = 0.4909
      print(q)
    }else if(analysis() == "Correlation"){
      ggData <- melt(subData, id.vars = c(sel_mirna(), category))
      print("hehe")
      q <- p + geom_point(mapping = aes_string(x =ggData[,sel_mirna()], y = "value"),color = "darkcyan", data = ggData) + facet_grid( as.formula(paste0(category,"~variable"))) +xlab(sel_mirna())
      print(q)
    }
  })
  
  output$plot.ui <- renderUI({
    plotOutput("plot",height = heat.height())
  })
  
  output$testTab <- renderTable({
    category <- category()
    if(category == "N_C_TYPES")
      subData <- subset(brcaData,subset = !is.na(get(category)) , c(category, sel_gene(), sel_mirna()))
    else
      subData <- subset(brcaData,subset = !is.na(get(category)) & N_C_TYPES =="Tumor" , c(category, sel_gene(), sel_mirna()))

    if(analysis() == "Differential expression"){
      diffTab <- matrix(NA, 2*length(unique(subData[,1])),length(unique(subData[,1])))
      colnames(diffTab) <- unique(subData[,1])
      rownames(diffTab) <- paste(rep(colnames(subData)[2:3],each = length(unique(subData[,1]))),rep(unique(subData[,1]),2),sep = ".")
      for(i in colnames(diffTab) ){
        for(j in rownames(diffTab)){
          tmp <- strsplit(j, split = "\\.")[[1]]
          if(tmp[2] != i){
            testRes <-t.test(subData[subData[,category] == tmp[2], tmp[1]],subData[subData[,category] == i, tmp[1]])
            diffTab[j,i] <- testRes$p.value
          }
        }
      }
      
      return(diffTab)
    }else if(analysis() == "Correlation"){
      corrTab <- matrix(NA,length(unique(subData[,1])),1)
      colnames(corrTab) <- colnames(subData)[2]
      rownames(corrTab) <- unique(subData[,1])
      for(gene in colnames(corrTab) ){
        for(group in rownames(corrTab)){
          subData.cor <- na.omit(subData[subData[,category] == group,c(gene,sel_mirna())])
          testRes <- cor.test(subData.cor[,1], subData.cor[,2],method = "spearman")
          corrTab[group,gene] <- paste(" r = ",testRes$estimate, ", p = ", testRes$p.value,sep = "")
        }
      }
      return(corrTab)
    }
  })
})