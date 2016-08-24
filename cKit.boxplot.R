options(stringsAsFactors = F)
load(file = "./gastrulation/data.Rdata")
gene <- "Kit"

library(ggplot2)
library(reshape2)
expr <- t(exprMat[gene,])
for(i in c("embryoStage","cellCategory", "cluster")){
  data <- melt(cbind(expr,annot_col[rownames(expr),i, drop =F]), id.vars = i, measure.vars = gene, variable.names = i)
  ggplot()+ geom_boxplot(data = data, mapping = aes(x = get(i), y = value, fill = get(i))) + 
            scale_fill_manual(
              values = annot_colors[[i]],
             name=i
            ) + 
            theme_bw() + theme(
              plot.background = element_blank()
              ,panel.grid.major = element_blank()
              ,panel.grid.minor = element_blank()
              #,panel.border = element_blank()
              ,axis.text.x = element_text(angle = 45, hjust = 1)
            ) +
            xlab(i) + ylab(paste0("Expression level")) + ggtitle(gene)
    ggsave(filename = paste(gene,i,"boxplot.pdf",sep = "."),width = 6,height = 4)
}