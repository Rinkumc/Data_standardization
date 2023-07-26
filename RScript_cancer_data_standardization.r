a=read.table("C:/Users/HP/OneDrive/Documents/GSE200146_counts.Cell.Lines.txt/GSE200146_counts.Cell.Lines.txt",header = T, row.names = 1,sep="\t")
head(a)
a=a[,-1]

sizefactors=colSums(a)
columnSums=apply(a,2,sum)
cpm=apply(a,2,function(x)(x/sum(x))*100000)
cpm["ENSG00000188976","PDC_4"]
cpm
cpm1 = cpm
cpm = cpm[1:100, ]

cpm=log2(cpm+1)

temp = cpm

calculate_zscore = function(temp){
  z_scores=temp
  men=rowMeans(temp)
  std=apply(temp,1,sd)
  for(i in 1:nrow(temp)){
    z_scores[i, 1:ncol(z_scores)] = (temp[i, 1:ncol(temp)] - men[i])/std[i]
  }
  return(z_scores)
}
temp1=calculate_zscore(temp)
newmat=temp1[!rowSums(is.na(temp1)),]
library(ComplexHeatmap)
library(circlize)
Heatmap(newmat,col = colorRamp2(c(-2,0,2),c("orange","white","purple")))
#saveRDS(cpm,"calculate_zscore.rds")

variance=apply(cpm1,1,var)
sorted=sort(variance, decreasing=TRUE)
top100=sorted[1:100]
top100
topnames=names(top100)
topnames
newcpm=cpm1[topnames,]
newtemp = newcpm
temp2=calculate_zscore(newtemp)
Heatmap(temp2,col = colorRamp2(c(-2,0,2),c("orange","white","purple")))

data=temp2
head(data)

metadata = read.csv("C:/Users/HP/OneDrive/Documents/meta.csv")
head(df)
#Set annotation
annotations=HeatmapAnnotation(cell_line=metadata$cell_line,age=metadata$Age,gender=metadata$Gender)
Heatmap(temp2,col = colorRamp2(c(-2,0,2),c("orange","white","purple")),top_annotation=annotations,width = 80, height = 45)
