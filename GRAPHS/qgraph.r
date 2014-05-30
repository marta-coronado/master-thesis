con objeto kk intentar instalar el qgraph library(qgraph)
qgraph(cor(kk),vsize=2,minimum=0.2,filetype="png")

http://s122.codeinspot.com/q/1622668

> library(ggplot2)
> p <- ggplot(data=melt(kk), aes(x=X1, y=X2, color=value))
> p + geom_point(size=5, alpha=0.7) + scale_color_gradient2()

> qgraph(cor(kk),label.font=14,color="gray90",border.width=2.5,border.color="gray40",label.color="gray5",filetype="png")
> qgraph(cor(kk),label.font=14,color="gray90",border.width=2.5,border.color="gray40",label.color="gray5",directed=TRUE,arrows=TRUE,layout="circular",asize=2,filename="kk", filetype="png")


> qgraph(cor(kk),label.font=6,color="gray90",border.width=2.5,border.color="gray40",label.color="gray5",nodeNames = names,labels=c("trn","exon","size","dist","brd","expr","AAs","mut","rec"),title="Correlations",details=TRUE,filetype = "sgv")

> qgraph(cor(kk),label.font=6,color="gray90",border.width=2.5,border.color="gray40",label.color="gray5",nodeNames = names,labels=c("trn","exon","size","dist","brd","expr","AAs","mut","rec"),title="Correlations",details=TRUE,borders=FALSE)

qgraph(cor(kk),label.font=6,color="gray90",border.width=2.5,border.color="gray40",label.color="gray5",labels=c("  transcripts  ","      exon      ","      size      ","   distance   ","    breadth    "," expression ","substitutions","   mutation   ","recombination"),title="Correlations",details=TRUE,borders=FALSE,shape="ellipse",vsize=12,vsize2=5,label.prop=0.95)
qgraph(cor(kk),label.font=6,color="gray90",border.width=2.5,border.color="gray40",label.color="gray5",labels=c("  transcripts  ","      exon      ","      size      ","   distance   ","    breadth    "," expression ","recombination"),title="Correlations",details=TRUE,borders=FALSE,shape="ellipse",vsize=12,vsize2=5,label.prop=0.95)


qgraph(cor(kk),label.font=6,color="gray90",border.width=2.5,border.color="gray40",label.color="gray5",labels=c("  transcripts  ","inclusion levels", "      size     ", "   distance   ", "substitutions", "  mutation  ", "     order     ","recombination"),title="Correlations",details=TRUE,borders=FALSE,shape="ellipse",vsize=12,vsize2=5,label.prop=0.95)
qgraph(cor(kk),label.font=6,color="gray90",border.width=2.5,border.color="gray40",label.color="gray5",labels=c("  transcripts  ","inclusion levels", "      size     ", "   distance   ", "     order     ","recombination"),title="Correlations",details=TRUE,borders=FALSE,shape="ellipse",vsize=12,vsize2=5,label.prop=0.95)
