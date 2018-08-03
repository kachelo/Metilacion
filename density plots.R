#Para leer archivos csv, seprated values

island<- read.table(file.choose(), header=TRUE)

#####Generate a new table to plot######
tabla<-matrix(nrow=nrow(island)*ncol(island),ncol=2)
colnames(tabla)<-c("Bvalues","samples" )
tabla<-as.data.frame(tabla)


tabla[,1]<-island[,1]
tabla[1:nrow(island),2]<-"a"

tabla[nrow(island)+1:(2*nrow(island)),1]<-island[,2]
tabla[nrow(island)+1:(2*nrow(island)),2]<-"b"

tabla[(2*nrow(island)+1):(3*nrow(island)),1]<-island[,3]
tabla[(2*nrow(island)+1):(3*nrow(island)),2]<-"c"

tabla[(3*nrow(island)+1):(4*nrow(island)),1]<-island[,4]
tabla[(3*nrow(island)+1):(4*nrow(island)),2]<-"d"

ggplot(tabla, aes(x = Bvalues, fill = samples)) + geom_density(alpha = 0.5)
