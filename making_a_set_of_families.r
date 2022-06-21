setwd('C:/Users/topalw/Desktop/PhD/Analyses/5.ref_panel/pedigree/')
d <- read.csv('succint_ped.csv')
library(kinship2)
d1 <- d[d$seq==1,]
dadf <- unique(d1$dadid[which(! d1$dadid %in% d1$id & ! is.na(d1$dadid))])
momf <- unique(d1$momid[which(! d1$momid %in% d1$id & ! is.na(d1$momid))])
add.founder <- function(id,sex=1){
  return(data.frame('id'=id,'dadid'=rep(NA,length(id)),'momid'=rep(NA,length(id)),'sex'=rep(sex,length(id))))
}
dadf.df <- add.founder(dadf,1)
momf.df <- add.founder(momf,2)
founders <- rbind(dadf.df,momf.df)
founders$seq <- rep(0,nrow(founders))
d1 <- rbind(d1,founders)
d1$famid <- makefamid(d1$id,d1$dadid,d1$momid)
h <- read.delim('hubs.list',h=F)
d1$alive <- 0
d1$alive[d1$id %in% h$V2] <- 1
ped <- pedigree(d1$id,d1$dadid,d1$momid,d1$sex,d1$seq,d1$alive,famid=d1$famid)

pdf('succint_ped_per_family_withHUB.pdf',height=12,width=16)

 for(i in c(1:4,6:42)){
  plot.pedigree(ped[i],cex=0.8)
  }
dev.off()
