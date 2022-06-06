setwd('C:/Users/topalw/Desktop/pedigree/')
d <- read.csv('succint_ped.csv')
library(kinship2)
d1 <- d[d$seq==1,]
dadf <- unique(d1$dadid[which(! d1$dadid %in% d1$id & ! is.na(d1$dadid))])
momf <- unique(d1$momid[which(! d1$momid %in% d1$id & ! is.na(d1$momid))])
dadf.df <- add.founder(dadf,1)
momf.df <- add.founder(momf,2)
founders <- rbind(dadf.df,momf.df)
founders$seq <- rep(0,nrow(founders))
d1 <- rbind(d1,founders)
famid <- makefamid(d1$id,d1$dadid,d1$momid)
ped <- pedigree(d1$id,d1$dadid,d1$momid,d1$sex,d1$seq,famid=famid)
pdf('succint_ped_per_family.pdf',height=12,width=16)
for(i in c(1:4,6:42)){
  plot(ped[i],cex=0.6)
  }
dev.off()
