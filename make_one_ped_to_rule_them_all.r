#0.0 wd / packages / files 
birds <- read.csv('BarnOwls_Legacy_20220420154011_Bird.csv')
clutch <- read.csv('BarnOwls_Legacy_20220420154011_Clutch.csv')
metadata <- read.csv('refpanel_metadata.csv')
library(kinship2)
library(pedigree)

#0.1 - make empty dataset
ids <- data.frame('id'=unique(birds$RingId),'dadid'=rep('',length(unique(birds$RingId))),
                  'momid'=rep('',length(unique(birds$RingId))),
                  'famid'=rep(0,length(unique(birds$RingId))),
                  'sex'=rep(0,length(unique(birds$RingId))))
ids <- ids[ids$id != '' | ! is.na(ids$id),]

#1.0 - FAMILIES
ids$famid <- birds$BornClutchId[match(ids$id,birds$RingId)] # get clutch id
ids$momid <- clutch$FemaleRing[match(ids$famid,clutch$Id)] # get momid
ids$dadid <- clutch$MaleRing[match(ids$famid,clutch$Id)] # get dadid

#2.0 - SEX 
#2.1 - make merged sex in birds 
birds$sex <- rep('',nrow(birds))
for(i in 1:nrow(birds)){
  if(birds$GeneticSex[i] != ''){
    birds$sex[i] <- birds$GeneticSex[i]
  } else {birds$sex[i] <- birds$PhenotypeSex[i]}
}
#2.2 - put sex in df and digitize 
ids$sex <- birds$sex[match(ids$id,birds$RingId)]
ids$sex[ids$sex==''] <- 3
ids$sex[ids$sex=='Male'] <- 1
ids$sex[ids$sex=='Female'] <- 2

#2.3 
ids$sex[ids$id %in% ids$dadid] <- 1
ids$sex[ids$id %in% ids$momid] <- 2

#2.4 - use WGS corrected sex for some individuals 
for(id in metadata$NEWname){
  ids$sex[ids$id==id] <- unique(metadata$sex[metadata$NEWname == id])
}

#2.5 - check if parents have correct sex WHAT ABOUT THESE GUYS?
 hefemales <- which(ids$sex == 1 & ids$id %in% ids$momid) # male mothers - 7041/7201 
 ids[hefemales,] 
 ids[ids$momid %in% ids$id[hefemales],] # one pair mislabeled either in sex or parental role + 1 hefemales
 ids[hefemales,5] <- c(2,2) # FIX
 shemales<-  which(ids$sex == 2 & ids$id %in% ids$dadid) # female fathers - 7794
 ids[shemales,] 
 ids[shemales,5] <- c(1)
 

#3.0 - make ped
 
#3.1 - make corrections 
ids$dadid[ids$dadid == ''] <- NA
ids$momid[ids$momid == ''] <- NA
ids <- ids[! ids$id %in% c('',NA,' '),] # remove empty shit
ids$sex <- as.numeric(ids$sex)
# check for missing individuals that exist in mother/father col
ids[! ids$momid %in% ids$id & ! is.na(ids$momid),]
ids[! ids$dadid %in% ids$id & ! is.na(ids$dadid),] # lel
ids <- rbind(ids,data.frame('id'='M031195','dadid'=NA,'momid'=NA,
                            'famid'=NA,'sex'=1))
#3.2 - add single parents function
# Structure is [,1] = IndId; [,2] = DadId; [,3] = MomId; [,4] = famid; [,5] = sex
fix.single.parents <- function(pedigree){
  for(i in 1:nrow(pedigree)){
    # make fake dad
    if(is.na(pedigree[i,2]) & !is.na(pedigree[i,3])){
      # make a fake parent ID (_1 for mother, _2 for father)
      fake_parent <- paste(pedigree[i,3],'1',sep='_')
      # make a dad line (assume founder)
      fake_df <- data.frame('id' = fake_parent,
                            'dadid' = NA, #dadid
                            'momid' = NA, #momid
                            'famid' = NA, # clutchId we do not care about this any more
                            'sex' = 1) #sex 
      colnames(fake_df) <- colnames(pedigree) # sanity check on names
      # replace id for dad for all kids of the clutch to make them a nice family
      if(!is.na(pedigree[i,4])){
        pedigree[which(pedigree[,4] == pedigree[i,4]),2] <- fake_parent
      } else{
        pedigree[i,2] <- fake_parent
      }
      # add mother line 
      pedigree <- rbind(pedigree , fake_df)
      # make fake mom
    } else if(!is.na(pedigree[i,2])& is.na(pedigree[i,3])){
      fake_parent <- paste(pedigree[i,2],'2',sep='_')
      fake_df <- data.frame('id' = fake_parent,
                            'dadid' = NA, #dadid
                            'momid' = NA, #momid
                            'famid' = NA, # clutchId we do not care about this any more
                            'sex' = 2) #sex 
      colnames(fake_df) <- colnames(pedigree)
      if(!is.na(pedigree[i,4])){
        pedigree[which(pedigree[,4] == pedigree[i,4]),3] <- fake_parent
      } else{
        pedigree[i,3] <- fake_parent
      }
      pedigree <- rbind(pedigree , fake_df)
    } 
  }
  return(pedigree[!duplicated(pedigree[,1]),])
}

ids <- fix.single.parents(ids)

#3.3 - find sequenced 
sequenced <- unique(metadata$NEWname[metadata$pop=='CH'])
ids$affected <- rep(0,nrow(ids))
ids$affected[match(sequenced,ids$id)] <- 1    

#3.3 - make actual ped
ids$ped_famid <- makefamid(ids$id, ids$dadid,ids$momid)
ped <- pedigree(id=ids$id,dadid= ids$dadid,momid =  ids$momid, sex=ids$sex , famid=ids$ped_famid,
                affected = ids$affected)
seqs <- unique(ids$ped_famid[ids$affected==1])
for(i in unique(ids$ped_famid)){print(ped[i+1])}
fam <- 206
fam <- 187
fam <- 88
plot.pedigree(ped[fam+1],affected=ids$affected[ids$ped_famid==fam],cex=0.6,symbolsize = 0.6,
              align=T)
#save csv ped
#write.csv(ids,'total_pedigree_052022.csv',quote=F,row.names = F)

#get and save kinship matrix for sequenced swiss
swiss <- match(unique(metadata$NEWname[metadata$pop=='CH']), ids$id)
k1 <- kinship(ids$id,ids$dadid,ids$momid)[swiss,swiss]
#write.table(k1,'ch_kinship_from_ped.table')

tmp <- pedigree.shrink(ped[3],avail=ids$affected[ids$ped_famid==3],
                       affected=ids$affected[ids$ped_famid==3])
pdf('trimmed_ped_fam3.pdf',width=16,height=12)
plot.pedigree(tmp$pedObj, cex=0.6 ,symbolsize = 0.6)
dev.off()
