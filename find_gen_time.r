setwd('C:/Users/topalw/Desktop/pedigree/')
birds <- read.csv('BarnOwls_Legacy_20220420154011_Bird.csv')
clutches <- read.csv('BarnOwls_Legacy_20220420154011_Clutch.csv')
head(birds)
#as.numeric(as.Date(clutches$LayingDate[1]) - as.Date('07/02/2016'))
mom.age <- rep(NA,nrow(clutches))
for(i in 1:nrow(clutches)){
  mom.ring <- clutches$FemaleRing[i]
  mom.birth <- strsplit(birds$HatchDate[match(mom.ring,birds$RingId)],split=' ')[[1]][1]
  mom.age[i] <- as.numeric(as.Date(clutches$LayingDate[i]) - as.Date(mom.birth))
}
mom.age <- mom.age[mom.age > 0]
mom.age <- mom.age[! is.na(mom.age) ]
range(mom.age,na.rm = T)
length(mom.age[!is.na(mom.age)])
