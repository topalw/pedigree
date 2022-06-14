setwd('C:/Users/topalw/Desktop/PhD/Analyses/5.ref_panel/pedigree/')
lc <- read.csv('Metadata 3Kowls - 3kOwls_LowCoverage.csv', h=T)
ped1 <- read.csv('genetically_fixed_pedigree.csv',h=T)

ped1$seq_lc <- rep(0,nrow(ped1))
ped1$seq_lc[match(lc$New_name, ped1$id)] <- 1
# duplicates 
dups <- lc$New_name[grep('a',lc$New_name)]
dups2 <- rep('',length(dups))
for(i in 1:length(dups)){
  dups2[i] <- strsplit(dups[i],split='a')[[1]][1]
}
ped1$seq_lc[match(dups2, ped1$id)] <- 1
colnames(ped1)[5] <- 'seq_hc'
# check
sum(ped1$seq_lc)

# sex
d <- read.table('CoverageLowCov_Samples.txt',h=T)
d$sample <- lc$New_name[match(d$LibName,lc$Library_name)]
head(d)
ped1$lc_g_sex <- rep('',nrow(ped1))
ped1$lc_g_sex[match(d$sample,ped1$id)[!is.na(match(d$sample,ped1$id))]] <- 
  d$Sex[!is.na(match(d$sample,ped1$id))]
# match name for dups 
tmp <- d$sample[grep('a',d$sample)]
dups.d <- cbind(dups,dups2)
tmp2 <- ped1$id[ped1$seq_lc==1 & ped1$lc_g_sex=='']
tmpd <- match(dups.d[match(tmp2,dups.d[,2]),1],d$sample)
ped1$lc_g_sex[match(tmp2[-2],ped1$id)]<- d$Sex[tmpd[!is.na(tmpd)]]
pdf('sexxx_diff_3k.pdf',height=12,width=12)
plot(as.numeric(ped1$lc_g_sex[ped1$seq_lc==1 & ped1$id !='889830']) - as.numeric(ped1$sex[ped1$seq_lc==1& ped1$id !='889830']),
     pch=20,cex=3,
     col=c('blue','orange','grey')[as.factor(ped1$sex[ped1$seq_lc ==1& ped1$id !='889830'])],
     ylab='genetic - pedigree sex',yaxt='n')
axis(2,at = c(-2,-1,1),labels = c('p.u->g.m','u->f / f->m',
                                    'p.m->g.f'),cex=3,padj=0)
legend('bottomright',pch = rep(20,3),col=c('blue','orange','grey'),
       legend=c('ped m','ped f','ped u'),cex=3)
 dev.off()

