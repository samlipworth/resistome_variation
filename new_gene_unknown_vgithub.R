
######load libraries#######
library(tidyverse)
library(igraph)
library(ggnewscale)
library(patchwork)

#########Supplementary dataset#################


guuids<-read_csv('./data/guuids',col_names = "Isolate")
guuids$Study<-ifelse(grepl('-',guuids$Isolate),'Oxford','Not Oxford')
guuids$Study<-ifelse(grepl('ERR403',guuids$Isolate),'Norway',guuids$Study)
guuids$Study<-ifelse(grepl('SRR796',guuids$Isolate),'Australia',guuids$Study)
guuids$Study<-ifelse(grepl('ERR439',guuids$Isolate),'BSAC',guuids$Study)
guuids$Study<-ifelse(grepl('ERR434',guuids$Isolate),'BSAC',guuids$Study)
guuids$Study<-ifelse(grepl('ERR121',guuids$Isolate),'Thailand',guuids$Study)
guuids$Study<-ifelse(grepl('ERR306',guuids$Isolate),'Sweden',guuids$Study)
guuids$Study<-ifelse(grepl('ERR435',guuids$Isolate),'BSAC',guuids$Study)
guuids$Study<-ifelse(grepl('ERR485',guuids$Isolate),'BSAC',guuids$Study)
guuids$Study<-ifelse(grepl('ERS199',guuids$Isolate),'Sweden',guuids$Study)
guuids$Study<-ifelse(grepl('ERR460',guuids$Isolate),'BSAC',guuids$Study)

guuids$Accession<-ifelse(guuids$Study=='Norway','PRJEB32059','unknown')
guuids$Accession<-ifelse(guuids$Study=='BSAC','PRJEB4681',guuids$Accession)
guuids$Accession<-ifelse(guuids$Study=='Thailand','PRJEB11403',guuids$Accession)
guuids$Accession<-ifelse(guuids$Study=='Sweden','PRJEB23294',guuids$Accession)
guuids$Accession<-ifelse(guuids$Study=='Oxford','PRJNA604975',guuids$Accession)


t<-filter(guuids,Study =='Not Oxford')
table(guuids$Study)
table(guuids$Accession)

res<-read_tsv('./data/global_pheno.tsv')
res<-select(res,guuid,Coamox,Ceftriaxone,Ciprofloxacin,Gentamicin,Trimethoprim,Piptaz,Ampicillin,Fosfomycin)
guuids<-left_join(guuids,res,by=c("Isolate"="guuid"))


guuids<-janitor::remove_empty(guuids[,4:11],which = "rows")
#8586 E. coli isolates with linked WGS and pheno

#####load data#####
#filter for genes in the classes we are interested in
gene_names<-read_tsv('./data/amrfinder_new.tsv') 
gene_names<-filter(gene_names,`Element type` =='AMR')
gene_names<-filter(gene_names,grepl('BETA-LACTAM',Class) | grepl('GENTAMICIN',Subclass)| grepl('QUINOLONE',Class)|grepl('FOSFOMYCIN',Class)|grepl('TRIMETHOPRIM',Class))
gene_names_keep<-select(gene_names,gene=`Gene symbol`,Class)
gene_names_keep$gene<-str_replace_all(gene_names_keep$gene,'_.*','')
gene_names$`Gene symbol`<-str_replace_all(gene_names$`Gene symbol`,'_.*','')
gene_names<-select(gene_names,gene=`Gene symbol`) %>% distinct()

###the files required for this section are too large to upload to github
#files <- list.files(path="./", pattern="*.dist", full.names=TRUE, recursive=FALSE)

################ only do this if really needed
#out=NULL
#out_family=NULL
#for(file in files){
#   print(file)
#   dists<-read_tsv(file,col_names = F)
#   dists$X1<-paste0(dists$X1,'_',dists$X2)
#   dists$X3<-paste0(dists$X3,'_',dists$X4)
#   guuids_store<-unique(unlist(dists[c(1,3)]))
#   dists$mash_dist<-as.numeric(dists$X7)/as.numeric(dists$X8)
#   dists_family<-filter(dists,mash_dist >=0.7)
#   dists<-filter(dists,mash_dist ==1)
#   dists<-select(dists,X1,X3,mash_dist)
#   names(dists)<-c("X1","X2","weight")
#   dists<-filter(dists,X1!=X2)
#   
#   g<-graph_from_data_frame(dists,directed = F,vertices = guuids_store)
#   
#   g<-components(g)
#   g<-g$membership
#   g<-data.frame(g)
#   g$guuid<-rownames(g)
#   g$guuid<-str_replace_all(g$guuid,'[.].*','')
#   
#   names(g)<-c("allele","guuid")
#   
#   out<-rbind(out,data.frame(file,g))
#   
#   dists_family<-select(dists_family,X1,X3,mash_dist)
#   names(dists_family)<-c("X1","X2","weight")
#   dists_family<-filter(dists_family,X1!=X2)
#   
#   g<-graph_from_data_frame(dists_family,directed = F,vertices = guuids_store)
#   
#   g<-components(g)
#   g<-g$membership
#   g<-data.frame(g)
#   g$guuid<-rownames(g)
#   g$guuid<-str_replace_all(g$guuid,'[.].*','')
#   
#   names(g)<-c("allele","guuid")
#   
#   out_family<-rbind(out_family,data.frame(file,g))
# }
# 
# out$file<-str_replace_all(out$file,".//","")
# out$file<-str_replace_all(out$file,"amrfinder_","")
# out$file<-str_replace_all(out$file,".dist","")
# 
# names(out)<-c("class","allele","guuid")
# out<-separate(out,guuid,into=c('guuid','gene'),sep = '_')
# 
# out_family$file<-str_replace_all(out_family$file,".//","")
# out_family$file<-str_replace_all(out_family$file,"amrfinder_","")
# out_family$file<-str_replace_all(out_family$file,".dist","")
# 
# names(out_family)<-c("class","allele","guuid")
# out_family<-separate(out_family,guuid,into=c('guuid','gene'),sep = '_')
# 
# write_tsv(out,"out.tsv")
# write_tsv(out_family,"out_family.tsv")

#out.tsv is simply the output from the above

out<-read_tsv('./data/out.tsv') %>% filter(!gene=='GlpT')
out_family<-read_tsv('./data/out_family.tsv')

length(unique(out$gene))
out2<-out

guuids<-data.frame(unique(out$guuid))
names(guuids)<-c('guuid')

amrfinder<-read_tsv('./data/amrfinder_gladstone_merger.tsv')
amrfinder<-filter(amrfinder,grepl('AMINOGLYCOSIDE',Class) | grepl('AMINOGLYCOSIDE',Subclass) |
                    grepl('QUINOLONE',Class) | grepl('QUINOLONE',Subclass) |
                    grepl('BETA-LACTAM',Class) | grepl('BETA-LACTAM',Subclass) |
                    grepl('FOSFOMYCIN', Class) | grepl('FOSFOMYCIN',Subclass) |
                    grepl('TRIMETHOPRIM',Class) | grepl('TRIMETHOPRIM',Subclass))
point<-filter(amrfinder,`Element subtype` =='POINT')

amrfinder<-filter(amrfinder,guuid %in% guuids$guuid)

length(unique(out$guuid))
#7153 isolates with at least one AMRFinder hit

core<-c("AcrR","GyrA","ParA","ParC","ParE","PtsI","SoxR","SoxS","UhpT","MarR","GlpT","AmpC","emrR","ftsI")
out_family$allele<-ifelse(out_family$gene %in% core,paste0(out_family$gene,'_1'),out_family$allele)
out_family$amr_associated_gene_family<-paste0(out_family$class,'_',out_family$allele)

length(unique(out_family$amr_associated_gene_family))
#161 unique AMR-associated gene families

length(unique(amrfinder$`Gene symbol`))
#202 ARGs



out$new_allele<-paste0(out$class,'_',out$allele)
length(unique(out$new_allele))

#1292 new alleles

length(unique(out$guuid))

#8883 isolates
out_save<-out



ag_ecol<-out_save
count_table<-ag_ecol %>% distinct(guuid,new_allele,.keep_all = T) %>%  group_by(new_allele) %>% count() 
count_table$Singleton<-ifelse(count_table$n<=1,1,0)
count_table$Double<-ifelse(count_table$n >=2 & count_table$n <10,1,0)
count_table$Tenner<-ifelse(count_table$n >=10,1,0)

##############what prop are already known
amr<-read_tsv('./data/amrfinder.tsv') %>% distinct(`Contig id`,.keep_all = T)
out_single<-distinct(out_save,new_allele,.keep_all = T)
out_single$gene<-ifelse(out_single$gene=="FosA7","FosA7.5",out_single$gene)
out_single$guuid_full<-paste0(out_single$guuid,'_',out_single$gene)
out_single<-left_join(out_single,amr,by=c("guuid_full"="Contig id"))
out_single$exact<-ifelse(out_single$`% Identity to reference sequence`==100,"exact","not exact")
out_single<-select(out_single,new_allele,exact)
out_single<-left_join(out_single,count_table,by=c("new_allele"))
table(out_single$exact,out_single$Singleton)
253+507
253/760

197+335
197/532
chisq.test(out_single$exact,out_single$Singleton)

table(out_single$exact,out_single$Double)
chisq.test(out_single$exact,out_single$Double)
table(out_single$exact,out_single$Tenner)
chisq.test(out_single$exact,out_single$Tenner)
out_single$which<-ifelse(out_single$Singleton==1,'Single',ifelse(out_single$Double==1,'Double',ifelse(out_single$Tenner==1,"Tenner","Other")))
table(out_single$which,out_single$exact)
chisq.test(out_single$which,out_single$exact)
table(out_single$which)

table(count_table$Tenner)
158/1292
table(count_table$Double)
374/1292
table(count_table$Singleton)
#760 are singletons 
760/1292

plot_table<-count_table
plot_table<-separate(plot_table,new_allele,into=c('class','allele'),sep = '_')

plot_table$class<-ifelse(plot_table$class == 'betalactam','beta-lactam',plot_table$class)
Fig1<-ggplot(plot_table) +
  aes(x=n) +
  geom_histogram(bins = 100) + theme_minimal() + 
  xlab("Number of times ARG observed in dataset") + ylab("Count") + scale_x_log10() +
  facet_wrap(~class)

count_table<-select(count_table,-n)

ag_ecol<-left_join(ag_ecol,count_table,by=c("new_allele"="new_allele"))


guuids<-unique(ag_ecol$guuid)
#dates<-read_tsv('~/gn/data/ecoli/Escherichia_guuids_match') %>% select(guuid,date)
#guuids<-data.frame(guuids)
#guuids<-left_join(guuids,dates,by=c("guuids"="guuid"))
#guuids<-arrange(guuids,date)
#guuids<-filter(guuids,!is.na(date))
guuids<-sample(guuids)
collection=NULL
out=NULL

for(g in guuids){
  print(g)
  p<-filter(ag_ecol,guuid==g)
  collection<-rbind(collection,data.frame(p))
  alleles<-length(unique(collection$gene))
  new_alleles<-length(unique(collection$new_allele))
  singletons<-filter(collection,Singleton==1)
  singletons<-length(unique(singletons$new_allele))
  confirmed<-filter(collection,Double==1)
  confirmed<-length(unique(confirmed$new_allele))
  tenner<-filter(collection,Tenner==1)
  tenner<-length(unique(tenner$new_allele))
  out<-rbind(out,data.frame(g,alleles,singletons,confirmed,tenner))
}

out$isolate<-1:nrow(out)
out<-select(out,isolate,g,everything())

out<-pivot_longer(out,names_to = "catalogue",values_to = "n",cols = 3:6)

out$catalogue<-case_when(out$catalogue == 'alleles' ~ 'AMR Finder',
                         out$catalogue == 'confirmed' ~ 'Allele observed at least twice',
                         out$catalogue == 'singletons' ~ 'Allele observed at least once',
                         out$catalogue == 'tenner' ~ 'Allele observed at least ten times')



ecoli_working<-ag_ecol %>% select(guuid,new_allele)
names(ecoli_working)<-c('guuid','gene')

ecoli_working<-data.frame(table(ecoli_working)) %>% spread(gene,Freq)



ecoli_working<-select(ecoli_working,-guuid)

library(micropan)
double<-filter(ag_ecol,Double==1)
ecoli_working_double<-select(ecoli_working,double$new_allele)
single<-filter(ag_ecol,Singleton==1)
ecoli_working_single<-select(ecoli_working,single$new_allele)
tenner<-filter(ag_ecol,Tenner==1)
ecoli_working_tenner<-select(ecoli_working,tenner$new_allele)

tenner_stat<-filter(ag_ecol,Double==1) %>% select(new_allele) %>% distinct()


ecoli_rare_tenner<-rarefaction(ecoli_working_tenner,n.perm = 100)
ecoli_rare_double<-rarefaction(ecoli_working_double,n.perm = 100)
ecoli_rare_single<-rarefaction(ecoli_working_single,n.perm = 100)


for(i in 1:nrow(ecoli_rare_double)){
  print(i)
  ecoli_rare_double$b1[i]<-quantile(ecoli_rare_double[i,2:101],probs=c(0.025))
  ecoli_rare_double$b2[i]<-quantile(ecoli_rare_double[i,2:101],probs=c(0.975))
  
  ecoli_rare_single$b1[i]<-quantile(ecoli_rare_single[i,2:101],probs=c(0.025))
  ecoli_rare_single$b2[i]<-quantile(ecoli_rare_single[i,2:101],probs=c(0.975))
  
  ecoli_rare_tenner$b1[i]<-quantile(ecoli_rare_tenner[i,2:101],probs=c(0.025))
  ecoli_rare_tenner$b2[i]<-quantile(ecoli_rare_tenner[i,2:101],probs=c(0.975))
  
}

ecoli_rare_double$isolate<-NA
ecoli_rare_double$catalogue<-"Double_Rarefaction"
ecoli_rare_double$n<-NA
ecoli_rare_double<-select(ecoli_rare_double,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_double)<-c('isolate','g','catalogue','n','b1','b2')

ecoli_rare_single$isolate<-NA
ecoli_rare_single$catalogue<-"Single_Rarefaction"
ecoli_rare_single$n<-NA
ecoli_rare_single<-select(ecoli_rare_single,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_single)<-c('isolate','g','catalogue','n','b1','b2')

ecoli_rare_tenner$isolate<-NA
ecoli_rare_tenner$catalogue<-"Tenner_Rarefaction"
ecoli_rare_tenner$n<-NA
ecoli_rare_tenner<-select(ecoli_rare_tenner,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_tenner)<-c('isolate','g','catalogue','n','b1','b2')



out2<-out
out2$b1<-NA
out2$b2<-NA

out2<-rbind(out2,ecoli_rare_double,ecoli_rare_single,ecoli_rare_tenner)

out2$b1<-unlist(out2$b1)
out2$b2<-unlist(out2$b2)
out2<-filter(out2,catalogue != 'AMR Finder')

colors<-c( "#4BACCC","#239135","#DE1414")
out2$catalogue<-ifelse(out2$catalogue=="Allele observed at least once","Allele observed once",
                       ifelse(out2$catalogue=="Allele observed at least twice","Allele observed 2-9 times",
                              ifelse(out2$catalogue =="Allele observed at least ten times","Allele observed at least ten times",out2$catalogue)))
names(colors)<-c('Allele observed once','Allele observed 2-9 times','Allele observed at least ten times')

ecol_aminog_boot<-ggplot(out2) +
  aes(x=isolate,y=n,color=catalogue,group=catalogue) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 resistance associated allele")  + labs(color="Category") +
  scale_color_manual(values=colors) 

library(ggnewscale)
ecol_aminog_boot<-ecol_aminog_boot + new_scale_color()  


ecol_aminog_boot<-ecol_aminog_boot +
  geom_ribbon(data = filter(out2,catalogue=="Single_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2) + 
  geom_ribbon(data = filter(out2,catalogue=="Double_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2) +
  geom_ribbon(data = filter(out2,catalogue=="Tenner_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2)

ecol_aminog_boot_save<-ecol_aminog_boot

#######now we create the same thing but for drug classes, multiple drug classes on the same plot with one category


gene_classes<-read_tsv('./data/amino_gene_classes.tsv',col_names=F) 

###################amino##############
       
ag_out<-filter(out_save,gene %in% gene_classes$X1)



ag_ecol<-ag_out


count_table<-ag_ecol %>% distinct(guuid,new_allele,.keep_all = T) %>%  group_by(new_allele) %>% count() 
count_table$Singleton<-ifelse(count_table$n==1,1,0)
count_table$Double<-ifelse(count_table$n >=2 & count_table$n <10,1,0)
count_table$Tenner<-ifelse(count_table$n >=10,1,0)

count_table<-select(count_table,-n)

ag_ecol<-left_join(ag_ecol,count_table,by=c("new_allele"="new_allele"))


guuids<-unique(ag_ecol$guuid)
#dates<-read_tsv('~/gn/data/ecoli/Escherichia_guuids_match') %>% select(guuid,date)
#guuids<-data.frame(guuids)
#guuids<-left_join(guuids,dates,by=c("guuids"="guuid"))
#guuids<-arrange(guuids,date)
#guuids<-filter(guuids,!is.na(date))
guuids<-sample(guuids)
collection=NULL
out=NULL

for(g in guuids){
  print(g)
  p<-filter(ag_ecol,guuid==g)
  collection<-rbind(collection,data.frame(p))
  alleles<-length(unique(collection$gene))
  new_alleles<-length(unique(collection$new_allele))
  singletons<-filter(collection,Singleton==1)
  singletons<-length(unique(singletons$new_allele))
  confirmed<-filter(collection,Double==1)
  confirmed<-length(unique(confirmed$new_allele))
  tenner<-filter(collection,Tenner==1)
  tenner<-length(unique(tenner$new_allele))
  out<-rbind(out,data.frame(g,alleles,singletons,confirmed,tenner))
}

out$isolate<-1:nrow(out)
out<-select(out,isolate,g,everything())

out<-pivot_longer(out,names_to = "catalogue",values_to = "n",cols = 3:6)

out$catalogue<-case_when(out$catalogue == 'alleles' ~ 'AMR Finder',
                         out$catalogue == 'confirmed' ~ 'Allele observed at least twice',
                         out$catalogue == 'singletons' ~ 'Allele observed at least once',
                         out$catalogue == 'tenner' ~ 'Allele observed at least ten times')

colors<-c( "#4BACCC","#239135","#DE1414")
names(colors)<-c('Allele observed at least once','Allele observed at least twice','Allele observed at least ten times')
out<-filter(out,catalogue != 'AMR Finder')
ecol_aminog<-ggplot(out) +
  aes(x=isolate,y=n,color=catalogue,group=catalogue) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 resistance associated allele") + ggtitle("E. coli - Aminoglycosides")  + labs(color="Category") +
  scale_color_manual(values=colors) + geom_vline(xintercept = 300,linetype="dashed")

############bootstrap###############


ag_ecol<-filter(ag_ecol,guuid %in% guuids)

ecoli_working<-ag_ecol %>% select(guuid,new_allele)
names(ecoli_working)<-c('guuid','gene')

ecoli_working<-data.frame(table(ecoli_working)) %>% spread(gene,Freq)



ecoli_working<-select(ecoli_working,-guuid)

library(micropan)
double<-filter(ag_ecol,Double==1)
ecoli_working_double<-select(ecoli_working,double$new_allele)
single<-filter(ag_ecol,Singleton==1)
ecoli_working_single<-select(ecoli_working,single$new_allele)
tenner<-filter(ag_ecol,Tenner==1)
ecoli_working_tenner<-select(ecoli_working,tenner$new_allele)



ecoli_rare_tenner<-rarefaction(ecoli_working_tenner,n.perm = 100)
ecoli_rare_double<-rarefaction(ecoli_working_double,n.perm = 100)
ecoli_rare_single<-rarefaction(ecoli_working_single,n.perm = 100)



for(i in 1:nrow(ecoli_rare_double)){
  print(i)
  ecoli_rare_double$b1[i]<-quantile(ecoli_rare_double[i,2:101],probs=c(0.025))
  ecoli_rare_double$b2[i]<-quantile(ecoli_rare_double[i,2:101],probs=c(0.975))
  
  ecoli_rare_single$b1[i]<-quantile(ecoli_rare_single[i,2:101],probs=c(0.025))
  ecoli_rare_single$b2[i]<-quantile(ecoli_rare_single[i,2:101],probs=c(0.975))
  
  ecoli_rare_tenner$b1[i]<-quantile(ecoli_rare_tenner[i,2:101],probs=c(0.025))
  ecoli_rare_tenner$b2[i]<-quantile(ecoli_rare_tenner[i,2:101],probs=c(0.975))
  
}

ecoli_rare_double$isolate<-NA
ecoli_rare_double$catalogue<-"Double_Rarefaction"
ecoli_rare_double$n<-NA
ecoli_rare_double<-select(ecoli_rare_double,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_double)<-c('isolate','g','catalogue','n','b1','b2')

ecoli_rare_single$isolate<-NA
ecoli_rare_single$catalogue<-"Single_Rarefaction"
ecoli_rare_single$n<-NA
ecoli_rare_single<-select(ecoli_rare_single,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_single)<-c('isolate','g','catalogue','n','b1','b2')

ecoli_rare_tenner$isolate<-NA
ecoli_rare_tenner$catalogue<-"Tenner_Rarefaction"
ecoli_rare_tenner$n<-NA
ecoli_rare_tenner<-select(ecoli_rare_tenner,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_tenner)<-c('isolate','g','catalogue','n','b1','b2')



out2<-out
out2$b1<-NA
out2$b2<-NA

out2<-rbind(out2,ecoli_rare_double,ecoli_rare_single,ecoli_rare_tenner)

out2$b1<-unlist(out2$b1)
out2$b2<-unlist(out2$b2)
ecol_aminog_boot<-ggplot(out2) +
  aes(x=isolate,y=n,color=catalogue,group=catalogue) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 resistance associated allele") + ggtitle("E. coli - Aminoglycosides - 74% (61-86)")  + labs(color="Category") +
  scale_color_manual(values=colors) + geom_vline(xintercept = 300,linetype="dashed") 

library(ggnewscale)
ecol_aminog_boot<-ecol_aminog_boot + new_scale_color()  


ecol_aminog_boot<-ecol_aminog_boot +
  geom_ribbon(data = filter(out2,catalogue=="Single_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2) + 
  geom_ribbon(data = filter(out2,catalogue=="Double_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2) +
  geom_ribbon(data = filter(out2,catalogue=="Tenner_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2)

aminoglycosides<-out2


########beta-lactam############

gene_classes<-read_tsv('./data/betalactam_gene_classes.tsv',col_names=F) 



ag_out<-filter(out_save,gene %in% gene_classes$X1)


ag_ecol<-ag_out


count_table<-ag_ecol %>% distinct(guuid,new_allele,.keep_all = T) %>%  group_by(new_allele) %>% count() 
count_table$Singleton<-ifelse(count_table$n==1,1,0)
count_table$Double<-ifelse(count_table$n >=2 & count_table$n <10,1,0)
count_table$Tenner<-ifelse(count_table$n >=10,1,0)


count_table<-select(count_table,-n)

ag_ecol<-left_join(ag_ecol,count_table,by=c("new_allele"="new_allele"))


guuids<-unique(ag_ecol$guuid)
#dates<-read_tsv('~/gn/data/ecoli/Escherichia_guuids_match') %>% select(guuid,date)
#guuids<-data.frame(guuids)
#guuids<-left_join(guuids,dates,by=c("guuids"="guuid"))
#guuids<-arrange(guuids,date)
#guuids<-filter(guuids,!is.na(date))
guuids<-sample(guuids)
collection=NULL
out=NULL

for(g in guuids){
  print(g)
  p<-filter(ag_ecol,guuid==g)
  collection<-rbind(collection,data.frame(p))
  alleles<-length(unique(collection$gene))
  new_alleles<-length(unique(collection$new_allele))
  singletons<-filter(collection,Singleton==1)
  singletons<-length(unique(singletons$new_allele))
  confirmed<-filter(collection,Double==1)
  confirmed<-length(unique(confirmed$new_allele))
  tenner<-filter(collection,Tenner==1)
  tenner<-length(unique(tenner$new_allele))
  out<-rbind(out,data.frame(g,alleles,singletons,confirmed,tenner))
}

out$isolate<-1:nrow(out)
out<-select(out,isolate,g,everything())

out<-pivot_longer(out,names_to = "catalogue",values_to = "n",cols = 3:6)

out$catalogue<-case_when(out$catalogue == 'alleles' ~ 'AMR Finder',
                         out$catalogue == 'confirmed' ~ 'Allele observed at least twice',
                         out$catalogue == 'singletons' ~ 'Allele observed at least once',
                         out$catalogue == 'tenner' ~ 'Allele observed at least ten times')

colors<-c( "#4BACCC","#239135","#DE1414")
names(colors)<-c('Allele observed at least once','Allele observed at least twice','Allele observed at least ten times')
out<-filter(out,catalogue != 'AMR Finder')
ecol_aminog<-ggplot(out) +
  aes(x=isolate,y=n,color=catalogue,group=catalogue) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 resistance associated allele") + ggtitle("E. coli - Aminoglycosides")  + labs(color="Category") +
  scale_color_manual(values=colors) + geom_vline(xintercept = 300,linetype="dashed")

############bootstrap###############


ag_ecol<-filter(ag_ecol,guuid %in% guuids)

ecoli_working<-ag_ecol %>% select(guuid,new_allele)
names(ecoli_working)<-c('guuid','gene')

ecoli_working<-data.frame(table(ecoli_working)) %>% spread(gene,Freq)



ecoli_working<-select(ecoli_working,-guuid)

library(micropan)
double<-filter(ag_ecol,Double==1)
ecoli_working_double<-select(ecoli_working,double$new_allele)
single<-filter(ag_ecol,Singleton==1)
ecoli_working_single<-select(ecoli_working,single$new_allele)
tenner<-filter(ag_ecol,Tenner==1)
ecoli_working_tenner<-select(ecoli_working,tenner$new_allele)



ecoli_rare_tenner<-rarefaction(ecoli_working_tenner,n.perm = 100)
ecoli_rare_double<-rarefaction(ecoli_working_double,n.perm = 100)
ecoli_rare_single<-rarefaction(ecoli_working_single,n.perm = 100)



for(i in 1:nrow(ecoli_rare_double)){
  print(i)
  ecoli_rare_double$b1[i]<-quantile(ecoli_rare_double[i,2:101],probs=c(0.025))
  ecoli_rare_double$b2[i]<-quantile(ecoli_rare_double[i,2:101],probs=c(0.975))
  
  ecoli_rare_single$b1[i]<-quantile(ecoli_rare_single[i,2:101],probs=c(0.025))
  ecoli_rare_single$b2[i]<-quantile(ecoli_rare_single[i,2:101],probs=c(0.975))
  
  ecoli_rare_tenner$b1[i]<-quantile(ecoli_rare_tenner[i,2:101],probs=c(0.025))
  ecoli_rare_tenner$b2[i]<-quantile(ecoli_rare_tenner[i,2:101],probs=c(0.975))
  
}

ecoli_rare_double$isolate<-NA
ecoli_rare_double$catalogue<-"Double_Rarefaction"
ecoli_rare_double$n<-NA
ecoli_rare_double<-select(ecoli_rare_double,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_double)<-c('isolate','g','catalogue','n','b1','b2')

ecoli_rare_single$isolate<-NA
ecoli_rare_single$catalogue<-"Single_Rarefaction"
ecoli_rare_single$n<-NA
ecoli_rare_single<-select(ecoli_rare_single,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_single)<-c('isolate','g','catalogue','n','b1','b2')

ecoli_rare_tenner$isolate<-NA
ecoli_rare_tenner$catalogue<-"Tenner_Rarefaction"
ecoli_rare_tenner$n<-NA
ecoli_rare_tenner<-select(ecoli_rare_tenner,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_tenner)<-c('isolate','g','catalogue','n','b1','b2')



out2<-out
out2$b1<-NA
out2$b2<-NA

out2<-rbind(out2,ecoli_rare_double,ecoli_rare_single,ecoli_rare_tenner)

out2$b1<-unlist(out2$b1)
out2$b2<-unlist(out2$b2)
ecol_aminog_boot<-ggplot(out2) +
  aes(x=isolate,y=n,color=catalogue,group=catalogue) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 resistance associated allele") + ggtitle("E. coli - Aminoglycosides - 74% (61-86)")  + labs(color="Category") +
  scale_color_manual(values=colors) + geom_vline(xintercept = 300,linetype="dashed") 

library(ggnewscale)
ecol_aminog_boot<-ecol_aminog_boot + new_scale_color()  


ecol_aminog_boot<-ecol_aminog_boot +
  geom_ribbon(data = filter(out2,catalogue=="Single_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2) + 
  geom_ribbon(data = filter(out2,catalogue=="Double_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2) +
  geom_ribbon(data = filter(out2,catalogue=="Tenner_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2)

beta_lactam<-out2


############quinolone##########
gene_classes<-read_tsv('./data/quinolone_gene_classes.tsv',col_names=F) 



  ag_out<-filter(out_save,gene %in% gene_classes$X1)

ag_ecol<-ag_out


count_table<-ag_ecol %>% distinct(guuid,new_allele,.keep_all = T) %>%  group_by(new_allele) %>% count() 
count_table$Singleton<-ifelse(count_table$n==1,1,0)
count_table$Double<-ifelse(count_table$n >=2 & count_table$n <10,1,0)
count_table$Tenner<-ifelse(count_table$n >=10,1,0)


count_table<-select(count_table,-n)

ag_ecol<-left_join(ag_ecol,count_table,by=c("new_allele"="new_allele"))


guuids<-unique(ag_ecol$guuid)
#dates<-read_tsv('~/gn/data/ecoli/Escherichia_guuids_match') %>% select(guuid,date)
#guuids<-data.frame(guuids)
#guuids<-left_join(guuids,dates,by=c("guuids"="guuid"))
#guuids<-arrange(guuids,date)
#guuids<-filter(guuids,!is.na(date))
guuids<-sample(guuids)
collection=NULL
out=NULL

for(g in guuids){
  print(g)
  p<-filter(ag_ecol,guuid==g)
  collection<-rbind(collection,data.frame(p))
  alleles<-length(unique(collection$gene))
  new_alleles<-length(unique(collection$new_allele))
  singletons<-filter(collection,Singleton==1)
  singletons<-length(unique(singletons$new_allele))
  confirmed<-filter(collection,Double==1)
  confirmed<-length(unique(confirmed$new_allele))
  tenner<-filter(collection,Tenner==1)
  tenner<-length(unique(tenner$new_allele))
  out<-rbind(out,data.frame(g,alleles,singletons,confirmed,tenner))
}

out$isolate<-1:nrow(out)
out<-select(out,isolate,g,everything())

out<-pivot_longer(out,names_to = "catalogue",values_to = "n",cols = 3:6)

out$catalogue<-case_when(out$catalogue == 'alleles' ~ 'AMR Finder',
                         out$catalogue == 'confirmed' ~ 'Allele observed at least twice',
                         out$catalogue == 'singletons' ~ 'Allele observed at least once',
                         out$catalogue == 'tenner' ~ 'Allele observed at least ten times')

colors<-c( "#4BACCC","#239135","#DE1414")
names(colors)<-c('Allele observed at least once','Allele observed at least twice','Allele observed at least ten times')
out<-filter(out,catalogue != 'AMR Finder')
ecol_aminog<-ggplot(out) +
  aes(x=isolate,y=n,color=catalogue,group=catalogue) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 resistance associated allele") + ggtitle("E. coli - Aminoglycosides")  + labs(color="Category") +
  scale_color_manual(values=colors) + geom_vline(xintercept = 300,linetype="dashed")

ecol_aminog
############bootstrap###############


ag_ecol<-filter(ag_ecol,guuid %in% guuids)

ecoli_working<-ag_ecol %>% select(guuid,new_allele)
names(ecoli_working)<-c('guuid','gene')

ecoli_working<-data.frame(table(ecoli_working)) %>% spread(gene,Freq)



ecoli_working<-select(ecoli_working,-guuid)

library(micropan)
double<-filter(ag_ecol,Double==1)
ecoli_working_double<-select(ecoli_working,double$new_allele)
single<-filter(ag_ecol,Singleton==1)
ecoli_working_single<-select(ecoli_working,single$new_allele)
tenner<-filter(ag_ecol,Tenner==1)
ecoli_working_tenner<-select(ecoli_working,tenner$new_allele)



ecoli_rare_tenner<-rarefaction(ecoli_working_tenner,n.perm = 100)
ecoli_rare_double<-rarefaction(ecoli_working_double,n.perm = 100)
ecoli_rare_single<-rarefaction(ecoli_working_single,n.perm = 100)



for(i in 1:nrow(ecoli_rare_double)){
  print(i)
  ecoli_rare_double$b1[i]<-quantile(ecoli_rare_double[i,2:101],probs=c(0.025))
  ecoli_rare_double$b2[i]<-quantile(ecoli_rare_double[i,2:101],probs=c(0.975))
  
  ecoli_rare_single$b1[i]<-quantile(ecoli_rare_single[i,2:101],probs=c(0.025))
  ecoli_rare_single$b2[i]<-quantile(ecoli_rare_single[i,2:101],probs=c(0.975))
  
  ecoli_rare_tenner$b1[i]<-quantile(ecoli_rare_tenner[i,2:101],probs=c(0.025))
  ecoli_rare_tenner$b2[i]<-quantile(ecoli_rare_tenner[i,2:101],probs=c(0.975))
  
}

ecoli_rare_double$isolate<-NA
ecoli_rare_double$catalogue<-"Double_Rarefaction"
ecoli_rare_double$n<-NA
ecoli_rare_double<-select(ecoli_rare_double,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_double)<-c('isolate','g','catalogue','n','b1','b2')

ecoli_rare_single$isolate<-NA
ecoli_rare_single$catalogue<-"Single_Rarefaction"
ecoli_rare_single$n<-NA
ecoli_rare_single<-select(ecoli_rare_single,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_single)<-c('isolate','g','catalogue','n','b1','b2')

ecoli_rare_tenner$isolate<-NA
ecoli_rare_tenner$catalogue<-"Tenner_Rarefaction"
ecoli_rare_tenner$n<-NA
ecoli_rare_tenner<-select(ecoli_rare_tenner,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_tenner)<-c('isolate','g','catalogue','n','b1','b2')



out2<-out
out2$b1<-NA
out2$b2<-NA

out2<-rbind(out2,ecoli_rare_double,ecoli_rare_single,ecoli_rare_tenner)

out2$b1<-unlist(out2$b1)
out2$b2<-unlist(out2$b2)
ecol_aminog_boot<-ggplot(out2) +
  aes(x=isolate,y=n,color=catalogue,group=catalogue) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 resistance associated allele") + ggtitle("E. coli - Aminoglycosides - 74% (61-86)")  + labs(color="Category") +
  scale_color_manual(values=colors) + geom_vline(xintercept = 300,linetype="dashed") 

library(ggnewscale)
ecol_aminog_boot<-ecol_aminog_boot + new_scale_color()  


ecol_aminog_boot<-ecol_aminog_boot +
  geom_ribbon(data = filter(out2,catalogue=="Single_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2) + 
  geom_ribbon(data = filter(out2,catalogue=="Double_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2) +
  geom_ribbon(data = filter(out2,catalogue=="Tenner_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2)

quinolone<-out2


###########Cephalosporin###########

gene_classes<-read_tsv('./data/cephalosporin_gene_classes.tsv',col_names=F) 



ag_out<-filter(out_save,gene %in% gene_classes$X1)


ag_ecol<-ag_out


count_table<-ag_ecol %>% distinct(guuid,new_allele,.keep_all = T) %>%  group_by(new_allele) %>% count() 
count_table$Singleton<-ifelse(count_table$n==1,1,0)
count_table$Double<-ifelse(count_table$n >=2 & count_table$n <10,1,0)
count_table$Tenner<-ifelse(count_table$n >=10,1,0)

count_table<-select(count_table,-n)

ag_ecol<-left_join(ag_ecol,count_table,by=c("new_allele"="new_allele"))


guuids<-unique(ag_ecol$guuid)
#dates<-read_tsv('~/gn/data/ecoli/Escherichia_guuids_match') %>% select(guuid,date)
#guuids<-data.frame(guuids)
#guuids<-left_join(guuids,dates,by=c("guuids"="guuid"))
#guuids<-arrange(guuids,date)
#guuids<-filter(guuids,!is.na(date))
guuids<-sample(guuids)
collection=NULL
out=NULL

for(g in guuids){
  print(g)
  p<-filter(ag_ecol,guuid==g)
  collection<-rbind(collection,data.frame(p))
  alleles<-length(unique(collection$gene))
  new_alleles<-length(unique(collection$new_allele))
  singletons<-filter(collection,Singleton==1)
  singletons<-length(unique(singletons$new_allele))
  confirmed<-filter(collection,Double==1)
  confirmed<-length(unique(confirmed$new_allele))
  tenner<-filter(collection,Tenner==1)
  tenner<-length(unique(tenner$new_allele))
  out<-rbind(out,data.frame(g,alleles,singletons,confirmed,tenner))
}

out$isolate<-1:nrow(out)
out<-select(out,isolate,g,everything())

out<-pivot_longer(out,names_to = "catalogue",values_to = "n",cols = 3:6)

out$catalogue<-case_when(out$catalogue == 'alleles' ~ 'AMR Finder',
                         out$catalogue == 'confirmed' ~ 'Allele observed at least twice',
                         out$catalogue == 'singletons' ~ 'Allele observed at least once',
                         out$catalogue == 'tenner' ~ 'Allele observed at least ten times')

colors<-c( "#4BACCC","#239135","#DE1414")
names(colors)<-c('Allele observed at least once','Allele observed at least twice','Allele observed at least ten times')
out<-filter(out,catalogue != 'AMR Finder')
ecol_aminog<-ggplot(out) +
  aes(x=isolate,y=n,color=catalogue,group=catalogue) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 resistance associated allele") + ggtitle("E. coli - Aminoglycosides")  + labs(color="Category") +
  scale_color_manual(values=colors) + geom_vline(xintercept = 300,linetype="dashed")

ecol_aminog
############bootstrap###############


ag_ecol<-filter(ag_ecol,guuid %in% guuids)

ecoli_working<-ag_ecol %>% select(guuid,new_allele)
names(ecoli_working)<-c('guuid','gene')

ecoli_working<-data.frame(table(ecoli_working)) %>% spread(gene,Freq)



ecoli_working<-select(ecoli_working,-guuid)

library(micropan)
double<-filter(ag_ecol,Double==1)
ecoli_working_double<-select(ecoli_working,double$new_allele)
single<-filter(ag_ecol,Singleton==1)
ecoli_working_single<-select(ecoli_working,single$new_allele)
tenner<-filter(ag_ecol,Tenner==1)
ecoli_working_tenner<-select(ecoli_working,tenner$new_allele)



ecoli_rare_tenner<-rarefaction(ecoli_working_tenner,n.perm = 100)
ecoli_rare_double<-rarefaction(ecoli_working_double,n.perm = 100)
ecoli_rare_single<-rarefaction(ecoli_working_single,n.perm = 100)



for(i in 1:nrow(ecoli_rare_double)){
  print(i)
  ecoli_rare_double$b1[i]<-quantile(ecoli_rare_double[i,2:101],probs=c(0.025))
  ecoli_rare_double$b2[i]<-quantile(ecoli_rare_double[i,2:101],probs=c(0.975))
  
  ecoli_rare_single$b1[i]<-quantile(ecoli_rare_single[i,2:101],probs=c(0.025))
  ecoli_rare_single$b2[i]<-quantile(ecoli_rare_single[i,2:101],probs=c(0.975))
  
  ecoli_rare_tenner$b1[i]<-quantile(ecoli_rare_tenner[i,2:101],probs=c(0.025))
  ecoli_rare_tenner$b2[i]<-quantile(ecoli_rare_tenner[i,2:101],probs=c(0.975))
  
}

ecoli_rare_double$isolate<-NA
ecoli_rare_double$catalogue<-"Double_Rarefaction"
ecoli_rare_double$n<-NA
ecoli_rare_double<-select(ecoli_rare_double,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_double)<-c('isolate','g','catalogue','n','b1','b2')

ecoli_rare_single$isolate<-NA
ecoli_rare_single$catalogue<-"Single_Rarefaction"
ecoli_rare_single$n<-NA
ecoli_rare_single<-select(ecoli_rare_single,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_single)<-c('isolate','g','catalogue','n','b1','b2')

ecoli_rare_tenner$isolate<-NA
ecoli_rare_tenner$catalogue<-"Tenner_Rarefaction"
ecoli_rare_tenner$n<-NA
ecoli_rare_tenner<-select(ecoli_rare_tenner,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_tenner)<-c('isolate','g','catalogue','n','b1','b2')



out2<-out
out2$b1<-NA
out2$b2<-NA

out2<-rbind(out2,ecoli_rare_double,ecoli_rare_single,ecoli_rare_tenner)

out2$b1<-unlist(out2$b1)
out2$b2<-unlist(out2$b2)
ecol_aminog_boot<-ggplot(out2) +
  aes(x=isolate,y=n,color=catalogue,group=catalogue) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 resistance associated allele") + ggtitle("E. coli - Aminoglycosides - 74% (61-86)")  + labs(color="Category") +
  scale_color_manual(values=colors) + geom_vline(xintercept = 300,linetype="dashed") 

library(ggnewscale)
ecol_aminog_boot<-ecol_aminog_boot + new_scale_color()  


ecol_aminog_boot<-ecol_aminog_boot +
  geom_ribbon(data = filter(out2,catalogue=="Single_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2) + 
  geom_ribbon(data = filter(out2,catalogue=="Double_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2) +
  geom_ribbon(data = filter(out2,catalogue=="Tenner_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2)

cephalosporin<-out2

#########fosfomycin###########

gene_classes<-read_tsv('./data/fosfomycin_gene_classes.tsv',col_names=F) 



ag_out<-filter(out_save,gene %in% gene_classes$X1)


ag_ecol<-ag_out


count_table<-ag_ecol %>% distinct(guuid,new_allele,.keep_all = T) %>%  group_by(new_allele) %>% count() 
count_table$Singleton<-ifelse(count_table$n==1,1,0)
count_table$Double<-ifelse(count_table$n >=2 & count_table$n <10,1,0)
count_table$Tenner<-ifelse(count_table$n >=10,1,0)

count_table<-select(count_table,-n)

ag_ecol<-left_join(ag_ecol,count_table,by=c("new_allele"="new_allele"))


guuids<-unique(ag_ecol$guuid)
#dates<-read_tsv('~/gn/data/ecoli/Escherichia_guuids_match') %>% select(guuid,date)
#guuids<-data.frame(guuids)
#guuids<-left_join(guuids,dates,by=c("guuids"="guuid"))
#guuids<-arrange(guuids,date)
#guuids<-filter(guuids,!is.na(date))
guuids<-sample(guuids)
collection=NULL
out=NULL

for(g in guuids){
  print(g)
  p<-filter(ag_ecol,guuid==g)
  collection<-rbind(collection,data.frame(p))
  alleles<-length(unique(collection$gene))
  new_alleles<-length(unique(collection$new_allele))
  singletons<-filter(collection,Singleton==1)
  singletons<-length(unique(singletons$new_allele))
  confirmed<-filter(collection,Double==1)
  confirmed<-length(unique(confirmed$new_allele))
  tenner<-filter(collection,Tenner==1)
  tenner<-length(unique(tenner$new_allele))
  out<-rbind(out,data.frame(g,alleles,singletons,confirmed,tenner))
}

out$isolate<-1:nrow(out)
out<-select(out,isolate,g,everything())

out<-pivot_longer(out,names_to = "catalogue",values_to = "n",cols = 3:6)

out$catalogue<-case_when(out$catalogue == 'alleles' ~ 'AMR Finder',
                         out$catalogue == 'confirmed' ~ 'Allele observed at least twice',
                         out$catalogue == 'singletons' ~ 'Allele observed at least once',
                         out$catalogue == 'tenner' ~ 'Allele observed at least ten times')

colors<-c( "#4BACCC","#239135","#DE1414")
names(colors)<-c('Allele observed at least once','Allele observed at least twice','Allele observed at least ten times')
out<-filter(out,catalogue != 'AMR Finder')
ecol_aminog<-ggplot(out) +
  aes(x=isolate,y=n,color=catalogue,group=catalogue) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 resistance associated allele") + ggtitle("E. coli - Aminoglycosides")  + labs(color="Category") +
  scale_color_manual(values=colors) + geom_vline(xintercept = 300,linetype="dashed")

ecol_aminog
############bootstrap###############


ag_ecol<-filter(ag_ecol,guuid %in% guuids)

ecoli_working<-ag_ecol %>% select(guuid,new_allele)
names(ecoli_working)<-c('guuid','gene')

ecoli_working<-data.frame(table(ecoli_working)) %>% spread(gene,Freq)



ecoli_working<-select(ecoli_working,-guuid)

library(micropan)
double<-filter(ag_ecol,Double==1)
ecoli_working_double<-select(ecoli_working,double$new_allele)
single<-filter(ag_ecol,Singleton==1)
ecoli_working_single<-select(ecoli_working,single$new_allele)
tenner<-filter(ag_ecol,Tenner==1)
ecoli_working_tenner<-select(ecoli_working,tenner$new_allele)



ecoli_rare_tenner<-rarefaction(ecoli_working_tenner,n.perm = 100)
ecoli_rare_double<-rarefaction(ecoli_working_double,n.perm = 100)
ecoli_rare_single<-rarefaction(ecoli_working_single,n.perm = 100)



for(i in 1:nrow(ecoli_rare_double)){
  print(i)
  ecoli_rare_double$b1[i]<-quantile(ecoli_rare_double[i,2:101],probs=c(0.025))
  ecoli_rare_double$b2[i]<-quantile(ecoli_rare_double[i,2:101],probs=c(0.975))
  
  ecoli_rare_single$b1[i]<-quantile(ecoli_rare_single[i,2:101],probs=c(0.025))
  ecoli_rare_single$b2[i]<-quantile(ecoli_rare_single[i,2:101],probs=c(0.975))
  
  ecoli_rare_tenner$b1[i]<-quantile(ecoli_rare_tenner[i,2:101],probs=c(0.025))
  ecoli_rare_tenner$b2[i]<-quantile(ecoli_rare_tenner[i,2:101],probs=c(0.975))
  
}

ecoli_rare_double$isolate<-NA
ecoli_rare_double$catalogue<-"Double_Rarefaction"
ecoli_rare_double$n<-NA
ecoli_rare_double<-select(ecoli_rare_double,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_double)<-c('isolate','g','catalogue','n','b1','b2')

ecoli_rare_single$isolate<-NA
ecoli_rare_single$catalogue<-"Single_Rarefaction"
ecoli_rare_single$n<-NA
ecoli_rare_single<-select(ecoli_rare_single,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_single)<-c('isolate','g','catalogue','n','b1','b2')

ecoli_rare_tenner$isolate<-NA
ecoli_rare_tenner$catalogue<-"Tenner_Rarefaction"
ecoli_rare_tenner$n<-NA
ecoli_rare_tenner<-select(ecoli_rare_tenner,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_tenner)<-c('isolate','g','catalogue','n','b1','b2')



out2<-out
out2$b1<-NA
out2$b2<-NA

out2<-rbind(out2,ecoli_rare_double,ecoli_rare_single,ecoli_rare_tenner)

out2$b1<-unlist(out2$b1)
out2$b2<-unlist(out2$b2)
ecol_aminog_boot<-ggplot(out2) +
  aes(x=isolate,y=n,color=catalogue,group=catalogue) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 resistance associated allele") + ggtitle("E. coli - Aminoglycosides - 74% (61-86)")  + labs(color="Category") +
  scale_color_manual(values=colors) + geom_vline(xintercept = 300,linetype="dashed") 

library(ggnewscale)
ecol_aminog_boot<-ecol_aminog_boot + new_scale_color()  


ecol_aminog_boot<-ecol_aminog_boot +
  geom_ribbon(data = filter(out2,catalogue=="Single_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2) + 
  geom_ribbon(data = filter(out2,catalogue=="Double_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2) +
  geom_ribbon(data = filter(out2,catalogue=="Tenner_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2)

fosfomycin<-out2

##############trimethoprim#########

gene_classes<-read_tsv('./data/trimethoprim_gene_classes.tsv',col_names=F) 



  ag_out<-filter(out_save,gene %in% gene_classes$X1)


ag_ecol<-ag_out


count_table<-ag_ecol %>% distinct(guuid,new_allele,.keep_all = T) %>%  group_by(new_allele) %>% count() 
count_table$Singleton<-ifelse(count_table$n==1,1,0)
count_table$Double<-ifelse(count_table$n >=2 & count_table$n <10,1,0)
count_table$Tenner<-ifelse(count_table$n >=10,1,0)

count_table<-select(count_table,-n)

ag_ecol<-left_join(ag_ecol,count_table,by=c("new_allele"="new_allele"))


guuids<-unique(ag_ecol$guuid)
#dates<-read_tsv('~/gn/data/ecoli/Escherichia_guuids_match') %>% select(guuid,date)
#guuids<-data.frame(guuids)
#guuids<-left_join(guuids,dates,by=c("guuids"="guuid"))
#guuids<-arrange(guuids,date)
#guuids<-filter(guuids,!is.na(date))
guuids<-sample(guuids)
collection=NULL
out=NULL

for(g in guuids){
  print(g)
  p<-filter(ag_ecol,guuid==g)
  collection<-rbind(collection,data.frame(p))
  alleles<-length(unique(collection$gene))
  new_alleles<-length(unique(collection$new_allele))
  singletons<-filter(collection,Singleton==1)
  singletons<-length(unique(singletons$new_allele))
  confirmed<-filter(collection,Double==1)
  confirmed<-length(unique(confirmed$new_allele))
  tenner<-filter(collection,Tenner==1)
  tenner<-length(unique(tenner$new_allele))
  out<-rbind(out,data.frame(g,alleles,singletons,confirmed,tenner))
}

out$isolate<-1:nrow(out)
out<-select(out,isolate,g,everything())

out<-pivot_longer(out,names_to = "catalogue",values_to = "n",cols = 3:6)

out$catalogue<-case_when(out$catalogue == 'alleles' ~ 'AMR Finder',
                         out$catalogue == 'confirmed' ~ 'Allele observed at least twice',
                         out$catalogue == 'singletons' ~ 'Allele observed at least once',
                         out$catalogue == 'tenner' ~ 'Allele observed at least ten times')

colors<-c( "#4BACCC","#239135","#DE1414")
names(colors)<-c('Allele observed at least once','Allele observed at least twice','Allele observed at least ten times')
out<-filter(out,catalogue != 'AMR Finder')
ecol_aminog<-ggplot(out) +
  aes(x=isolate,y=n,color=catalogue,group=catalogue) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 resistance associated allele") + ggtitle("E. coli - Aminoglycosides")  + labs(color="Category") +
  scale_color_manual(values=colors) + geom_vline(xintercept = 300,linetype="dashed")

ecol_aminog
############bootstrap###############


ag_ecol<-filter(ag_ecol,guuid %in% guuids)

ecoli_working<-ag_ecol %>% select(guuid,new_allele)
names(ecoli_working)<-c('guuid','gene')

ecoli_working<-data.frame(table(ecoli_working)) %>% spread(gene,Freq)



ecoli_working<-select(ecoli_working,-guuid)

library(micropan)
double<-filter(ag_ecol,Double==1)
ecoli_working_double<-select(ecoli_working,double$new_allele)
single<-filter(ag_ecol,Singleton==1)
ecoli_working_single<-select(ecoli_working,single$new_allele)
tenner<-filter(ag_ecol,Tenner==1)
ecoli_working_tenner<-select(ecoli_working,tenner$new_allele)



ecoli_rare_tenner<-rarefaction(ecoli_working_tenner,n.perm = 100)
ecoli_rare_double<-rarefaction(ecoli_working_double,n.perm = 100)
ecoli_rare_single<-rarefaction(ecoli_working_single,n.perm = 100)



for(i in 1:nrow(ecoli_rare_double)){
  print(i)
  ecoli_rare_double$b1[i]<-quantile(ecoli_rare_double[i,2:101],probs=c(0.025))
  ecoli_rare_double$b2[i]<-quantile(ecoli_rare_double[i,2:101],probs=c(0.975))
  
  ecoli_rare_single$b1[i]<-quantile(ecoli_rare_single[i,2:101],probs=c(0.025))
  ecoli_rare_single$b2[i]<-quantile(ecoli_rare_single[i,2:101],probs=c(0.975))
  
  ecoli_rare_tenner$b1[i]<-quantile(ecoli_rare_tenner[i,2:101],probs=c(0.025))
  ecoli_rare_tenner$b2[i]<-quantile(ecoli_rare_tenner[i,2:101],probs=c(0.975))
  
}

ecoli_rare_double$isolate<-NA
ecoli_rare_double$catalogue<-"Double_Rarefaction"
ecoli_rare_double$n<-NA
ecoli_rare_double<-select(ecoli_rare_double,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_double)<-c('isolate','g','catalogue','n','b1','b2')

ecoli_rare_single$isolate<-NA
ecoli_rare_single$catalogue<-"Single_Rarefaction"
ecoli_rare_single$n<-NA
ecoli_rare_single<-select(ecoli_rare_single,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_single)<-c('isolate','g','catalogue','n','b1','b2')

ecoli_rare_tenner$isolate<-NA
ecoli_rare_tenner$catalogue<-"Tenner_Rarefaction"
ecoli_rare_tenner$n<-NA
ecoli_rare_tenner<-select(ecoli_rare_tenner,Genome,isolate,catalogue,n,b1,b2)
names(ecoli_rare_tenner)<-c('isolate','g','catalogue','n','b1','b2')



out2<-out
out2$b1<-NA
out2$b2<-NA

out2<-rbind(out2,ecoli_rare_double,ecoli_rare_single,ecoli_rare_tenner)

out2$b1<-unlist(out2$b1)
out2$b2<-unlist(out2$b2)
ecol_aminog_boot<-ggplot(out2) +
  aes(x=isolate,y=n,color=catalogue,group=catalogue) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 resistance associated allele") + ggtitle("E. coli - Aminoglycosides - 74% (61-86)")  + labs(color="Category") +
  scale_color_manual(values=colors) + geom_vline(xintercept = 300,linetype="dashed") 

library(ggnewscale)
ecol_aminog_boot<-ecol_aminog_boot + new_scale_color()  


ecol_aminog_boot<-ecol_aminog_boot +
  geom_ribbon(data = filter(out2,catalogue=="Single_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2) + 
  geom_ribbon(data = filter(out2,catalogue=="Double_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2) +
  geom_ribbon(data = filter(out2,catalogue=="Tenner_Rarefaction"), aes(x=isolate,ymin=b1,ymax=b2),alpha=0.2)

trimethoprim<-out2

###########put it together#########

trimethoprim$drug<-"Trimethoprim"
fosfomycin$drug<-"Fosfomycin"
cephalosporin$drug<-"Cephalosporin"
beta_lactam$drug<-"Beta-lactam"
aminoglycosides$drug<-"Aminoglycoside"
quinolone$drug<-"Quinolone"

all<-rbind(beta_lactam,aminoglycosides,quinolone,cephalosporin,fosfomycin,trimethoprim)
all$catalogue<-ifelse(all$catalogue=="Allele observed at least once","Allele observed once",
                      ifelse(all$catalogue == "Allele observed at least twice","Allele observed 2-9 times",
                             ifelse(all$catalogue == "Allele observed at least ten times","Allele observed at least ten times",all$catalogue)))
double<-filter(all,catalogue =='Allele observed 2-9 times' | catalogue == 'Double_Rarefaction')

colors<-c("#AD3B3B" ,"#205719","#2D48E0", "#32C2B1" ,"#F79C34", "#6D5EAB")
names(colors)<-c('Aminoglycoside','Beta-lactam','Quinolone','Cephalosporin','Fosfomycin','Trimethoprim')
ecol_double<-ggplot(filter(double,catalogue=="Allele observed 2-9 times")) +
  aes(x=isolate,y=n,color=drug,group=drug) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 ARG allele") + ggtitle("ARG allele observed 2-9 times")  + labs(color="Category") +
  scale_color_manual(values=colors) 

ten<-filter(all,catalogue =='Allele observed at least ten times' | catalogue == 'Tenner_Rarefaction')
ecol_ten<-ggplot(filter(ten,catalogue=="Allele observed at least ten times")) +
  aes(x=isolate,y=n,color=drug,group=drug) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 ARG allele") + ggtitle("ARG allele observed at least ten times")  + labs(color="Category") +
  scale_color_manual(values=colors) 

single<-filter(all,catalogue == 'Allele observed once' | catalogue == 'Single_Rarefacetion')
ecol_single<-ggplot(filter(single,catalogue=="Allele observed once")) +
  aes(x=isolate,y=n,color=drug,group=drug) +
  geom_line(size=2) + theme_minimal() + ylab("Number of unique alleles") + xlab("Number of isolates with >=1 ARG allele") + ggtitle("ARG allele observed once")  + labs(color="Category") +
  scale_color_manual(values=colors) 


ecol_aminog_boot_save /(ecol_single + ecol_double + ecol_ten) + plot_layout(guides = "collect")
############explainable res###########

############ ER gentamicin ############


all_guuids<-read_tsv('./data/guuids',col_names = c("guuid"))

res<-read_tsv('./data/global_pheno.tsv')

res<-filter(res,guuid %in% all_guuids$guuid)
all_guuids<-filter(all_guuids,guuid %in% res$guuid)

guuids<-left_join(all_guuids,res,by=c("guuid"))

#comment out as desired to get exact match or default settings

#amrfinder<-read_tsv("./data/prediction/amrfinder_gladstone_merger.tsv") %>% filter(`% Identity to reference sequence` ==100 & `% Coverage of reference sequence` ==100)
amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv")

ag_res<-filter(guuids,Gentamicin==1)
ag_counts<-filter(amrfinder,grepl('GENTAMICIN',Subclass)) %>% filter(guuid %in% ag_res$guuid) %>% group_by(`Gene symbol`) %>% count() %>% arrange(desc(n))
ag_genes<-unique(ag_counts$`Gene symbol`)
ag_genes<-ag_genes[!grepl('blaEC',ag_genes)]

ag_amrfinder<-filter(amrfinder,guuid %in% guuids$guuid)


###### here we calculate stats
res_genes<-filter(ag_amrfinder,`Gene symbol` %in% ag_genes)
guuids$predict_res<-ifelse(guuids$guuid %in% res_genes$guuid,1,0)
cm<-caret::confusionMatrix(as.factor(guuids$predict_res),reference=as.factor(guuids$Gentamicin),positive="1")


## prevelance
650/(7803+650)

tn<-cm$table[1,1]
fn<-cm$table[1,2]
tp<-cm$table[2,2]
fp<-cm$table[2,1]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
#sens
prop.test(x=(tp),n=(tp+fn))
#spec
prop.test(x=(tn),n=(tn+fp))
#npv
prop.test(x=(tn),n=(tn+fn))
#ppv
prop.test(x=(tp),n=(tp+fp))
#Maj
prop.test(x=(fp),n=(tn+fp))
#Vmj
prop.test(x=(fn),n=(tp+fn))


############ done with stats


out=NULL


ag_res<-filter(guuids,Gentamicin==1)
ag_sens<-filter(guuids, Gentamicin==0)
ag_sens_running=NULL
ag_res_running=NULL
for(i in ag_genes){
  
  found<-filter(ag_amrfinder,`Gene symbol` == i)
  ag_res_<-filter(ag_res,guuid %in% found$guuid)
  ag_res_running<-rbind(ag_res_,ag_res_running)
  ag_res_running<-distinct(ag_res_running,guuid,.keep_all = T)
  explained<-nrow(ag_res_running)/nrow(ag_res)
  ag_sens_<-filter(ag_sens,guuid %in% found$guuid)
  ag_sens_running<-rbind(ag_sens_,ag_sens_running)
  ag_sens_running<-distinct(ag_sens_running,.keep_all = T)
  wrong<-nrow(ag_sens_running)/nrow(ag_sens)
  
  s_with_gene<-nrow(ag_sens_)/nrow(ag_sens)
  r_with_gene<-nrow(ag_res_)/nrow(ag_res)
  out<-rbind(out,data.frame(i,explained,wrong,s_with_gene,r_with_gene))
}


out$explained<-round((out$explained *100),1)
out$wrong<-round((out$wrong * 100),1)
out$s_with_gene<-round((out$s_with_gene * 100),1)
out$r_with_gene<-round((out$r_with_gene * 100),1)
out$index<-1:nrow(out)
out$index<-paste0(out$index,'   (', out$i,')')
out$index<-factor(out$index,levels=(out$index))

out2<-select(out,index,s_with_gene,r_with_gene)
out2<-pivot_longer(out2,names_to = "which",values_to = "pc",cols = 2:3)
out2$which<- case_when(out2$which =='r_with_gene' ~ "Resistant isolates with gene",
                       out2$which == 's_with_gene' ~ "Sensistive isolates with gene")
plot<-ggplot(out2) +
  geom_bar(aes(x=index,y=pc,fill=which),position="dodge",stat="identity", alpha=0.7) + 
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  scale_y_continuous(breaks=seq(0,100,10),limits = c(0,100)) +
  xlab("") +
  ylab("% isolates") +
  theme(axis.text.x = element_text(angle=270,hjust=0,vjust=0.2)) +
  labs(fill="") 


plot<-plot + new_scale_color()

Gentamicin<-plot+
  geom_point(data=out,aes(x=index,y=explained,color="Resistance explained")) +
  geom_path(data=out,aes(x=index,y=explained,color="Resistance explained"),group=explained) +
  geom_point(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified")) +
  geom_path(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified"),group=wrong) +
  scale_color_brewer(palette = "Set1")  +
  labs(color="") + ggtitle("Gentamicin")


#######ER coamox###########


all_guuids<-read_tsv('./data/guuids',col_names = c("guuid"))

res<-read_tsv('./data/global_pheno.tsv')

res<-filter(res,guuid %in% all_guuids$guuid)
all_guuids<-filter(all_guuids,guuid %in% res$guuid)

guuids<-left_join(all_guuids,res,by=c("guuid"))

amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv")
#amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv") %>% filter(`% Identity to reference sequence` ==100 & `% Coverage of reference sequence` ==100)


resfinder<-read_tsv('./data/resfinder_phenotypes.txt')
resfinder<-filter(resfinder,grepl('Clav',Phenotype))
resfinder$`Gene_accession no.`<-str_replace_all(resfinder$`Gene_accession no.`,'_.*','')
ag_res<-filter(guuids,Coamox==1)
ag_counts<-filter(amrfinder,grepl('BETA-LACTAM',Class)) %>% filter(`Gene symbol` %in% resfinder$`Gene_accession no.`) %>%  filter(guuid %in% ag_res$guuid) %>% group_by(`Gene symbol`) %>% count() %>% arrange(desc(n))
ag_genes<-unique(ag_counts$`Gene symbol`)
ag_genes<-ag_genes[!grepl('blaEC',ag_genes)]



ag_amrfinder<-filter(amrfinder,guuid %in% guuids$guuid)

###### here we calculate stats
res_genes<-filter(ag_amrfinder,`Gene symbol` %in% ag_genes)
guuids$predict_res<-ifelse(guuids$guuid %in% res_genes$guuid,1,0)
cm<-caret::confusionMatrix(as.factor(guuids$predict_res),reference=as.factor(guuids$Coamox),positive="1")

## prevelance
1621/(1621+3454)
tn<-cm$table[1,1]
fn<-cm$table[1,2]
tp<-cm$table[2,2]
fp<-cm$table[2,1]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
#sens
prop.test(x=(tp),n=(tp+fn))
#spec
prop.test(x=(tn),n=(tn+fp))
#npv
prop.test(x=(tn),n=(tn+fn))
#ppv
prop.test(x=(tp),n=(tp+fp))
#Maj
prop.test(x=(fp),n=(tn+fp))
#Vmj
prop.test(x=(fn),n=(tp+fn))
          
############ done with stats


out=NULL


ag_res<-filter(guuids,Coamox==1)
ag_sens<-filter(guuids, Coamox==0)
ag_sens_running=NULL
ag_res_running=NULL
tmp=NULL
for(i in ag_genes){
  
  found<-filter(ag_amrfinder,`Gene symbol` == i)
  ag_res_<-filter(ag_res,guuid %in% found$guuid)
  ag_res_running<-rbind(ag_res_,ag_res_running)
  ag_res_running<-distinct(ag_res_running,guuid,.keep_all = T)
  explained<-nrow(ag_res_running)/nrow(ag_res)
  ag_sens_<-filter(ag_sens,guuid %in% found$guuid)
  ag_sens_running<-rbind(ag_sens_,ag_sens_running)
  ag_sens_running<-distinct(ag_sens_running,.keep_all = T)
  wrong<-nrow(ag_sens_running)/nrow(ag_sens)
  
  s_with_gene<-nrow(ag_sens_)/nrow(ag_sens)
  r_with_gene<-nrow(ag_res_)/nrow(ag_res)
  out<-rbind(out,data.frame(i,explained,wrong,s_with_gene,r_with_gene))
  tmp<-rbind(tmp,data.frame(i,nrow(ag_sens),nrow(ag_sens_),nrow(ag_res),nrow(ag_res_)))
  
}


out$explained<-round((out$explained *100),1)
out$wrong<-round((out$wrong * 100),1)
out$s_with_gene<-round((out$s_with_gene * 100),1)
out$r_with_gene<-round((out$r_with_gene * 100),1)
out$index<-1:nrow(out)
out$index<-paste0(out$index,'   (', out$i,')')
out$index<-factor(out$index,levels=(out$index))

out2<-select(out,index,s_with_gene,r_with_gene)
out2<-pivot_longer(out2,names_to = "which",values_to = "pc",cols = 2:3)
out2$which<- case_when(out2$which =='r_with_gene' ~ "Resistant isolates with gene",
                       out2$which == 's_with_gene' ~ "Sensistive isolates with gene")
plot<-ggplot(out2) +
  geom_bar(aes(x=index,y=pc,fill=which),position="dodge",stat="identity", alpha=0.7) + 
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  scale_y_continuous(breaks=seq(0,100,10),limits = c(0,100)) +
  xlab("") +
  ylab("% isolates") +
  theme(axis.text.x = element_text(angle=270,hjust=0,vjust=0.2)) +
  labs(fill="") 


plot<-plot + new_scale_color()

Coamox<-plot+
  geom_point(data=out,aes(x=index,y=explained,color="Resistance explained")) +
  geom_path(data=out,aes(x=index,y=explained,color="Resistance explained"),group=explained) +
  geom_point(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified")) +
  geom_path(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified"),group=wrong) +
  scale_color_brewer(palette = "Set1")  +
  labs(color="") + ggtitle("Co-amoxiclav")

#sensitivty gap
100-max(out$explained)

##########ER Amp ###############

all_guuids<-read_tsv('./data/guuids',col_names = c("guuid"))

res<-read_tsv('./data/global_pheno.tsv')

res<-filter(res,guuid %in% all_guuids$guuid)
all_guuids<-filter(all_guuids,guuid %in% res$guuid)

guuids<-left_join(all_guuids,res,by=c("guuid"))

amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv")
#amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv") %>% filter(`% Identity to reference sequence` ==100 & `% Coverage of reference sequence` ==100)


ag_res<-filter(guuids,Ampicillin==1)
ag_counts<-filter(amrfinder,Class=='BETA-LACTAM') %>% filter(guuid %in% ag_res$guuid) %>% group_by(`Gene symbol`) %>% count() %>% arrange(desc(n))
ag_genes<-unique(ag_counts$`Gene symbol`)
ag_genes<-ag_genes[!grepl('blaEC',ag_genes)]

ag_amrfinder<-filter(amrfinder,guuid %in% guuids$guuid)

###### here we calculate stats
res_genes<-filter(ag_amrfinder,`Gene symbol` %in% ag_genes)
guuids$predict_res<-ifelse(guuids$guuid %in% res_genes$guuid,1,0)
cm<-caret::confusionMatrix(as.factor(guuids$predict_res),reference=as.factor(guuids$Ampicillin),positive="1")

## prevelance
3842/(3842+3639)

tn<-cm$table[1,1]
fn<-cm$table[1,2]
tp<-cm$table[2,2]
fp<-cm$table[2,1]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
#sens
prop.test(x=(tp),n=(tp+fn))
#spec
prop.test(x=(tn),n=(tn+fp))
#npv
prop.test(x=(tn),n=(tn+fn))
#ppv
prop.test(x=(tp),n=(tp+fp))
#Maj
prop.test(x=(fp),n=(tn+fp))
#Vmj
prop.test(x=(fn),n=(tp+fn))
          


############ done with stats


out=NULL


ag_res<-filter(guuids,Ampicillin==1)
ag_sens<-filter(guuids, Ampicillin==0)
ag_sens_running=NULL
ag_res_running=NULL
tmp=NULL
for(i in ag_genes){
  
  found<-filter(ag_amrfinder,`Gene symbol` == i)
  ag_res_<-filter(ag_res,guuid %in% found$guuid)
  ag_res_running<-rbind(ag_res_,ag_res_running)
  ag_res_running<-distinct(ag_res_running,guuid,.keep_all = T)
  explained<-nrow(ag_res_running)/nrow(ag_res)
  ag_sens_<-filter(ag_sens,guuid %in% found$guuid)
  ag_sens_running<-rbind(ag_sens_,ag_sens_running)
  ag_sens_running<-distinct(ag_sens_running,.keep_all = T)
  wrong<-nrow(ag_sens_running)/nrow(ag_sens)
  
  s_with_gene<-nrow(ag_sens_)/nrow(ag_sens)
  r_with_gene<-nrow(ag_res_)/nrow(ag_res)
  out<-rbind(out,data.frame(i,explained,wrong,s_with_gene,r_with_gene))
  tmp<-rbind(tmp,data.frame(i,nrow(ag_sens),nrow(ag_sens_),nrow(ag_res),nrow(ag_res_)))
}


out$explained<-round((out$explained *100),1)
out$wrong<-round((out$wrong * 100),1)
out$s_with_gene<-round((out$s_with_gene * 100),1)
out$r_with_gene<-round((out$r_with_gene * 100),1)
out$index<-1:nrow(out)
out$index<-paste0(out$index,'   (', out$i,')')
out$index<-factor(out$index,levels=(out$index))

out2<-select(out,index,s_with_gene,r_with_gene)
out2<-pivot_longer(out2,names_to = "which",values_to = "pc",cols = 2:3)
out2$which<- case_when(out2$which =='r_with_gene' ~ "Resistant isolates with gene",
                       out2$which == 's_with_gene' ~ "Sensistive isolates with gene")
plot<-ggplot(out2) +
  geom_bar(aes(x=index,y=pc,fill=which),position="dodge",stat="identity", alpha=0.7) + 
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  scale_y_continuous(breaks=seq(0,100,10),limits = c(0,100)) +
  xlab("") +
  ylab("% isolates") +
  theme(axis.text.x = element_text(angle=270,hjust=0,vjust=0.2)) +
  labs(fill="") 

  
plot<-plot + new_scale_color()
  
Ampicillin<-plot+
  geom_point(data=out,aes(x=index,y=explained,color="Resistance explained")) +
  geom_path(data=out,aes(x=index,y=explained,color="Resistance explained"),group=explained) +
  geom_point(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified")) +
  geom_path(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified"),group=wrong) +
  scale_color_brewer(palette = "Set1")  +
  labs(color="") +
  ggtitle("Ampicillin")


##############ER Cipro###############

all_guuids<-read_tsv('./data/guuids',col_names = c("guuid"))

res<-read_tsv('./data/global_pheno.tsv')

res<-filter(res,guuid %in% all_guuids$guuid)
all_guuids<-filter(all_guuids,guuid %in% res$guuid)

guuids<-left_join(all_guuids,res,by=c("guuid"))

amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv")
#amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv") %>% filter(`% Identity to reference sequence` ==100 & `% Coverage of reference sequence` ==100)


ag_res<-filter(guuids,Ciprofloxacin==1)
ag_counts<-filter(amrfinder,Class=='QUINOLONE') %>% filter(guuid %in% ag_res$guuid) %>% group_by(`Gene symbol`) %>% count() %>% arrange(desc(n))
ag_genes<-unique(ag_counts$`Gene symbol`)
ag_genes<-ag_genes[!grepl('blaEC',ag_genes)]

ag_amrfinder<-filter(amrfinder,guuid %in% guuids$guuid)

###### here we calculate stats
res_genes<-filter(ag_amrfinder,`Gene symbol` %in% ag_genes)
guuids$predict_res<-ifelse(guuids$guuid %in% res_genes$guuid,1,0)
cm<-caret::confusionMatrix(as.factor(guuids$predict_res),reference=as.factor(guuids$Ciprofloxacin),positive="1")

## prevelance
1165/(1165+7331)

tn<-cm$table[1,1]
fn<-cm$table[1,2]
tp<-cm$table[2,2]
fp<-cm$table[2,1]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
#sens
prop.test(x=(tp),n=(tp+fn))
#spec
prop.test(x=(tn),n=(tn+fp))
#npv
prop.test(x=(tn),n=(tn+fn))
#ppv
prop.test(x=(tp),n=(tp+fp))
#Maj
prop.test(x=(fp),n=(tn+fp))
#Vmj
prop.test(x=(fn),n=(tp+fn))


############ done with stats


out=NULL


ag_res<-filter(guuids,Ciprofloxacin==1)
ag_sens<-filter(guuids, Ciprofloxacin==0)
ag_sens_running=NULL
ag_res_running=NULL
for(i in ag_genes){
  
  found<-filter(ag_amrfinder,`Gene symbol` == i)
  ag_res_<-filter(ag_res,guuid %in% found$guuid)
  ag_res_running<-rbind(ag_res_,ag_res_running)
  ag_res_running<-distinct(ag_res_running,guuid,.keep_all = T)
  explained<-nrow(ag_res_running)/nrow(ag_res)
  ag_sens_<-filter(ag_sens,guuid %in% found$guuid)
  ag_sens_running<-rbind(ag_sens_,ag_sens_running)
  ag_sens_running<-distinct(ag_sens_running,.keep_all = T)
  wrong<-nrow(ag_sens_running)/nrow(ag_sens)
  
  s_with_gene<-nrow(ag_sens_)/nrow(ag_sens)
  r_with_gene<-nrow(ag_res_)/nrow(ag_res)
  out<-rbind(out,data.frame(i,explained,wrong,s_with_gene,r_with_gene))
}


out$explained<-round((out$explained *100),1)
out$wrong<-round((out$wrong * 100),1)
out$s_with_gene<-round((out$s_with_gene * 100),1)
out$r_with_gene<-round((out$r_with_gene * 100),1)
out$index<-1:nrow(out)
out$index<-paste0(out$index,'   (', out$i,')')
out$index<-factor(out$index,levels=(out$index))

out2<-select(out,index,s_with_gene,r_with_gene)
out2<-pivot_longer(out2,names_to = "which",values_to = "pc",cols = 2:3)
out2$which<- case_when(out2$which =='r_with_gene' ~ "Resistant isolates with gene",
                       out2$which == 's_with_gene' ~ "Sensistive isolates with gene")
plot<-ggplot(out2) +
  geom_bar(aes(x=index,y=pc,fill=which),position="dodge",stat="identity", alpha=0.7) + 
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  scale_y_continuous(breaks=seq(0,100,10),limits = c(0,100)) +
  xlab("") +
  ylab("% isolates") +
  theme(axis.text.x = element_text(angle=270,hjust=0,vjust=0.2)) +
  labs(fill="") 


plot<-plot + new_scale_color()

Ciprofloxacin<-plot+
  geom_point(data=out,aes(x=index,y=explained,color="Resistance explained")) +
  geom_path(data=out,aes(x=index,y=explained,color="Resistance explained"),group=explained) +
  geom_point(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified")) +
  geom_path(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified"),group=wrong) +
  scale_color_brewer(palette = "Set1")  +
  labs(color="") + 
  ggtitle("Ciprofloxacin")

#sensitivity gap
100-max(out$explained)

#######################ER Ceftriaxone##############


all_guuids<-read_tsv('./data/guuids',col_names = c("guuid"))

res<-read_tsv('./data/global_pheno.tsv')

res<-filter(res,guuid %in% all_guuids$guuid)
all_guuids<-filter(all_guuids,guuid %in% res$guuid)

guuids<-left_join(all_guuids,res,by=c("guuid"))

amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv")
#amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv") %>% filter(`% Identity to reference sequence` ==100 & `% Coverage of reference sequence` ==100)


ag_res<-filter(guuids,Ceftriaxone==1)
ag_counts<-filter(amrfinder,Subclass=='CEPHALOSPORIN') %>% filter(guuid %in% ag_res$guuid) %>% group_by(`Gene symbol`) %>% count() %>% arrange(desc(n))
ag_genes<-unique(ag_counts$`Gene symbol`)
ag_genes<-ag_genes[!grepl('blaEC',ag_genes)]

ag_amrfinder<-filter(amrfinder,guuid %in% guuids$guuid)

###### here we calculate stats
res_genes<-filter(ag_amrfinder,`Gene symbol` %in% ag_genes)
guuids$predict_res<-ifelse(guuids$guuid %in% res_genes$guuid,1,0)
cm<-caret::confusionMatrix(as.factor(guuids$predict_res),reference=as.factor(guuids$Ceftriaxone),positive="1")

## prevelance
270/(270+3122)

tn<-cm$table[1,1]
fn<-cm$table[1,2]
tp<-cm$table[2,2]
fp<-cm$table[2,1]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
#sens
prop.test(x=(tp),n=(tp+fn))
#spec
prop.test(x=(tn),n=(tn+fp))
#npv
prop.test(x=(tn),n=(tn+fn))
#ppv
prop.test(x=(tp),n=(tp+fp))
#Maj
prop.test(x=(fp),n=(tn+fp))
#Vmj
prop.test(x=(fn),n=(tp+fn))


############ done with stats


out=NULL


ag_res<-filter(guuids,Ceftriaxone==1)
ag_sens<-filter(guuids, Ceftriaxone==0)
ag_sens_running=NULL
ag_res_running=NULL
for(i in ag_genes){
  
  found<-filter(ag_amrfinder,`Gene symbol` == i)
  ag_res_<-filter(ag_res,guuid %in% found$guuid)
  ag_res_running<-rbind(ag_res_,ag_res_running)
  ag_res_running<-distinct(ag_res_running,guuid,.keep_all = T)
  explained<-nrow(ag_res_running)/nrow(ag_res)
  ag_sens_<-filter(ag_sens,guuid %in% found$guuid)
  ag_sens_running<-rbind(ag_sens_,ag_sens_running)
  ag_sens_running<-distinct(ag_sens_running,.keep_all = T)
  wrong<-nrow(ag_sens_running)/nrow(ag_sens)
  
  s_with_gene<-nrow(ag_sens_)/nrow(ag_sens)
  r_with_gene<-nrow(ag_res_)/nrow(ag_res)
  out<-rbind(out,data.frame(i,explained,wrong,s_with_gene,r_with_gene))
}


out$explained<-round((out$explained *100),1)
out$wrong<-round((out$wrong * 100),1)
out$s_with_gene<-round((out$s_with_gene * 100),1)
out$r_with_gene<-round((out$r_with_gene * 100),1)
out$index<-1:nrow(out)
out$index<-paste0(out$index,'   (', out$i,')')
out$index<-factor(out$index,levels=(out$index))

out2<-select(out,index,s_with_gene,r_with_gene)
out2<-pivot_longer(out2,names_to = "which",values_to = "pc",cols = 2:3)
out2$which<- case_when(out2$which =='r_with_gene' ~ "Resistant isolates with gene",
                       out2$which == 's_with_gene' ~ "Sensistive isolates with gene")
plot<-ggplot(out2) +
  geom_bar(aes(x=index,y=pc,fill=which),position="dodge",stat="identity", alpha=0.7) + 
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  scale_y_continuous(breaks=seq(0,100,10),limits = c(0,100)) +
  xlab("") +
  ylab("% isolates") +
  theme(axis.text.x = element_text(angle=270,hjust=0,vjust=0.2)) +
  labs(fill="") 


plot<-plot + new_scale_color()

Ceftriaxone<-plot+
  geom_point(data=out,aes(x=index,y=explained,color="Resistance explained")) +
  geom_path(data=out,aes(x=index,y=explained,color="Resistance explained"),group=explained) +
  geom_point(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified")) +
  geom_path(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified"),group=wrong) +
  scale_color_brewer(palette = "Set1")  +
  labs(color="") +
  ggtitle("Ceftriaxone")

##################ER Fosfomycin#######


all_guuids<-read_tsv('./data/guuids',col_names = c("guuid"))

res<-read_tsv('./data/global_pheno.tsv')

res<-filter(res,guuid %in% all_guuids$guuid)
all_guuids<-filter(all_guuids,guuid %in% res$guuid)

guuids<-left_join(all_guuids,res,by=c("guuid"))

amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv")
#amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv") %>% filter(`% Identity to reference sequence` ==100 & `% Coverage of reference sequence` ==100)


ag_res<-filter(guuids,Fosfomycin==1)
ag_counts<-filter(amrfinder,grepl('FOSFOMYCIN',Class)) %>% filter(guuid %in% ag_res$guuid) %>% group_by(`Gene symbol`) %>% count() %>% arrange(desc(n))
ag_genes<-unique(ag_counts$`Gene symbol`)
ag_genes<-ag_genes[!grepl('blaEC',ag_genes)]

ag_amrfinder<-filter(amrfinder,guuid %in% guuids$guuid)

###### here we calculate stats
res_genes<-filter(ag_amrfinder,`Gene symbol` %in% ag_genes)
guuids$predict_res<-ifelse(guuids$guuid %in% res_genes$guuid,1,0)
cm<-caret::confusionMatrix(as.factor(guuids$predict_res),reference=as.factor(guuids$Fosfomycin),positive="1")

## prevelance
(15/(15+1708))*100

tn<-cm$table[1,1]
fn<-cm$table[1,2]
tp<-cm$table[2,2]
fp<-cm$table[2,1]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
#sens
prop.test(x=(tp),n=(tp+fn))
#spec
prop.test(x=(tn),n=(tn+fp))
#npv
prop.test(x=(tn),n=(tn+fn))
#ppv
prop.test(x=(tp),n=(tp+fp))
#Maj
prop.test(x=(fp),n=(tn+fp))
#Vmj
prop.test(x=(fn),n=(tp+fn))

############ done with stats


out=NULL


ag_res<-filter(guuids,Fosfomycin==1)
ag_sens<-filter(guuids, Fosfomycin==0)
ag_sens_running=NULL
ag_res_running=NULL
for(i in ag_genes){
  
  found<-filter(ag_amrfinder,`Gene symbol` == i)
  ag_res_<-filter(ag_res,guuid %in% found$guuid)
  ag_res_running<-rbind(ag_res_,ag_res_running)
  ag_res_running<-distinct(ag_res_running,guuid,.keep_all = T)
  explained<-nrow(ag_res_running)/nrow(ag_res)
  ag_sens_<-filter(ag_sens,guuid %in% found$guuid)
  ag_sens_running<-rbind(ag_sens_,ag_sens_running)
  ag_sens_running<-distinct(ag_sens_running,.keep_all = T)
  wrong<-nrow(ag_sens_running)/nrow(ag_sens)
  
  s_with_gene<-nrow(ag_sens_)/nrow(ag_sens)
  r_with_gene<-nrow(ag_res_)/nrow(ag_res)
  out<-rbind(out,data.frame(i,explained,wrong,s_with_gene,r_with_gene))
}


out$explained<-round((out$explained *100),1)
out$wrong<-round((out$wrong * 100),1)
out$s_with_gene<-round((out$s_with_gene * 100),1)
out$r_with_gene<-round((out$r_with_gene * 100),1)
out$index<-1:nrow(out)
out$index<-paste0(out$index,'   (', out$i,')')
out$index<-factor(out$index,levels=(out$index))

out2<-select(out,index,s_with_gene,r_with_gene)
out2<-pivot_longer(out2,names_to = "which",values_to = "pc",cols = 2:3)
out2$which<- case_when(out2$which =='r_with_gene' ~ "Resistant isolates with gene",
                       out2$which == 's_with_gene' ~ "Sensistive isolates with gene")
plot<-ggplot(out2) +
  geom_bar(aes(x=index,y=pc,fill=which),position="dodge",stat="identity", alpha=0.7) + 
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  scale_y_continuous(breaks=seq(0,100,10),limits = c(0,100)) +
  xlab("") +
  ylab("% isolates") +
  theme(axis.text.x = element_text(angle=270,hjust=0,vjust=0.2)) +
  labs(fill="") 


plot<-plot + new_scale_color()

Fosfomycin<-plot+
  geom_point(data=out,aes(x=index,y=explained,color="Resistance explained")) +
  geom_path(data=out,aes(x=index,y=explained,color="Resistance explained"),group=explained) +
  geom_point(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified")) +
  geom_path(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified"),group=wrong) +
  scale_color_brewer(palette = "Set1")  +
  labs(color="") +
  ggtitle("Fosfomycin")

##########ER Trimethoprim########



all_guuids<-read_tsv('./data/guuids',col_names = c("guuid"))

res<-read_tsv('./data/global_pheno.tsv')

res<-filter(res,guuid %in% all_guuids$guuid)
all_guuids<-filter(all_guuids,guuid %in% res$guuid)

guuids<-left_join(all_guuids,res,by=c("guuid"))

amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv")
#amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv") %>% filter(`% Identity to reference sequence` ==100 & `% Coverage of reference sequence` ==100)


ag_res<-filter(guuids,Trimethoprim==1)
ag_counts<-filter(amrfinder,grepl('TRIMETHOPRIM',Class)) %>% filter(guuid %in% ag_res$guuid) %>% group_by(`Gene symbol`) %>% count() %>% arrange(desc(n))
ag_genes<-unique(ag_counts$`Gene symbol`)
ag_genes<-ag_genes[!grepl('blaEC',ag_genes)]

ag_amrfinder<-filter(amrfinder,guuid %in% guuids$guuid)

###### here we calculate stats
res_genes<-filter(ag_amrfinder,`Gene symbol` %in% ag_genes)
guuids$predict_res<-ifelse(guuids$guuid %in% res_genes$guuid,1,0)
cm<-caret::confusionMatrix(as.factor(guuids$predict_res),reference=as.factor(guuids$Trimethoprim),positive="1")

## prevelance
1515/(1515+2592)

tn<-cm$table[1,1]
fn<-cm$table[1,2]
tp<-cm$table[2,2]
fp<-cm$table[2,1]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
#sens
prop.test(x=(tp),n=(tp+fn))
#spec
prop.test(x=(tn),n=(tn+fp))
#npv
prop.test(x=(tn),n=(tn+fn))
#ppv
prop.test(x=(tp),n=(tp+fp))
#Maj
prop.test(x=(fp),n=(tn+fp))
#Vmj
prop.test(x=(fn),n=(tp+fn))


############ done with stats

out=NULL


ag_res<-filter(guuids,Trimethoprim==1)
ag_sens<-filter(guuids, Trimethoprim==0)
ag_sens_running=NULL
ag_res_running=NULL
for(i in ag_genes){
  
  found<-filter(ag_amrfinder,`Gene symbol` == i)
  ag_res_<-filter(ag_res,guuid %in% found$guuid)
  ag_res_running<-rbind(ag_res_,ag_res_running)
  ag_res_running<-distinct(ag_res_running,guuid,.keep_all = T)
  explained<-nrow(ag_res_running)/nrow(ag_res)
  ag_sens_<-filter(ag_sens,guuid %in% found$guuid)
  ag_sens_running<-rbind(ag_sens_,ag_sens_running)
  ag_sens_running<-distinct(ag_sens_running,.keep_all = T)
  wrong<-nrow(ag_sens_running)/nrow(ag_sens)
  
  s_with_gene<-nrow(ag_sens_)/nrow(ag_sens)
  r_with_gene<-nrow(ag_res_)/nrow(ag_res)
  out<-rbind(out,data.frame(i,explained,wrong,s_with_gene,r_with_gene))
}


out$explained<-round((out$explained *100),1)
out$wrong<-round((out$wrong * 100),1)
out$s_with_gene<-round((out$s_with_gene * 100),1)
out$r_with_gene<-round((out$r_with_gene * 100),1)
out$index<-1:nrow(out)
out$index<-paste0(out$index,'   (', out$i,')')
out$index<-factor(out$index,levels=(out$index))

out2<-select(out,index,s_with_gene,r_with_gene)
out2<-pivot_longer(out2,names_to = "which",values_to = "pc",cols = 2:3)
out2$which<- case_when(out2$which =='r_with_gene' ~ "Resistant isolates with gene",
                       out2$which == 's_with_gene' ~ "Sensistive isolates with gene")
plot<-ggplot(out2) +
  geom_bar(aes(x=index,y=pc,fill=which),position="dodge",stat="identity", alpha=0.7) + 
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  scale_y_continuous(breaks=seq(0,100,10),limits = c(0,100)) +
  xlab("") +
  ylab("% isolates") +
  theme(axis.text.x = element_text(angle=270,hjust=0,vjust=0.2)) +
  labs(fill="") 


plot<-plot + new_scale_color()

Trimethoprim<-plot+
  geom_point(data=out,aes(x=index,y=explained,color="Resistance explained")) +
  geom_path(data=out,aes(x=index,y=explained,color="Resistance explained"),group=explained) +
  geom_point(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified")) +
  geom_path(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified"),group=wrong) +
  scale_color_brewer(palette = "Set1")  +
  labs(color="") +
  ggtitle("Trimethoprim")

##########ER AMIKACIN



all_guuids<-read_tsv('./data/guuids',col_names = c("guuid"))

res<-read_tsv('./data/global_pheno.tsv')

res<-filter(res,guuid %in% all_guuids$guuid)
all_guuids<-filter(all_guuids,guuid %in% res$guuid)

guuids<-left_join(all_guuids,res,by=c("guuid"))

amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv")

ag_res<-filter(guuids,Amikacin==1)
ag_counts<-filter(amrfinder,grepl('AMIKACIN',Subclass)) %>% group_by(`Gene symbol`) %>% count() %>% arrange(desc(n))
ag_genes<-unique(ag_counts$`Gene symbol`)
ag_genes<-ag_genes[!grepl('blaEC',ag_genes)]

ag_amrfinder<-filter(amrfinder,guuid %in% guuids$guuid)
out=NULL


ag_res<-filter(guuids,Amikacin==1)
ag_sens<-filter(guuids, Amikacin==0)
ag_sens_running=NULL
ag_res_running=NULL
for(i in ag_genes){
  
  found<-filter(ag_amrfinder,`Gene symbol` == i)
  ag_res_<-filter(ag_res,guuid %in% found$guuid)
  ag_res_running<-rbind(ag_res_,ag_res_running)
  ag_res_running<-distinct(ag_res_running,guuid,.keep_all = T)
  explained<-nrow(ag_res_running)/nrow(ag_res)
  ag_sens_<-filter(ag_sens,guuid %in% found$guuid)
  ag_sens_running<-rbind(ag_sens_,ag_sens_running)
  ag_sens_running<-distinct(ag_sens_running,.keep_all = T)
  wrong<-nrow(ag_sens_running)/nrow(ag_sens)
  
  s_with_gene<-nrow(ag_sens_)/nrow(ag_sens)
  r_with_gene<-nrow(ag_res_)/nrow(ag_res)
  out<-rbind(out,data.frame(i,explained,wrong,s_with_gene,r_with_gene))
}


out$explained<-round((out$explained *100),1)
out$wrong<-round((out$wrong * 100),1)
out$s_with_gene<-round((out$s_with_gene * 100),1)
out$r_with_gene<-round((out$r_with_gene * 100),1)
out$index<-1:nrow(out)
out$index<-paste0(out$index,'   (', out$i,')')
out$index<-factor(out$index,levels=(out$index))

out2<-select(out,index,s_with_gene,r_with_gene)
out2<-pivot_longer(out2,names_to = "which",values_to = "pc",cols = 2:3)
out2$which<- case_when(out2$which =='r_with_gene' ~ "Resistant isolates with gene",
                       out2$which == 's_with_gene' ~ "Sensistive isolates with gene")
plot<-ggplot(out2) +
  geom_bar(aes(x=index,y=pc,fill=which),position="dodge",stat="identity", alpha=0.7) + 
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  scale_y_continuous(breaks=seq(0,100,10),limits = c(0,100)) +
  xlab("") +
  ylab("% isolates") +
  theme(axis.text.x = element_text(angle=270,hjust=0,vjust=0.2)) +
  labs(fill="") 


plot<-plot + new_scale_color()

Amikacin<-plot+
  geom_point(data=out,aes(x=index,y=explained,color="Resistance explained")) +
  geom_path(data=out,aes(x=index,y=explained,color="Resistance explained"),group=explained) +
  geom_point(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified")) +
  geom_path(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified"),group=wrong) +
  scale_color_brewer(palette = "Set1")  +
  labs(color="")

################## ER Piptaz ####################



all_guuids<-read_tsv('./data/guuids',col_names = c("guuid"))

res<-read_tsv('./data/global_pheno.tsv')

res<-filter(res,guuid %in% all_guuids$guuid)
all_guuids<-filter(all_guuids,guuid %in% res$guuid)

guuids<-left_join(all_guuids,res,by=c("guuid"))

amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv")
#amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv") %>% filter(`% Identity to reference sequence` ==100 & `% Coverage of reference sequence` ==100)


resfinder<-read_tsv('./data/resfinder_phenotypes.txt')
resfinder<-filter(resfinder,grepl('Piperacillin',Phenotype) & grepl('Tazo',Phenotype))
resfinder<-distinct(resfinder,`Gene_accession no.`)
resfinder$`Gene_accession no.`<-str_replace_all(resfinder$`Gene_accession no.`,'_.*','')
ag_res<-filter(guuids,Piptaz==1)
ag_counts<-filter(amrfinder,grepl('BETA-LACTAM',Class)) %>% filter(guuid %in% ag_res$guuid) %>% filter(`Gene symbol` %in% resfinder$`Gene_accession no.`) %>% group_by(`Gene symbol`) %>% count() %>% arrange(desc(n))
ag_genes<-unique(ag_counts$`Gene symbol`)
ag_genes<-ag_genes[!grepl('blaEC',ag_genes)]

ag_amrfinder<-filter(amrfinder,guuid %in% guuids$guuid)

###### here we calculate stats
res_genes<-filter(ag_amrfinder,`Gene symbol` %in% ag_genes)
guuids$predict_res<-ifelse(guuids$guuid %in% res_genes$guuid,1,0)
cm<-caret::confusionMatrix(as.factor(guuids$predict_res),reference=as.factor(guuids$Piptaz),positive="1")

## prevelance
329/(329+7832)
tn<-cm$table[1,1]
fn<-cm$table[1,2]
tp<-cm$table[2,2]
fp<-cm$table[2,1]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
#sens
prop.test(x=(tp),n=(tp+fn))
#spec
prop.test(x=(tn),n=(tn+fp))
#npv
prop.test(x=(tn),n=(tn+fn))
#ppv
prop.test(x=(tp),n=(tp+fp))
#Maj
prop.test(x=(fp),n=(tn+fp))
#Vmj
prop.test(x=(fn),n=(tp+fn))

############ done with stats


out=NULL


ag_res<-filter(guuids,Piptaz==1)
ag_sens<-filter(guuids, Piptaz==0)
ag_sens_running=NULL
ag_res_running=NULL
tmp=NULL
for(i in ag_genes){
  
  found<-filter(ag_amrfinder,`Gene symbol` == i)
  ag_res_<-filter(ag_res,guuid %in% found$guuid)
  ag_res_running<-rbind(ag_res_,ag_res_running)
  ag_res_running<-distinct(ag_res_running,guuid,.keep_all = T)
  explained<-nrow(ag_res_running)/nrow(ag_res)
  ag_sens_<-filter(ag_sens,guuid %in% found$guuid)
  ag_sens_running<-rbind(ag_sens_,ag_sens_running)
  ag_sens_running<-distinct(ag_sens_running,.keep_all = T)
  wrong<-nrow(ag_sens_running)/nrow(ag_sens)
  
  s_with_gene<-nrow(ag_sens_)/nrow(ag_sens)
  r_with_gene<-nrow(ag_res_)/nrow(ag_res)
  out<-rbind(out,data.frame(i,explained,wrong,s_with_gene,r_with_gene))
  tmp<-rbind(tmp,data.frame(i,nrow(ag_sens),nrow(ag_sens_),nrow(ag_res),nrow(ag_res_)))
  
}


out$explained<-round((out$explained *100),1)
out$wrong<-round((out$wrong * 100),1)
out$s_with_gene<-round((out$s_with_gene * 100),1)
out$r_with_gene<-round((out$r_with_gene * 100),1)
out$index<-1:nrow(out)
out$index<-paste0(out$index,'   (', out$i,')')
out$index<-factor(out$index,levels=(out$index))

out2<-select(out,index,s_with_gene,r_with_gene)
out2<-pivot_longer(out2,names_to = "which",values_to = "pc",cols = 2:3)
out2$which<- case_when(out2$which =='r_with_gene' ~ "Resistant isolates with gene",
                       out2$which == 's_with_gene' ~ "Sensistive isolates with gene")
plot<-ggplot(out2) +
  geom_bar(aes(x=index,y=pc,fill=which),position="dodge",stat="identity", alpha=0.7) + 
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  scale_y_continuous(breaks=seq(0,100,10),limits = c(0,100)) +
  xlab("") +
  ylab("% isolates") +
  theme(axis.text.x = element_text(angle=270,hjust=0,vjust=0.2)) +
  labs(fill="") 


plot<-plot + new_scale_color()

Piptaz<-plot+
  geom_point(data=out,aes(x=index,y=explained,color="Resistance explained")) +
  geom_path(data=out,aes(x=index,y=explained,color="Resistance explained"),group=explained) +
  geom_point(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified")) +
  geom_path(data=out,aes(x=index,y=wrong,color="Sensitive isolates wrongly classified"),group=wrong) +
  scale_color_brewer(palette = "Set1")  +
  labs(color="") +
  ggtitle("Piperacillin-Tazobactam")

Fig1<-(Ampicillin + Coamox) / (Piptaz + Ceftriaxone+ Ciprofloxacin) +(Gentamicin + Trimethoprim + Fosfomycin) +
  plot_layout(guides="collect")



###################app1 amox stat##################

all_guuids<-read_tsv('./data/guuids',col_names = c("guuid"))

res<-read_tsv('./data/global_pheno.tsv')

res<-filter(res,guuid %in% all_guuids$guuid)
all_guuids<-filter(all_guuids,guuid %in% res$guuid)

guuids<-left_join(all_guuids,res,by=c("guuid"))
guuids$amp_amox<-case_when(guuids$Ampicillin ==1 ~ 1,
                           guuids$Ampicillin ==0 ~ 0,
                           guuids$Amoxicillin ==1 ~ 1,
                           guuids$Amoxicillin ==0 ~ 0,
                           TRUE ~ 9)


guuids$amp_amox[guuids$amp_amox ==9]<-NA
amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv")

ag_res<-filter(guuids,amp_amox==1)
ag_counts<-filter(amrfinder,Class=='BETA-LACTAM') %>% filter(guuid %in% ag_res$guuid) %>% group_by(`Gene symbol`) %>% count() %>% arrange(desc(n))
ag_genes<-unique(ag_counts$`Gene symbol`)
ag_genes<-ag_genes[!grepl('blaEC',ag_genes)]

ag_amrfinder<-filter(amrfinder,guuid %in% guuids$guuid)
###### here we calculate stats
res_genes<-filter(ag_amrfinder,`Gene symbol` %in% ag_genes)
guuids$predict_res<-ifelse(guuids$guuid %in% res_genes$guuid,1,0)
cm<-caret::confusionMatrix(as.factor(guuids$predict_res),reference=as.factor(guuids$amp_amox),positive="1")

## prevelance
650/(7803+650)

concordance<-((cm$table[2,2] + cm$table[1,1])/(cm$table[2,2] + cm$table[1,1] + cm$table[1,2] + cm$table[2,1]))

conc.fun<-function(data,idx){
  g<-data[idx,]
  cmb<-caret::confusionMatrix(as.factor(g$predict_res),reference=as.factor(g$Gentamicin),positive="1")
  ((cmb$table[2,2] + cmb$table[1,1])/(cmb$table[2,2] + cmb$table[1,1] + cmb$table[1,2] + cmb$table[2,1]))*100
}

conc_boot<-boot::boot(guuids,conc.fun,R=1000)
boot::boot.ci(boot.out=conc_boot,type="perc")

sens<-(cm$table[2,2]/(cm$table[2,2] + cm$table[1,2]))*100
spec<-(cm$table[1,1]/(cm$table[1,1] + cm$table[2,1]))*100

sens.fun<-function(data,idx){
  g<-data[idx,]
  cmb<-caret::confusionMatrix(as.factor(g$predict_res),reference=as.factor(g$Gentamicin),positive="1")
  (cmb$table[2,2]/(cmb$table[2,2]+cmb$table[1,2]))*100
}

spec.fun<-function(data,idx){
  g<-data[idx,]
  cmb<-caret::confusionMatrix(as.factor(g$predict_res),reference=as.factor(g$Gentamicin),positive="1")
  (cmb$table[1,1]/(cmb$table[1,1]+cmb$table[2,1]))*100
}

sens_boot<-boot::boot(guuids,sens.fun,R=1000)
boot::boot.ci(boot.out=sens_boot,type="perc")

spec_boot<-boot::boot(guuids,spec.fun,R=1000)
boot::boot.ci(boot.out=spec_boot,type="perc")

npv<-(cm$table[1,1]/(cm$table[1,1] + cm$table[2,1]))*100
ppv<-(cm$table[2,2]/(cm$table[2,2] + cm$table[1,2]))*100

npv.fun<-function(data,idx){
  g<-data[idx,]
  cmb<-caret::confusionMatrix(as.factor(g$predict_res),reference=as.factor(g$Gentamicin),positive="1")
  (cmb$table[1,1]/(cmb$table[1,1] + cmb$table[2,1]))*100
}

ppv.fun<-function(data,idx){
  g<-data[idx,]
  cmb<-caret::confusionMatrix(as.factor(g$predict_res),reference=as.factor(g$Gentamicin),positive="1")
  (cmb$table[2,2]/(cmb$table[2,2] + cmb$table[1,2]))*100
}

npv_boot<-boot::boot(guuids,npv.fun,R=1000)
ppv_boot<-boot::boot(guuids,ppv.fun,R=1000)

boot::boot.ci(boot.out=npv_boot,type="perc")
boot::boot.ci(boot.out=ppv_boot,type="perc")

ME<-(cm$table[2,1]/(cm$table[1,1]+cm$table[2,1]))*100

me.fun<-function(data,idx){
  g<-data[idx,]
  cmb<-caret::confusionMatrix(as.factor(g$predict_res),reference=as.factor(g$Gentamicin),positive="1")
  (cmb$table[2,1]/(cmb$table[1,1]+cmb$table[2,1]))*100
}
me_boot<-boot::boot(guuids,me.fun,R=1000)
boot::boot.ci(boot.out=me_boot,type="perc")


VME<-(cm$table[1,2]/(cm$table[2,2] + cm$table[1,2]))*100

vme.fun<-function(data,idx){
  g<-data[idx,]
  cmb<-caret::confusionMatrix(as.factor(g$predict_res),reference=as.factor(g$Gentamicin),positive="1")
  (cmb$table[1,2]/(cmb$table[2,2]+cmb$table[1,2]))*100
}

vme_boot<-boot::boot(guuids,vme.fun,R=1000)
boot::boot.ci(boot.out=vme_boot,type="perc")

############ done with stats


################## Piptaz with blaTEM-1 ####################



all_guuids<-read_tsv('./data/guuids',col_names = c("guuid"))

res<-read_tsv('./data/global_pheno.tsv')

res<-filter(res,guuid %in% all_guuids$guuid)
all_guuids<-filter(all_guuids,guuid %in% res$guuid)

guuids<-left_join(all_guuids,res,by=c("guuid"))

amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv")
amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv") %>% filter(`% Identity to reference sequence` ==100 & `% Coverage of reference sequence` ==100)


resfinder<-read_tsv('./data/resfinder_phenotypes.txt')
resfinder<-filter(resfinder,grepl('Piperacillin',Phenotype) & grepl('Tazo',Phenotype))
resfinder<-distinct(resfinder,`Gene_accession no.`)
resfinder$`Gene_accession no.`<-str_replace_all(resfinder$`Gene_accession no.`,'_.*','')
ag_res<-filter(guuids,Piptaz==1)
ag_counts<-filter(amrfinder,grepl('BETA-LACTAM',Class)) %>% filter(guuid %in% ag_res$guuid) %>% filter(`Gene symbol` %in% resfinder$`Gene_accession no.`) %>% group_by(`Gene symbol`) %>% count() %>% arrange(desc(n))
ag_genes<-unique(ag_counts$`Gene symbol`)
ag_genes<-ag_genes[!grepl('blaEC',ag_genes)]
ag_genes<-append(ag_genes,c('blaTEM-1','blaTEM'))

ag_amrfinder<-filter(amrfinder,guuid %in% guuids$guuid)

###### here we calculate stats
res_genes<-filter(ag_amrfinder,`Gene symbol` %in% ag_genes)
guuids$predict_res<-ifelse(guuids$guuid %in% res_genes$guuid,1,0)
cm<-caret::confusionMatrix(as.factor(guuids$predict_res),reference=as.factor(guuids$Piptaz),positive="1")

## prevelance
329/(329+7832)

tn<-cm$table[1,1]
fn<-cm$table[1,2]
tp<-cm$table[2,2]
fp<-cm$table[2,1]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
#sens
prop.test(x=(tp),n=(tp+fn))
#spec
prop.test(x=(tn),n=(tn+fp))
#npv
prop.test(x=(tn),n=(tn+fn))
#ppv
prop.test(x=(tp),n=(tp+fp))
#Maj
prop.test(x=(fp),n=(tn+fp))
#Vmj
prop.test(x=(fn),n=(tp+fn))



#######Coamox + blaTEM-1###########


all_guuids<-read_tsv('./data/guuids',col_names = c("guuid"))

res<-read_tsv('./data/global_pheno.tsv')

res<-filter(res,guuid %in% all_guuids$guuid)
all_guuids<-filter(all_guuids,guuid %in% res$guuid)

guuids<-left_join(all_guuids,res,by=c("guuid"))

amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv")
amrfinder<-read_tsv("./data/amrfinder_gladstone_merger.tsv") %>% filter(`% Identity to reference sequence` ==100 & `% Coverage of reference sequence` ==100)


resfinder<-read_tsv('./data/resfinder_phenotypes.txt')
resfinder<-filter(resfinder,grepl('Clav',Phenotype))
resfinder$`Gene_accession no.`<-str_replace_all(resfinder$`Gene_accession no.`,'_.*','')

ag_res<-filter(guuids,Coamox==1)
ag_counts<-filter(amrfinder,grepl('BETA-LACTAM',Class)) %>% filter(`Gene symbol` %in% resfinder$`Gene_accession no.`) %>%  filter(guuid %in% ag_res$guuid) %>% group_by(`Gene symbol`) %>% count() %>% arrange(desc(n))
ag_genes<-unique(ag_counts$`Gene symbol`)
ag_genes<-ag_genes[!grepl('blaEC',ag_genes)]
ag_genes<-append(ag_genes,c('blaTEM-1','blaTEM'))


ag_amrfinder<-filter(amrfinder,guuid %in% guuids$guuid)

###### here we calculate stats
res_genes<-filter(ag_amrfinder,`Gene symbol` %in% ag_genes)
guuids$predict_res<-ifelse(guuids$guuid %in% res_genes$guuid,1,0)
cm<-caret::confusionMatrix(as.factor(guuids$predict_res),reference=as.factor(guuids$Coamox),positive="1")

## prevelance
1621/(1621+3454)

tn<-cm$table[1,1]
fn<-cm$table[1,2]
tp<-cm$table[2,2]
fp<-cm$table[2,1]
#concordance
prop.test(x=(tp+tn),n=(tp+tn+fp+fn))
#sens
prop.test(x=(tp),n=(tp+fn))
#spec
prop.test(x=(tn),n=(tn+fp))
#npv
prop.test(x=(tn),n=(tn+fn))
#ppv
prop.test(x=(tp),n=(tp+fp))
#Maj
prop.test(x=(fp),n=(tn+fp))
#Vmj
prop.test(x=(fn),n=(tp+fn))



###########does TEM allele affect R/S#############
tem_fam<-filter(out_family,gene=='TEM-1')
table(tem_fam$allele)
tem_fam<-filter(out_family,allele %in% tem_fam$allele & class=='betalactam')
table(tem_fam$gene)


out_save$full_guuid<-paste0(out_save$guuid,'_',out_save$gene)
tem_fam$full_guuid<-paste0(tem_fam$guuid,'_',tem_fam$gene)

tem<-filter(out_save,full_guuid  %in% tem_fam$full_guuid & class =='betalactam')
table(tem$allele)
length(unique(tem$allele))
tem$allele<-ifelse(tem$allele ==1 | tem$allele ==3|tem$allele ==9 | tem$allele ==11,tem$allele,'other')
other<-filter(tem,allele=='other')
other_amrfinder<-read_tsv('./other_tem_resfinder')
other<-left_join(other,other_amrfinder,by=c("full_guuid"="SEQUENCE"))
other_exact<-filter(other,`%IDENTITY` ==100)
length(unique(other_exact$gene))
length(unique(other_exact$guuid))
table(other_exact$GENE)

pheno<-read_tsv('/data/global_pheno.tsv')
tem<-left_join(tem,pheno,by=c("guuid"="guuid"))

amrfinder<-read_tsv('./data/amrfinder_gladstone_merger.tsv')
tem<-filter(tem,guuid %in% amrfinder$guuid)
resfinder<-read_tsv('./data/resfinder_phenotypes.txt')
resfinder<-filter(resfinder,grepl('Clav',Phenotype))
resfinder$`Gene_accession no.`<-str_replace_all(resfinder$`Gene_accession no.`,'_.*','')
amrfinder<-filter(amrfinder,guuid %in% tem$guuid)
coamox_genes<-filter(amrfinder,grepl('BETA-LACTAM',Class)) %>% filter(`Gene symbol` %in% resfinder$`Gene_accession no.`) 
amrfinder<-filter(amrfinder,`Gene symbol` %in% coamox_genes$`Gene symbol`)
amrfinder<-select(amrfinder,guuid,gene=`Gene symbol`)
count<-amrfinder %>% group_by(gene) %>% count() %>% filter(n>=5)
amrfinder$gene<-ifelse(amrfinder$gene %in% count$gene,amrfinder$gene,'other')
amrfinder<-data.frame(table(amrfinder)) %>% spread(gene,Freq)

tem<-left_join(tem,amrfinder,by=c("guuid"))
tem<-filter(tem,!is.na(Coamox))
tem[is.na(tem)]<-0
write_tsv(tem,'tem_data1.tsv')
model<-logistf::logistf(Coamox ~ as.factor(allele) +`blaCMY-2` + `blaOXA-1` + other,data=tem)
model<-logistf::logistf(Coamox ~ as.factor(allele),data=tem)
model<-logistf::logistf(Coamox ~ `blaCMY-2`,data=tem)
model<-logistf::logistf(Coamox ~ `blaOXA-1`,data=tem)
model<-logistf::logistf(Coamox ~ other,data=tem)

summary(model)
exp(coef(model))
exp(confint(model))

amrfinder<-filter(amrfinder,grepl('bla',`Gene symbol`))
amrfinder<-filter(amrfinder,!grepl('blaEC',`Gene symbol`) & ! `Gene symbol` =='blaTEM-1')
tem$other<-ifelse(tem$guuid %in% amrfinder$guuid,'other','nother')
table(tem$allele,tem$Coamox,tem$other)

nother<-filter(tem,other=='nother')
t<-table(nother$Coamox,nother$allele)





###########

out_save$full_guuid<-paste0(out_save$guuid,'_',out_save$gene)
tem_fam$full_guuid<-paste0(tem_fam$guuid,'_',tem_fam$gene)

tem<-filter(out_save,full_guuid  %in% tem_fam$full_guuid)
table(tem$allele)
tem$allele<-ifelse(tem$allele ==1 | tem$allele ==3|tem$allele ==9 | tem$allele ==11,tem$allele,'other')
pheno<-read_tsv(''~/gn/gwas/predictio'./data/global_pheno.tsv')
tem<-left_join(tem,pheno,by=c("guuid"="guuid"))

amrfinder<-read_tsv('./data/amrfinder_gladstone_merger.tsv')
tem<-filter(tem,guuid %in% amrfinder$guuid)
resfinder<-read_tsv('./resfinder_phenotypes.txt')
resfinder<-filter(resfinder,grepl('Tazoba',Phenotype))
resfinder$`Gene_accession no.`<-str_replace_all(resfinder$`Gene_accession no.`,'_.*','')
amrfinder<-filter(amrfinder,guuid %in% tem$guuid)
coamox_genes<-filter(amrfinder,grepl('BETA-LACTAM',Class)) %>% filter(`Gene symbol` %in% resfinder$`Gene_accession no.`) 
amrfinder<-filter(amrfinder,`Gene symbol` %in% coamox_genes$`Gene symbol`)
amrfinder<-select(amrfinder,guuid,gene=`Gene symbol`)
count<-amrfinder %>% group_by(gene) %>% count() %>% filter(n>=5)
amrfinder$gene<-ifelse(amrfinder$gene %in% count$gene,amrfinder$gene,'other')
amrfinder<-data.frame(table(amrfinder)) %>% spread(gene,Freq)

tem<-left_join(tem,amrfinder,by=c("guuid"))
tem<-filter(tem,!is.na(Piptaz))
tem[is.na(tem)]<-0

model<-logistf::logistf(Piptaz ~ as.factor(allele) +`blaCMY-2` + `blaOXA-1` + other,data=tem)
model<-logistf::logistf(Piptaz ~ as.factor(allele),data=tem)
model<-logistf::logistf(Piptaz ~ `blaCMY-2`,data=tem)
model<-logistf::logistf(Piptaz ~ `blaOXA-1`,data=tem)
model<-logistf::logistf(Piptaz ~ other,data=tem)

summary(model)
exp(coef(model))
exp(confint(model))


###################

out<-out_save
out$location<-ifelse(grepl('-',out$guuid),'Oxford','Not Oxford')
t<-out %>% group_by(new_allele,location) %>% count()


t<-pivot_wider(t,names_from = location,values_from = n,id_cols = new_allele,values_fill = 0)
unique_oxford<-filter(t,`Not Oxford`==0)

six_times<-filter(t,  Oxford >=6)
fa<-ggplot(filter(t,Oxford >=1)) +
  aes(x=(Oxford),y=(`Not Oxford` )) +
  geom_point() +
  theme_minimal() + xlab("N times observed in Oxford") +
  ylab("N times observed in non-Oxford studies") + ggpubr::stat_cor(method = "spearman") #+ 
  scale_x_log10() + scale_y_log10() 
  

tt<-filter(t,Oxford>=1)
correlation<-cor.test(tt$Oxford,tt$`Not Oxford`,method = "spearman")
DescTools::SpearmanRho(tt$Oxford,tt$`Not Oxford`,conf.level=0.95)

out3=NULL
for(i in 1:10){
  m<-filter(t, Oxford >=i)
  unique<-filter(m,`Not Oxford`==0)
  out3<-rbind(out3,data.frame(i,(nrow(unique)/nrow(m))))
}

all<-out3


out<-filter(out_save,class== 'betalactam')
out$location<-ifelse(grepl('-',out$guuid),'Oxford','Not Oxford')


t<-out %>% group_by(new_allele,location) %>% count()

t<-pivot_wider(t,names_from = location,values_from = n,id_cols = new_allele,values_fill = 0)

plot(t$`Not Oxford`,t$Oxford)

out3=NULL
for(i in 1:10){
  m<-filter(t, Oxford >=i)
  unique<-filter(m,`Not Oxford`==0)
  out3<-rbind(out3,data.frame(i,(nrow(unique)/nrow(m))))
}

bl<-out3

out<-filter(out_save,class=='aminoglycoside')
out$location<-ifelse(grepl('-',out$guuid),'Oxford','Not Oxford')


t<-out %>% group_by(new_allele,location) %>% count()

t<-pivot_wider(t,names_from = location,values_from = n,id_cols = new_allele,values_fill = 0)

plot(t$`Not Oxford`,t$Oxford)

out3=NULL
for(i in 1:10){
  m<-filter(t, Oxford >=i)
  unique<-filter(m,`Not Oxford`==0)
  out3<-rbind(out3,data.frame(i,(nrow(unique)/nrow(m))))
}
aminog<-out3


out<-filter(out_save,class=='quinolone')
out$location<-ifelse(grepl('-',out$guuid),'Oxford','Not Oxford')


t<-out %>% group_by(new_allele,location) %>% count()

t<-pivot_wider(t,names_from = location,values_from = n,id_cols = new_allele,values_fill = 0)

plot(t$`Not Oxford`,t$Oxford)

out3=NULL
for(i in 1:10){
  m<-filter(t, Oxford >=i)
  unique<-filter(m,`Not Oxford`==0)
  out3<-rbind(out3,data.frame(i,(nrow(unique)/nrow(m))))
}
quinolone<-out3


out<-filter(out_save,class=='fosfomycin')
out$location<-ifelse(grepl('-',out$guuid),'Oxford','Not Oxford')


t<-out %>% group_by(new_allele,location) %>% count()

t<-pivot_wider(t,names_from = location,values_from = n,id_cols = new_allele,values_fill = 0)

plot(t$`Not Oxford`,t$Oxford)

out3=NULL
for(i in 1:10){
  m<-filter(t, Oxford >=i)
  unique<-filter(m,`Not Oxford`==0)
  out3<-rbind(out3,data.frame(i,(nrow(unique)/nrow(m))))
}
fosfo<-out3


out<-filter(out_save, class=="trimethoprim")
out$location<-ifelse(grepl('-',out$guuid),'Oxford','Not Oxford')


t<-out %>% group_by(new_allele,location) %>% count()

t<-pivot_wider(t,names_from = location,values_from = n,id_cols = new_allele,values_fill = 0)

plot(t$`Not Oxford`,t$Oxford)

out3=NULL
for(i in 1:10){
  m<-filter(t, Oxford >=i)
  unique<-filter(m,`Not Oxford`==0)
  out3<-rbind(out3,data.frame(i,(nrow(unique)/nrow(m))))
}
trim<-out3

out<-filter(out_save,class =='cephalosporin')
out$location<-ifelse(grepl('-',out$guuid),'Oxford','Not Oxford')


t<-out %>% group_by(new_allele,location) %>% count()

t<-pivot_wider(t,names_from = location,values_from = n,id_cols = new_allele,values_fill = 0)

plot(t$`Not Oxford`,t$Oxford)

out3=NULL
for(i in 1:10){
  m<-filter(t, Oxford >=i)
  unique<-filter(m,`Not Oxford`==0)
  out3<-rbind(out3,data.frame(i,(nrow(unique)/nrow(m))))
}
ceph<-out3

ceph$which<-"Cephalosporin"
trim$which<-"Trimethoprim"
fosfo$which<-"Fosfomycin"
quinolone$which<-"Quinolone"
all$which<-"All"
bl$which<-"Beta-lactam"
aminog$which<-"Aminoglycoside"

all<-rbind(all,bl,aminog,ceph,trim,fosfo,quinolone)

colors<-c( "#000000","#AD3B3B" ,"#205719","#2D48E0", "#32C2B1" ,"#F79C34", "#6D5EAB")
names(colors)<-c('All','Aminoglycoside','Beta-lactam','Quinolone','Cephalosporin','Fosfomycin','Trimethoprim')

names(all)<-c('i','prop','which')
fb<-ggplot(all) +
  aes(x=i,y=prop,color=which) +
  geom_point() + geom_path() + theme_minimal() +
  ylab("Proportion of ARG-alleles unique to Oxford") +
  xlab("ARG-allele observed at least N times in Oxford") +
  labs(color="Drug class") +
  scale_x_continuous(breaks = seq(0,10,1)) + 
  scale_color_manual(values = colors)

fa+fb


##########msa###########

library(ggmsa)
colors<-data.frame(names=c("A","C","T","G"), color=c("#FF0011" ,"#08FF29", "#0008FC", "#FFFF00"))
bseq<-ggmsa('./data/tem.fasta',start=133,end=143,custom_color = colors,seq_name = T)  + scale_x_continuous(breaks=seq(133,144,by=5))
aseq<-ggmsa('./data/tem.fasta',start=13,end=23,custom_color = colors,seq_name = T) + scale_x_continuous(breaks=seq(13,23,by=5))
cseq<-ggmsa('./data/tem.fasta',start=223,end=233,custom_color = colors,seq_name = T) + scale_x_continuous(breaks=seq(223,232,by=5))
dseq<-ggmsa('./data/tem.fasta',start=391,end=401,custom_color = colors,seq_name = T) + scale_x_continuous(breaks=seq(391,401,by=5))
eseq<-ggmsa('./data/tem.fasta',start=469,end=479,custom_color = colors,seq_name = T) + scale_x_continuous(breaks=seq(469,479,by=5))
fseq<-ggmsa('./data/tem.fasta',start=712,end=722,custom_color = colors,seq_name = T) + scale_x_continuous(breaks=seq(712,722,by=5))


aseq + bseq + cseq + dseq + eseq + fseq

#########check depth of classifications###########

singletons<-filter(count_table,Singleton==1)
double<-filter(count_table,Double==1)
tenner<-filter(count_table,Tenner==1)

out_save$which<-ifelse(out_save$new_allele %in% singletons$new_allele,"singletons",
                       ifelse(out_save$new_allele %in% double$new_allele, "doubles",
                              ifelse(out_save$new_allele %in% tenner$new_allele, "tenner","unknown" )))

all_depths<-read_tsv('./data/all_depth.tsv')

out_save<-left_join(out_save,all_depths,by=c("gene","guuid"))

quantile(out_save$depth,probs=c(0.025,0.975),na.rm = T)
out_save$singleton<-ifelse(out_save$which=="singletons","Singleton","Not Singleton")
#out_save$depth<-ifelse(out_save$depth >436,NA,out_save$depth)
kruskal.test(out_save$depth,out_save$singleton)

out_save %>% group_by(singleton) %>% summarise(d=median(depth,na.rm = T),q1=quantile(depth,probs=0.25,na.rm=T),q3=quantile(depth,probs=0.75,na.rm=T))

Fig_S1<-ggplot(out_save) +
  aes(x=depth) +
  geom_histogram() +
  facet_wrap(~singleton) + theme_minimal()  +
  xlim(c(0,437)) + 
  xlab("Sequencing depth") +
  ylab("Count")

########skesa shovill###############
ss_tem<-select(tem,full_guuid,allele)
singleton<-filter(count_table,Singleton==1)

out_sing<-out_save[out_save$new_allele %in% singleton$new_allele,]
out_sing$guuid_full<-paste0(out_sing$guuid,'_',out_sing$gene)

skesa_shovill<-read_tsv('./data/skesa_shovill.dist',col_names = F)
skesa_shovill<-filter(skesa_shovill,X1==X2)
skesa_shovill$single<-ifelse(skesa_shovill$X1 %in% out_sing$guuid_full & skesa_shovill$X2 %in% out_sing$guuid_full,1,0)
table(skesa_shovill$single)
skesa_shovill$exact<-ifelse(skesa_shovill$X3==0,"exact","not exact")
t<-table(skesa_shovill$single,skesa_shovill$exact)
t
92+12
12/104

116+6629
116/6745
fisher.test(t)

skesa_shovill<-left_join(skesa_shovill,ss_tem,by=c("X1"="full_guuid"))
skesa_shovill<-left_join(skesa_shovill,ss_tem,by=c("X2"="full_guuid"))

############# phenotypic variability#############

#The greatest phenotypic heterogeneity occured with cariage of blaTEM-1...
amrfinder<-read_tsv('./data/amrfinder_gladstone_merger.tsv')
tem<-filter(amrfinder,`Gene symbol` =='blaTEM-1')
tem<-filter(guuids,guuid %in% tem$guuid)
table(tem$Coamox)
1107+1023
1107/2130

table(tem$Piptaz)

table(tem$Piptaz)
197+2902
197/3099
