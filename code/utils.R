#function that takes BAT-produced taxonomy assignements for bacteria and reduces them to the specified rank

gtdb_get_clade = function(s, clade="p"){
  if (is.na(s)){
    return(NA)
  }else{
    s<-as.character(s)
    s = str_split(s, pattern = ";", simplify = F)[[1]] 
    s = str_split(s, pattern = "__", simplify = T) %>% as.data.frame %>%
      filter(V1 == clade)
    name = s[1,2]
    if(name == ""){
      return("Unknown")
    }else{
      return(as.character(name))
    }
  }
}




#old stuff
library(tidyverse)



#define function for handling reports from blast,  bbduk, and quast

compile_data<-function(blast_path,bbduk_path,quast_path){
#read the tables
blast<-read.delim(blast_path,header=F)
bbduk<-read.delim(bbduk_path,header=F)
quast<-read.delim(quast_path,header=F)
colnames(blast)<-c('nreads','seed','yeast','blast_num')
colnames(bbduk)<-c('nreads','seed','yeast','read_num')
colnames(quast)<-c('nreads','seed','genome','cov')

#spread quast
quast2<-quast %>% spread(genome,cov,fill=0)


#join the tables
df<-left_join(blast,bbduk)
df<-left_join(df,quast2)

##transform nreads into numeric form
df<-transform(df, nreads=sub(".*nreads_([0-9]+).*", "\\1", nreads)) 
df$nreads<-as.numeric(df$nreads)
df <- df %>% mutate(bp=nreads*250)
df <-df %>% arrange(bp)

#make binary columns
df<-df %>% mutate(blast_abs=ifelse(blast_num>0,1,0)) %>%
  mutate(read_abs=ifelse(read_num>0,1,0)) 
return(df)
}


###colors
col_pres_abs<-"#5d0f72"
col_cypho<-"#f6a746"
col_cypho_light<-'#face98'
col_tremella<-"#464fcc"
col_tremella_light<-"#848add"
col_lecanoro<-"#ac2323"
col_alga<-"#09703a"
col_X12<-"#0c8d69"
col_GT57<-"#d03610"
col_G1<-"#811c02"




