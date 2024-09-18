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




###colors
cols <- c("Beijerinckiaceae"="#D38E18",  "Acetobacteraceae"= "#CB1B06", "Sphingomonadaceae" = "#A37FBF", "Acidobacteriaceae" = "#2D4FB1", "Nostocaceae" = "#1DA6B5","UBA10450"="#8CB28D")




