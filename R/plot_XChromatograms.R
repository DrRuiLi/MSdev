plot_XChromatograms <- function(xchrom , norm = T,move = T){

extract_xchrom <- function(i,xchrom){
  x <- xchrom[i]

  rt <- rtime(x)
  intensity <- intensity(x)

  data.frame(adduct.id = i,rt , intensity)
}

if (norm) {
  xchrom <- normalise(xchrom)
  if (move) {
    chrom.data <- lapply(1:nrow(xchrom), extract_xchrom,xchrom)%>%
      data.table::rbindlist()%>%
      mutate(intensity = intensity*100,
             intensity = case_when(is.na(intensity)~ 0 ,
                                   T~intensity),
             rt = rt +adduct.id*3,
             intensity = intensity+adduct.id*3)%>%
      mutate(adduct.id = as.factor(adduct.id))
  }else{
    chrom.data <- lapply(1:nrow(xchrom), extract_xchrom,xchrom)%>%
      data.table::rbindlist()%>%
      mutate(intensity = intensity*100,
             intensity = case_when(is.na(intensity)~ 0 ,
                                   T~intensity))%>%
      mutate(adduct.id = as.factor(adduct.id))

  }

}else{
  chrom.data <- lapply(1:nrow(xchrom), extract_xchrom,xchrom)%>%
    data.table::rbindlist()%>%
    mutate(adduct.id = as.factor(adduct.id))
}



ggplot(chrom.data)+
  geom_line(aes(x = rt , y = intensity , col = adduct.id),size = 1)+
  theme_bw()



}
