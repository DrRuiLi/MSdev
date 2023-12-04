ggplot_sdf <- function(sdf,
                       cex = 1,
                       show_ele = F){


  sdf.formula <- MF(sdf,addH=T)
  sdf.mz <- chemform_mz(sdf.formula)%>%round(digits = 4)
  atom.data <- atomblock(sdf)[,1:2]%>%
    `colnames<-`(c("x","y"))%>%
    as.data.frame()%>%
    rownames_to_column("Atom_id" )%>%
    dplyr::mutate(element = str_extract(Atom_id,
                                        "[:alpha:]*"))
  bond.length.short <- ifelse(show_ele,0.1,0)
  bond.data <- bondblock(sdf)[,1:3]%>%
    `colnames<-`(c("from","to","bond_type"))%>%
    as.data.frame()%>%
    dplyr::mutate(
      bond_id = 1:n(),
      x = atom.data$x[from],
      xend = atom.data$x[to],
      y = atom.data$y[from],
      yend = atom.data$y[to]
    )%>%
    dplyr::mutate(xl = (xend-x),
                  yl = (yend - y),
                  x = x + bond.length.short*xl,
                  xend = xend - bond.length.short*xl,
                  y = y+bond.length.short*yl,
                  yend = yend - bond.length.short*yl)
  for (i in 1:nrow(bond.data)) {
    bond.data <- dplyr_copy_row(bond.data,
                                i,
                                bond.data$bond_type[i]-1)
  }

  lw <- 0.5*cex
  sw <- 0.5*cex
  col.bond <- "#666666"
  ggplot()+
    ### 3 bond
    geom_segment(data = filter(bond.data,bond_type == 3),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = col.bond,linewidth = lw+2*sw+2*lw)+
    geom_segment(data = filter(bond.data,bond_type == 3),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = "white",linewidth = lw+2*sw)+
    geom_segment(data = filter(bond.data,bond_type == 3),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = col.bond,linewidth = lw)+
    ### 2 bond
    geom_segment(data = filter(bond.data,bond_type == 2),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = col.bond,linewidth = 2*lw+sw)+
    geom_segment(data = filter(bond.data,bond_type == 2),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = "white",linewidth = sw)+
    ### 1 bond
    geom_segment(data = filter(bond.data,bond_type == 1),
                 aes(x = x,xend = xend ,
                     y = y,yend = yend),
                 col = col.bond,linewidth = lw)+
    geom_text(aes(x = median(range(atom.data$x)),
                  y = max(atom.data$y)+diff(range(atom.data$y))*0.5,
                  label = paste(sdf.formula,"\n",sdf.mz)),
              size = 2)+
    ylim(c(min(atom.data$y),max(atom.data$y)+diff(range(atom.data$y))*0.8))+
    xlim(expand_range(range(atom.data$x),multi = 0.2))+
    theme_void()->p


  if (show_ele) {
    p <- p+geom_text(data = atom.data,
                aes(x = x, y = y ,label = element),
                size = 2 *cex)
  }else{
    p <- p+geom_point(data = atom.data,
                     aes(x = x, y = y ),
                     size = 0.5*cex)
  }
  p
  return(p)
}
