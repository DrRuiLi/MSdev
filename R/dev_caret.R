#' @title ggplot_roc
#'
#' @description
#' plot ROC of train and test
#'
#'
#' @param roc.train roc
#' @param roc.test roc
#'
#' @return ggplot
#' @export
#'

ggplot_roc <- function(roc.train,roc.test){


  #roc.train <- roc(sample(letters[1:2],replace = T,1000),runif(1000))
  #roc.test <- roc(sample(letters[1:2],replace = T,1000),runif(1000))
  col.train <- "#DC0000"
  col.test <- "#3C5488"
  pROC::ggroc(list(train = roc.train,test = roc.test))+
    geom_abline(slope = 1,intercept = 1,lty = "dashed",col = "grey")+
    scale_color_manual(values = c(train = col.train , val = col.test))+
    annotate("text" ,x = 0.3,y=0.3, label = paste0("Train Set\nAUC = " , round(roc.train$auc,3),"(",
                                                   paste0(round(ci.auc(roc.train)[c(1,3)],3),collapse = "~"),")") ,
             size = 2.8,col = col.train)+
    annotate("text" ,x = 0.3,y=0.1, label = paste0("Test Set\nAUC = " , round(roc.test$auc,3),"(",
                                                   paste0(round(ci.auc(roc.test)[c(1,3)],3),collapse = "~"),")") ,
             size = 2.8,col = col.test)+
    coord_fixed(ratio = 1)+
    theme_bw()+
    theme(legend.key.size = unit(0.1,"inch"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.position = "none",
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 8),
          axis.ticks = element_line(size = 0.3),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.grid.major = element_blank(),
          panel.border = element_rect(fill= NA,size = 0.1),
          text = element_text(size=8))->p
  return(p)
}
