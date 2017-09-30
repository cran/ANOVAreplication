pooled.sd <- function(data){
  p <- length(table(data$g))
  n.g <- table(data$g)
  sd.g <- aggregate(data$y,by=list(data$g),sd)[,2]
  sum((n.g-1)*sd.g)/(sum(n.g)-p)
}
