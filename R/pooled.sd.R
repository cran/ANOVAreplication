pooled.sd <- function(data){
  names(data) <- c("y","g")
  p <- length(table(data$g))
  n.g <- table(data$g)
  sd.g <- aggregate(data$y,by=list(data$g),sd)[,2]
  sqrt(sum((n.g - 1) * sd.g^2)/(sum(n.g) - p))
}
