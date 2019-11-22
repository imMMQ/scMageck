dis_grna <- function(data){
  data <- read.delim(data)
  L <- unique(data$cell)
  output <- NULL
  for (i in L) {
    a <- nrow(subset(data, cell == i))
    m <- matrix(c(i, a), ncol = 2, dimnames = list(i, c("cell", "number_gRNA")))
    output <- rbind(output, m)
  }
  
  output <- as.data.frame(output)
  output$number_gRNA <- as.numeric(output$number_gRNA)
  
  N <- unique(output$number_gRNA)
  output_1 <- NULL
  for (t in N) {
    b <- nrow(subset(output, number_gRNA == t))
    c <- matrix(c(t, b), ncol = 2, dimnames = list(t, c("number_gRNA", "ncell")))
    output_1 <- rbind(output_1, c)
  }
  
  # cell0 <- matrix(c("0", est_cell - count(output)$n), ncol = 2, dimnames = list(0, c("number_gRNA", "ncell")))
  # output_1 <- rbind(cell0, output_1)
  
  output_1 <- as.data.frame(output_1)
  ggplot(output_1, aes(number_gRNA, ncell)) +
    geom_col(width = 0.4, fill = "#56B4E9", colour = "black") + 
    ggtitle("sgRNA distribution") + geom_text(aes(label = ncell), size = 3, hjust = 0.5, vjust = .01)
}

L <- unique(data$cell)
output <- NULL
for (i in L) {
  a <- nrow(subset(data, cell == i))
      m <- matrix(c(i, a), ncol = 2, dimnames = list(i, c("cell", "number_gRNA")))
      output <- rbind(output, m)
}

output <- as.data.frame(output)
output$number_gRNA <- as.numeric(output$number_gRNA)

N <- unique(output$number_gRNA)
output_1 <- NULL
for (t in N) {
  b <- nrow(subset(output, number_gRNA == t))
  c <- matrix(c(t, b), ncol = 2, dimnames = list(t, c("number_gRNA", "ncell")))
  output_1 <- rbind(output_1, c)
}

# cell0 <- matrix(c("0", est_cell - count(output)$n), ncol = 2, dimnames = list(0, c("number_gRNA", "ncell")))
# output_1 <- rbind(cell0, output_1)

output_1 <- as.data.frame(output_1)
ggplot(output_1, aes(number_gRNA, ncell)) +
  geom_col(width = 0.4, fill = "#56B4E9", colour = "black") + 
  ggtitle("sgRNA distribution") + geom_text(aes(label = ncell), size = 3, hjust = 0.5, vjust = .01)