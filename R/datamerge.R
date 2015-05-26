# A function I wrote to merge together two data frames where one or both have variables the other doesn't
# Note that this function assumes variables of the same name in different dataframes are the same variable
Datamerge <- function(data1, data2) {
  add1 <- setdiff(colnames(data2), colnames(data1))
  add2 <- setdiff(colnames(data1), colnames(data2))
  hold1 <- data.frame(matrix(data = NA, nrow = nrow(data1), ncol=length(add1)))
  hold2 <- data.frame(matrix(data = NA, nrow = nrow(data2), ncol=length(add2)))
  colnames(hold1) <- add1
  colnames(hold2) <- add2
  end1 <- cbind(data1, hold1)
  end2 <- cbind(data2, hold2)
  final <- rbind(end1, end2)
  final
}