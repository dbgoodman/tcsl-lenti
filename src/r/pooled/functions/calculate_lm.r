calculate_lm <- function(data.dt, group.1, group.2){
  
  col.names <- c(group.1, group.2)
  
  m <- lm(data.dt[, col.names, with = F])
  
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                   list(a = as.character(format(coef(m)[1], digits = 3)),
                        b = as.character(format(coef(m)[2], digits = 3)),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  
  return(as.character(as.expression(eq)))
  
}