equation <- function(alpha, d) {
  sum <- 0
  for (n in 0:(d-1)) {
    sum <- sum + (1/factorial(n)) * (1-alpha)^d * (-d * log(1-alpha))^n
  }
  return(1 - alpha - sum)
}
bisection <- function(func, a, b, tol, max_iter) {
  iter <- 0
  while (abs(b - a) > tol && iter < max_iter) {
    c <- (a + b) / 2
    if (func(a) * func(c) < 0) {
      b <- c
    } else {
      a <- c
    }
    iter <- iter + 1
  }
  return((a + b) / 2)
}

tol <- 1e-6  # 定义误差容限
max_iter <- 100  # 定义最大迭代次数

alpha_list <- list()

for (d in 1:100) {  # 设置d的范围
  alpha <- bisection(function(x) equation(x, d), 0, 1, tol, max_iter)
  alpha_list[[as.character(d)]] <- alpha
}

#创建一个包含d和对应alpha值的数据框
df <- data.frame(d = 1:100, alpha = rep(NA, 100))

for (d in 1:100) {
  alpha <- bisection(function(x) equation(x, d), 0, 1, tol, max_iter)
  df$alpha[d] <- alpha
}

## plot the alpha with respect to d
ggplot(df, aes(x = d, y = alpha)) +
  geom_line() +
  labs(title = "alpha with respect to d",
       x = "d",
       y = "alpha") +
  theme_minimal()

