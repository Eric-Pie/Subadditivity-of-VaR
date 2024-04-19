#install.packages("ggplot2")
library(ggplot2)

#定义显性函数
explicit_function <- function(F_s, d) {
  sum <- 0
  for (n in 0:(d-1)) {
    sum <- sum + (1/factorial(n)) * (1-alpha)^d * (-d * log(1-alpha))^n
  }
  return(1 - sum)
}

#创建数据框并计算F_s值
alpha <- 1-1/exp(1) # 设置alpha的值
df <- data.frame(d = 1:2000, F_s = rep(NA, 2000))

for (d in 1:2000) {
  F_s <- explicit_function(F_s, d)
  df$F_s[d] <- F_s
}

#使用以下代码使用ggplot2库绘制二维图
plot <- ggplot(df, aes(x = d, y = F_s)) +
  geom_line() +
  geom_point() +
  labs(x = "d", y = "F_s") +
  ggtitle("d vs. F_s") +
  scale_x_continuous(limits = c(0, 2000))

print(plot)