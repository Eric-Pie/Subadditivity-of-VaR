alpha <- 1-10^(seq(log10(1-0.6), log10(1-0.999), length.out = 201))
values <- 1:201
plot(1-alpha, values, log = "x", xlab = expression(1-alpha))

## when 1-alpha close to 0,then more value will be pluged in, therefore we can observe more information when value is close to 0
## In this case, 1-alpha close to 0 represents that alpha is close to 1.