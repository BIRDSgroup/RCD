n <- 0:100

temp1 <- dpois(n, lambda = 5)
df <- data.frame(temp1)
df$temp2 <- dpois(n, 10)
df$temp3 <- dpois(n, 15)

library(ggplot2)

p <- ggplot(df) + geom_line(aes(x=n, y=temp1), color='red') + geom_point(aes(x=n, y=temp1), color='red')

p <- p + geom_line(aes(x=n, y=temp2), color='blue') + geom_point(aes(x=n, y=temp2), color='blue')

p <- p + geom_line(aes(x=n, y=temp3), color='green') + geom_point(aes(x=n, y=temp3), color='green')

p <- p + xlab('k') + ylab('Pr(x=k)') + xlim(0,30) + theme(panel.background = element_blank())
p
ggsave('poisson.svg', p)
t# p <- ggplot(df) + geom_density( aes(x=temp1), color = 'red') 
# p <- p + geom_density(aes(x=temp2), color='blue') 
# p <- p + geom_density(aes(x=temp3), color='green')
# 
# p
