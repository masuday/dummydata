# Drawing a lactation curve
data <- read.table("lc1.txt")
colnames(data) <- c("x","y")

# add the transformed values to the data frame
data <- cbind(data, log(data$x), log(data$y))
colnames(data) <- c("x","y","logx","logy")

#
# log-transformed linear model
#
result <- lm(logy ~ 1 + logx + x, data)
a <- exp( unname(coefficients(result))[1] )
b <- unname(coefficients(result))[2]
c <- unname(coefficients(result))[3]

# plot
x <- 1:305
y <- a*(x^b)*exp(c*x)
plot(data$x, data$y, type="l")
points(data$x, data$y, col="red", pch="*")

#
# non-linear model
#
result2 <- nls(y ~ a*(x^b)*exp(c*x), data, start=list(a=20,b=0.2,c=-0.003))
A <- unname(coefficients(result2))[1]
B <- unname(coefficients(result2))[2]
C <- unname(coefficients(result2))[3]

# plot
Y <- A*(x^B)*exp(C*x)
plot(x, y, type="l", ylim=c(20,34))
par(new=TRUE)
plot(x, Y, type="l", ylim=c(20,34), col="blue")
points(data$x, data$y, col="red", pch="*", ylim=c(20,34))
