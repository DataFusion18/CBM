rm(list=ls())

mcmc <- function(n, k = 2, N = 5000)
{
  x <- 1:n;
  res <- numeric(N)
  for(i in 1:N)
  {
    swap <- sample(1:n, k)
    x[swap] <- sample(x[swap],k);
    res[i] <- x[1];
  }
  return(res);
}

n <- 99; mu <- sum(1:n)/n;

mcmc(n) -> r1
plot(cumsum(r1)/1:length(r1), type="l", ylim=c(0,n), ylab="mean")
abline(mu,0,lty=2)

mcmc(n,round(n/2)) -> r2
lines(1:length(r2), cumsum(r2)/1:length(r2), col="blue")

mcmc(n,n) -> r3
lines(1:length(r3), cumsum(r3)/1:length(r3), col="red")

legend("topleft", c("k = 2", paste("k =",round(n/2)), paste("k =",n)), col=c("black","blue","red"), lwd=1)


K <- 5000;
M1 <- numeric(K)
M2 <- numeric(K)
M3 <- numeric(K)
for(i in 1:K)
{
  M1[i] <- mean(mcmc(n,2,100));
  M2[i] <- mean(mcmc(n,round(n/2),100));
  M3[i] <- mean(mcmc(n,n,100));
}

dev.new()
par(mfrow=c(3,1))
hist(M1, xlim=c(0,n), freq=FALSE)
hist(M2, xlim=c(0,n), freq=FALSE)
hist(M3, xlim=c(0,n), freq=FALSE)



