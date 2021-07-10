library(volesti)

#-- number of assets --#
d=30

#-- number of assets' returns --#
N = 500000

#-- generate random covariance matrix and mean for the assets' returns --#
Sigma = replicate(d+10, rnorm(d))
Sigma = Sigma %*% t(Sigma)
mu = rnorm(d)

#-- sample N assets' returns --#
X = MASS::mvrnorm(N, mu, Sigma)

##-- sample a random portfolio --#
ptf = sample_points(N = 1, exact = TRUE, body = "canonical simplex", Parameters = list("dimension"=d))
scores = numeric(length = N)

##-- compute the scores of the portfolio for those N assets' returns --#
for (i in 1:N) {
  scores[i] = SliceOfSimplex(X[i,], X[i,]%*%ptf)
}

##-- plot the distribution of scores --#
hist(scores, breaks = 100,  col = 'skyblue3')
