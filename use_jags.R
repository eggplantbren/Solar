model = "model
{
    period <- 365.25
    pi <- 3.141592653589793

    sd ~ dunif(0, 1)

    C ~ dunif(-10, 10)
    A ~ dunif(0, 10)

    phi ~ dunif(0, 2.0*pi)
    k ~ dunif(0, 1)

    for(i in 1:length(y))
    {
        top[i] <- exp(C + A*sin(2*pi*t[i]/period + phi))
        y[i] ~ dnorm(k*top[i], 1/(sd*top[i])^2)T(0, top[i])
    }
}
"



data = read.table("data.txt", header=FALSE)
colnames(data) = c("t", "y")

# Variables to monitor
variable_names = c("A", "C", "sd", "phi", "k")

# How many proper steps?
steps = 100000

# Thinning?
thin = 10

# Import the rjags library
library("rjags")

# Create a JAGS model object
jm = jags.model(textConnection(model), data,
                inits=list(list(C=9, A=1, phi=1, sd=0.5)))

# Do some MCMC
results = coda.samples(jm, variable_names, steps, thin=thin)

# Extract chains as data frame
results = as.data.frame(results[[1]])
