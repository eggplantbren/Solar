data = read.csv("solar.csv")

t = as.integer(as.Date(data$date))
t = t - min(t) + 1
y = data$energy_kwh
plot(t, y, type="l", ylim=c(0, 1.05*max(y)))

