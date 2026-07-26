import dnest4.classic as dn4
dn4.postprocess()

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data.txt")
posterior_sample = np.loadtxt("posterior_sample.txt")

t = np.linspace(0.0, 365.25, 10001)
plt.plot(data[:,0], data[:,1])

tot = np.empty(posterior_sample.shape[0])
peak = np.empty(posterior_sample.shape[0])
for i in range(posterior_sample.shape[0]):
    C, A, T, phi, mu, beta, L = posterior_sample[i, :]
    top = np.exp(C + A*np.sin(2.0*np.pi*t/T + phi))
    middle = top/(1.0 + np.exp(-mu))
    if i < 100:
        plt.plot(t,    top, color="g", alpha=0.1)
        plt.plot(t, middle, color="r", alpha=0.1)
plt.show()

