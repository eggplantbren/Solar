import dnest4.classic as dn4
dn4.postprocess()

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data.txt")
posterior_sample = np.loadtxt("posterior_sample.txt")

t = np.linspace(0.0, 365.0, 10001)
plt.plot(data[:,0], data[:,1])

for i in range(min([100, posterior_sample.shape[0]])):
    C, A, T, phi, L = posterior_sample[i, :]
    top = C + A*np.sin(2.0*np.pi*t/T + phi)
    plt.plot(t, np.exp(top), color="g", alpha=0.1)
    plt.plot(t, np.exp(top - L), color="r", alpha=0.1)
plt.show()


top = np.empty(posterior_sample.shape[0])
for i in range(posterior_sample.shape[0]):
    C, A, T, phi, L = posterior_sample[i, :]
    top[i] = np.exp(C + A)
plt.hist(top, 100)
print(top.mean(), top.std())
plt.show()

