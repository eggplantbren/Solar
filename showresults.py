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
    C, A, T, phi, L = posterior_sample[i, :]
    top = np.exp(C + A*np.sin(2.0*np.pi*t/T + phi))
    avg = top/np.exp(L)
    tot[i] = np.trapz(avg, x=t)
    peak[i] = np.max(top)
    if i < 100:
        plt.plot(t, top, color="g", alpha=0.1)
        plt.plot(t, avg, color="r", alpha=0.1)
plt.show()

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.hist(tot, 100)
plt.xlabel("Annual total (kWh)")
plt.subplot(1, 2, 2)
plt.hist(peak, 100)
plt.xlabel("Peak")
plt.show()

