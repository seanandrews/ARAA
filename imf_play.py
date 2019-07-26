import numpy as np
import matplotlib.pyplot as plt

# masses
m = np.logspace(-2, 2, 1024)

print(m)

# chabrier IMF
imf_0 = 0.093 * np.exp(-0.5*(np.log10(m) - np.log10(0.2))**2/(0.55**2))
imf_0[m >= 1.] = 0.041*m[m >= 1.]**(-1.35)


plt.loglog(m, imf_0)
plt.show()
