import scipy
import scipy.integrate as integrate
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi
from find_W import W, rho_S

W_A1 = W([(1.45, -206e-9)],rho_S(2.196e6,206e-9))
W_A3 = W([(1.45, -197e-9), (1, -399e-9), (1.45, -197e-9)],rho_S(2.196e6,(197+197)*1e-9))
W_A5 = W([(1.45, -180e-9), (1, -421e-9), (1.45, -182e-9), (1, -421e-9), (1.45, -180e-9)], rho_S(2.196e6,(182+180*2)*1e-9))
W_A7 = W([(1.45, -156e-9), (1, -452e-9), (1.45, -161e-9), (1, -448e-9), (1.45, -161e-9), (1, -452e-9), (1.45, -156e-9)], rho_S(2.196e6,(156+161)*2e-9))

W_silica = (W_A1[0],W_A3[0],W_A5[0],W_A7[0])

print(W_silica)

from matplotlib.pyplot import figure
figure(num=None, figsize=(6, 6), dpi=80, facecolor='w', edgecolor='k')
plt.title('W for structures A1-A7 (Fig 2c)')
plt.xlabel('Number of layers')
plt.ylabel('W (\u221Ag/m)')
plt.ylim(0,0.18)
# plt.xlim(10**(-9),10**(-4))
plt.plot((1,3,5,7),W_silica,'bo--')
plt.xticks((1,3,5,7))
plt.yticks((0,0.05,0.10,0.15))
plt.grid(True)
plt.show()
