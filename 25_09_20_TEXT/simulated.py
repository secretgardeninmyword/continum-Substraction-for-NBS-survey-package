import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit 










def func2(x, a, b, c):                    
    return a*np.exp(-(x-b)**2/(2*c**2))
x = np.linspace(0, 0.1, 100)
x1=[0.0549,0.0513,0.0572,0.0443]
y1=[41.8,34.8,29.8,2.7]
popt, _ = curve_fit(func2, x1, y1)
a1=popt[0]
b1=popt[1]
c1=popt[2]
plt.plot(x,func2(x,a1,b1,c1))
plt.scatter(x1,y1)
plt.show()
