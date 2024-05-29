import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

mpl.use("pdf")

mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.family"] = "serif"

cm = 1 / 2.54
mpl.rcParams["figure.figsize"] = [12.00 * cm, 8.00 * cm]
mpl.rcParams["figure.constrained_layout.use"] = True

def func_exp(x, A, B, C):
    return A + B * np.exp(-(35.0 - x)/C)

fig = plt.figure()
data = np.loadtxt("002.txt")
index = -20 + np.argmax(data[-20:, 1]) - 25
plt.step(data[:,0], data[:,1]/data[index+25,1], label = "$P<0$", color = "#5a7ec8")
popt, pcov = curve_fit(func_exp, data[30:index:,0], data[30:index,1], p0=[0.5, 1.0, 1.5])
plt.plot(data[30:index,0], func_exp(data[30:index,0], popt[0], popt[1], popt[2])/data[index+25,1], linestyle = "dashed", color = "#002269")
data = np.loadtxt("001.txt")
index = -20 + np.argmax(data[-20:, 1])
plt.step(data[:,0], data[:,1]/data[index,1], label = "$P=0$", color = "#b96bcd")
popt, pcov = curve_fit(func_exp, data[60:index:,0], data[60:index,1], p0=[0.5, 1.0, 1.5])
plt.plot(data[60:index,0], func_exp(data[60:index,0], popt[0], popt[1], popt[2])/data[index,1], linestyle = "dashed", color = "#6a0084")
data = np.loadtxt("003.txt")
index = -20 + np.argmax(data[-20:, 1])
plt.step(data[:,0], data[:,1]/data[index,1], label = "$P>0$", color = "#d37b81")
popt, pcov = curve_fit(func_exp, data[85:index:,0], data[85:index,1], p0=[0.25, 10.0, 0.5])
plt.plot(data[85:index,0], func_exp(data[85:index,0], popt[0], popt[1], popt[2])/data[index,1], linestyle = "dashed", color = "#9f000b")

plt.legend(loc='upper left')
plt.xlabel("$R$ (a.u.)")
plt.ylabel("$\\rho_n(R)$ (a.u.)")
plt.xlim(17.485,34.97)
plt.ylim(0.0,1.08)

plt.savefig("plot-exp-fit.pdf")
