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
plt.step(data[:,0], data[:,1], label = "results for $P<0$", color = "#5a7ec8")
index = -20 + np.argmax(data[-20:, 1]) - 25
popt, pcov = curve_fit(func_exp, data[30:index:,0], data[30:index,1], p0=[0.5, 1.0, 1.5])
plt.plot(data[30:index,0], func_exp(data[30:index,0], popt[0], popt[1], popt[2]), label = "fit for $P<0$ ($\\xi={:.2f}$)".format(popt[2]), linestyle = "dashed", color = "#002269")
data = np.loadtxt("001.txt")
plt.step(data[:,0], data[:,1], label = "results for $P=0$", color = "#b96bcd")
index = -20 + np.argmax(data[-20:, 1])
popt, pcov = curve_fit(func_exp, data[60:index:,0], data[60:index,1], p0=[0.5, 1.0, 1.5])
plt.plot(data[60:index,0], func_exp(data[60:index,0], popt[0], popt[1], popt[2]), label = "fit for $P=0$ ($\\xi={:.2f}$)".format(popt[2]), linestyle = "dashed", color = "#6a0084")
data = np.loadtxt("003.txt")
plt.step(data[:,0], data[:,1], label = "results for $P>0$", color = "#d37b81")
index = -20 + np.argmax(data[-20:, 1])
popt, pcov = curve_fit(func_exp, data[85:index:,0], data[85:index,1], p0=[0.25, 10.0, 0.5])
plt.plot(data[85:index,0], func_exp(data[85:index,0], popt[0], popt[1], popt[2]), label = "fit for $P>0$ ($\\xi={:.2f}$)".format(popt[2]), linestyle = "dashed", color = "#9f000b")

plt.legend(loc='upper left')
plt.xlabel("$R$ (a.u.)")
plt.ylabel("$\\rho(R)$ (a.u.)")
plt.xlim(17.485,34.97)

plt.savefig("plot-exp-fit.pdf")
