import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

t,dx,dy,dz = np.loadtxt("CoM_displacement.dat", unpack=True, usecols=(0,1,2,3))

t_max = len(t)
sampling_window = int(0.01 * t_max)

ang_sq_avg_sampled = np.zeros(sampling_window)
t0_iter_list = np.arange(0, t_max- sampling_window)
print("Samples taken: {}".format(len(t0_iter_list)))

t_shortened = t[:sampling_window]
disp_sq_avg_sampled= np.zeros(len(t_shortened))
for t0_i in t0_iter_list:
    init_x = dx[t0_i]
    init_y = dy[t0_i]
    init_z = dz[t0_i]
    disp_sq_avg_sampled += (dx[t0_i:t0_i+sampling_window] - init_x)**2 + (dy[t0_i:t0_i+sampling_window] - init_y)**2 + (dz[t0_i:t0_i+sampling_window] - init_z)**2

disp_sq_avg_sampled /= len(t0_iter_list)

# plt.figure(tight_layout=True)
# plt.plot(t_shortened, disp_sq_avg_sampled)
# plt.xlabel('Avg Displacement squared')
# plt.ylabel('t')
# plt.savefig('single_window_displacement_sampled.pdf')

with open("window_displacement.dat", "w") as file:
    for t_i in range(sampling_window):
        file.write("{}\t{}\n".format(t_shortened[t_i], disp_sq_avg_sampled[t_i]))

t, ds_sq = np.loadtxt("window_displacement.dat", unpack=True)

Dimension= 3

def f(x, m, c):
    return m * x + c

popt, pcov = curve_fit(f, t, ds_sq)

D_measured = popt[0] / (2.0*Dimension)
#p	D_0	D_rot_0	D_perp	D_par	D_rot
p, D_0, D_rot_0, D_perp, D_par, D_rot = np.loadtxt("D_info.txt", unpack=True)

print("D Measured = " ,D_measured)

fitline = f(t, *popt)
D_m = D_measured * D_0
print("expected D =", D_m)

plt.figure(tight_layout=True)
plt.plot(t, ds_sq, 'k-' , label='data')
plt.plot(t, fitline, 'r--', label="Fit")
plt.xlabel(r"$t/\tau$", fontsize=18)
plt.ylabel(r"$\langle \Delta s^2 \rangle$", fontsize=18)
plt.legend()
plt.savefig("D_disp_sampled_Davg.pdf")
# plt.show() 
