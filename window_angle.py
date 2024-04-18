import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#windowing average for the angular displacement 
t, dtheta= np.loadtxt('angular_displacement.dat', unpack=True)
p, D_0, D_rot_0, D_perp, D_par, D_rot = np.loadtxt("D_info.txt", unpack=True)

t_max = len(t)
sample_window = int(0.01*t_max)
ang_sq_avg_sampled = np.zeros(sample_window)
t0_iter_list = np.arange(0, t_max-sample_window)
print("Samples taken: {}".format(len(t0_iter_list)))

t_shortened = t[:sample_window]
for t0_i in t0_iter_list:
    init_angle = dtheta[t0_i]
    ang_sq_avg_sampled += (dtheta[t0_i:t0_i+sample_window]-init_angle)**2

ang_sq_avg_sampled /= len(t0_iter_list)

#saving the averaged data 
with open('window_Drot_data.txt', 'w') as file:
    file.write("# t0 \t <theta^2>\n")
    for i in range(sample_window):
        file.write("{} \t {}\n".format(t_shortened[i], ang_sq_avg_sampled[i]))
#linear regression
def f(x, m, c):
    return m * x + c

t_rot, dtheta_sq = np.loadtxt('window_Drot_data.txt', unpack = True)
fit_params, fit_cov = curve_fit(f, t_rot, dtheta_sq) 

fitline = f(t_rot, *fit_params)

print('Fitted parameters: m = {:.4e}, c = {:.4e}'.format(fit_params[0], fit_params[1]))
print('Diffusion coefficient (measured): D = {:.4e}'.format(fit_params[0]/(2)))
print('Diffusion coefficient (actual): D = {:.4e}'.format((D_rot_0)*fit_params[0]/(2)))

#ploting 
plt.figure(tight_layout=True)
plt.plot(t_rot, dtheta_sq, 'k-', label='Data')
plt.plot(t_rot, fitline, 'r--', label='Fit')
plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel(r'$\langle \theta^2 \rangle$', fontsize=18)
plt.legend(fontsize=14)
plt.savefig('D_rot_sampled.pdf')

with open('D_rot_value.txt', 'w') as file:
    file.write('# D\n')
    file.write('{}\n'.format(fit_params[0]/(2)))
