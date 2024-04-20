import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# t, ang_disp, ang_disp_sq = np.loadtxt("data/single_window_theta.dat", unpack=True)

# plt.figure(tight_layout=True)
# plt.plot(t, ang_disp)
# plt.xlabel(r"$t/\tau$", fontsize=18)
# plt.ylabel(r"$\Delta \theta$", fontsize=18)
# plt.savefig("plots/single_window_theta.pdf")

# plt.figure(tight_layout=True)
# plt.plot(t, ang_disp_sq)
# plt.xlabel(r"$t/\tau$", fontsize=18)
# plt.ylabel(r"$\Delta \theta^2$", fontsize=18)
# plt.savefig("plots/single_window_theta_sq.pdf")

t, ds_sq = np.loadtxt("single_window_displacement_sampled.dat", unpack=True)

Dimension= 1

def f(x, m, c):
    return m * x + c

popt, pcov = curve_fit(f, t, ds_sq)

D_measured = popt[0] / (2.0*Dimension)
#p	D_0	D_rot_0	D_perp	D_par	D_rot
p, D_0, D_rot_0, D_perp, D_par, D_rot = np.loadtxt("D_info.txt", unpack=True)

fitline = f(t, *popt)
D_m = D_measured * D_0
print("expected D =", D_m)

plt.figure(tight_layout=True)
plt.plot(t, ds_sq, 'k-' , label='data')
plt.plot(t, fitline, 'r--', label="Fit")
plt.xlabel(r"$t/\tau$", fontsize=18)
plt.ylabel(r"$\langle \Delta s^2 \rangle$", fontsize=18)
plt.legend()
plt.savefig("single_window_displacement_sampled.pdf")


t_rot, dtheta_sq = np.loadtxt('samp.txt', unpack = True)
fit_params, fit_cov = curve_fit(f, t_rot, dtheta_sq)

print('Fitted parameters: m = {:.4e}, c = {:.4e}'.format(fit_params[0], fit_params[1]))
print('Diffusion coefficient: D = {:.4e}'.format(fit_params[0]/(2)))

fitline = f(t_rot, *fit_params)

plt.figure(tight_layout=True)
plt.plot(t_rot, dtheta_sq, 'k-', label='Data')
plt.plot(t_rot, fitline, 'r--', label='Fit')
plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel(r'$\langle \theta^2 \rangle$', fontsize=18)
plt.legend(fontsize=14)
plt.savefig('D_rot_sampled.pdf')

# with open('D_rot_value.txt', 'w') as file:
#     file.write('# D\n')
#     file.write('{}\n'.format(fit_params[0]/(2)))
