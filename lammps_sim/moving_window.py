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

plt.figure(tight_layout=True)
plt.plot(t_shortened, disp_sq_avg_sampled)
plt.xlabel('Avg Displacement squared')
plt.ylabel('t')
plt.savefig('single_window_displacement_sampled.pdf')

print("Writing data...")
with open("single_window_displacement_sampled.dat", "w") as file:
    for t_i in range(sampling_window):
        file.write("{}\t{}\n".format(t_shortened[t_i], disp_sq_avg_sampled[t_i]))
# moving window for angular displacement 
t, dtheta= np.loadtxt('angular_displacement.dat', unpack=True)

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

plt.figure(tight_layout=True)
plt.plot(t_shortened, ang_sq_avg_sampled, 'k-', label='Data')
plt.ylabel(r'$\langle \Delta s^2 \rangle$')
plt.xlabel(r'$t$')
plt.savefig('single_window_Drot_data.pdf')

with open('single_window_Drot_data.txt', 'w') as file:
    file.write("# t0 \t <theta^2>\n")
    for i in range(sample_window):
        file.write("{} \t {}\n".format(t_shortened[i], ang_sq_avg_sampled[i]))

