import numpy as np 
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#calculating translational displacement 

t,dx,dy,dz = np.loadtxt("CoM_displacement.dat", unpack=True, usecols=(0,1,2,3))

inital_pos = np.array([dx[0], dy[0], dz[0]])
disp_vector = np.zeros((len(t),3))
with open("dis_data.txt", 'w') as f:
    f.write('#t dx dy dz\n')
    for i in range(len(t)):
        new_pos = np.array([dx[i], dy[i], dz[i]])
        d = np.subtract(new_pos,inital_pos)
        disp_vector[i,0]= d[0]
        disp_vector[i,1]= d[1]
        disp_vector[i,2]= d[2]
        f.write("{} {:.4f} {:.4f} {:.4f}\n".format(t[i], disp_vector[i,0], disp_vector[i,1], disp_vector[i,2])) 

#getting  the head vector 
t,e = np.loadtxt("txt", unpack=True)
f_vector = np.zeros((len(t),3))
with open("data/h_data.txt", 'w') as f:
    f.write('#t hx hy hz\n')
    for i in range(len(t)):
        e_1 = np.array([e1x[i], e1y[i], e1z[i]])
        e_2 = np.array([e2x[i], e2y[i], e2z[i]])
        # numerator = np.subtract(e_2, e_1)
        numerator = e_2 - e_1
        denominator = np.sqrt(np.dot(numerator, numerator))
        print('d=', denominator)
        # denominator = 1
        # print(t, "n=", numerator, 'd=' , denominator)
        h = numerator / denominator
        head_vector[i,0]= h[0]
        head_vector[i,1]= h[1]
        head_vector[i,2]= h[2]
        f.write("{} {:.4f} {:.4f} {:.4f}\n".format(t[i], f_vector[i,0], f_vector[i,1], f_vector[i,2]))

#getting the parallel displacement and fitting  
disp_para = np.zeros((len(t),3))
with open('data/disp_para.txt', 'w') as f:
    f.write('#t disp para\n')
    for i in range(len(t)):
        f_vector = np.array([f_vector[i,0], f_vector[i,1], f_vector[i,2]])
        d = np.array([disp_vector[i,0], disp_vector[i,1], disp_vector[i,2]])
        p = np.dot(head, d)
        f.write("{} {} {}\n".format(t[i], p, p**2)) 

def f(x, m, c):
    return m*x + c


t, p, p2 = np.loadtxt('data/disp_para.txt', unpack=True)

t_max = len(t)
sample_window_fraction = 0.02
sample_window = int(t_max * sample_window_fraction)
print("Total iterations: {}".format(t_max))
print("Sample window: {} iterations".format(sample_window))

p_sq_avg_sampled = np.zeros(sample_window)
t0_iter_list = np.arange(0, t_max-sample_window)
print("Samples taken: {}".format(len(t0_iter_list)))

t_shortened = t[:sample_window]

for t0_i in t0_iter_list:
    init_angle = p[t0_i]
    p_sq_avg_sampled+= (p[t0_i:t0_i+sample_window]-init_angle)**2

p_sq_avg_sampled /= len(t0_iter_list)
Dimension = 1

fit_params, fit_cov = curve_fit(f, t_shortened, p_sq_avg_sampled)
print('Fitted parameters: m = {:.4e}, c = {:.4e}'.format(fit_params[0], fit_params[1]))
print('Diffusion coefficient: D = {:.4e}'.format(fit_params[0]/(2*Dimension)))
print('or, approximately D = {:.4f}'.format(fit_params[0]/(2*Dimension)))

fitline = f(t_shortened, *fit_params)

plt.figure(tight_layout=True)
plt.plot(t_shortened, p_sq_avg_sampled, 'k-', label='Data')
plt.plot(t_shortened, fitline, 'r--', label='Fit')
plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel(r'$\langle MSD_{parallel}^2 \rangle$', fontsize=18)
plt.title(r"$D_{{para}} = {:.4f}$".format(fit_params[0]/(2*Dimension)))
plt.legend(fontsize=14)
plt.savefig('plots/D_paralle_cl.pdf')


