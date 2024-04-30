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
t,e1x, e1y, e1z, e2x, e2y, e2z = np.loadtxt("e2e_data.txt", unpack=True)
f_vector = np.zeros((len(t),3))
with open("h_data.txt", 'w') as f:
    f.write('#t hx hy hz\n')
    for i in range(len(t)):
        e_1 = np.array([e1x[i], e1y[i], e1z[i]])
        e_2 = np.array([e2x[i], e2y[i], e2z[i]])
        numerator = e_2 - e_1
        denominator = np.sqrt(np.dot(numerator, numerator))
        # print('d=', denominator)
        h = numerator / denominator
        f_vector[i,0]= h[0]
        f_vector[i,1]= h[1]
        f_vector[i,2]= h[2]
        f.write("{} {:.4f} {:.4f} {:.4f}\n".format(t[i], f_vector[i,0], f_vector[i,1], f_vector[i,2]))

#getting the parallel displacement   
disp_para = np.zeros((len(t),3))
with open('disp_para.txt', 'w') as f:
    f.write('#t disp para\n')
    for i in range(len(t)):
        h_vector = np.array([f_vector[i,0], f_vector[i,1], f_vector[i,2]])
        d = np.array([disp_vector[i,0], disp_vector[i,1], disp_vector[i,2]])
        p = np.dot(h_vector, d)
        f.write("{} {} {}\n".format(t[i], p, p**2)) 


#getting the perp displacement   

with open('disp_perp.txt', 'w') as k:
    k.write('#t disp perp\n')
    for i in range(len(t)):
        h_vector = np.array([f_vector[i,0], f_vector[i,1], f_vector[i,2]])
        hx = h_vector[0]
        hy = h_vector[1]
        hz = h_vector[2]

        if hx != 0 :
            f = np.array([-(hy + hz)/hx,1 ,1]) 
        elif hx == 0 and hy != 0 :
            f = np.array([1, -(hz)/hy ,1])
        elif hx == 0 and hy == 0 :
            f = np.array([1,0 ,0])

        f = f/np.linalg.norm(f)


        f_prime = np.cross(h_vector,f)/np.linalg.norm(np.cross(h_vector, f))

        dy = np.dot(f, disp_vector[i] )

        dz = np.dot(f_prime, disp_vector[i])

      

        k.write("{} {} {}\n".format(t[i], dy, dz)) 



#windowing average 
t, dy ,dz = np.loadtxt('disp_perp.txt', unpack=True)

t_max = len(t)
sample_window_fraction = 0.02
sampling_window = int(t_max * sample_window_fraction)

t0_iter_list = np.arange(0, t_max-sampling_window)
t_shortened = t[:sampling_window]

disp_sq_avg_sampled= np.zeros(len(t_shortened))
for t0_i in t0_iter_list:
    init_y = dy[t0_i]
    init_z = dz[t0_i]
    disp_sq_avg_sampled += (dy[t0_i:t0_i+sampling_window] - init_y)**2 + (dz[t0_i:t0_i+sampling_window] - init_z)**2

disp_sq_avg_sampled /= len(t0_iter_list)


with open("window_displacement_perp.dat", "w") as file:
    for t_i in range(len(t_shortened)):
        file.write("{}\t{}\n".format(t_shortened[t_i], disp_sq_avg_sampled[t_i]))


#getting D parallel 
def f(x, m, c):
    return m*x + c

Dimension = 2
fit_params, fit_cov = curve_fit(f, t_shortened, disp_sq_avg_sampled)
print('Fitted parameters: m = {:.4e}, c = {:.4e}'.format(fit_params[0], fit_params[1]))
print('Diffusion coefficient: D = {:.4e}'.format(fit_params[0]/(2*Dimension)))
print('or, approximately D = {:.4f}'.format(fit_params[0]/(2*Dimension)))

fitline = f(t_shortened, *fit_params)

p, D_0, D_rot_0, D_perp, D_par, D_rot = np.loadtxt("D_info.txt", unpack=True) 
plt.figure(tight_layout=True)
plt.plot(t_shortened, disp_sq_avg_sampled, 'k-', label='Data')
plt.plot(t_shortened, fitline, 'r--', label='Fit')
plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel(r'$\langle MSD_{parallel}^2 \rangle$', fontsize=18)
plt.title(r"$D_{{perp}} = {:.4e}$".format(D_0*fit_params[0]/(2*Dimension)))
plt.legend(fontsize=14)
plt.savefig('D_perp_c.pdf')

# #getting D total
# t, p_sq = np.loadtxt("window_displacement.dat", unpack=True)

# Dimension= 2

# def f(x, m, c):
#     return m * x + c

# popt, pcov = curve_fit(f, t, p_sq)

# D_measured = popt[0] / (2.0*Dimension)

# p, D_0, D_rot_0, D_perp, D_par, D_rot = np.loadtxt("D_info.txt", unpack=True) 
# D_PA = (D_para*D_0)
# print('D p = {:.4e}'.format(D_PA))

# D_actual = (D_measured*D_0)
# print('D total = {:.4e}'.format(D_actual))
 
# print(' D perp = {:.4e}'.format(D_actual - D_PA))

 

