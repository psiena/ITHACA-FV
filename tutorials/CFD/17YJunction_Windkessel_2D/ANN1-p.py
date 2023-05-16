#import plotly.graph_objects as go
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import red_coeff_mat
import errL2P_mat
import errL2U_mat
import errL2P_proj_mat
import errL2U_proj_mat

legend_elements_P = [Line2D([0], [0], color='tab:blue', lw=1, label='FOM-ROM error P'),
        Line2D([0], [0], color='tab:orange', lw=1, label='Projection error P'),
                   ]
legend_elements_U = [Line2D([0], [0], color='tab:blue', lw=1, label='FOM-ROM error U'),
        Line2D([0], [0], color='tab:orange', lw=1, label='Projection error U'),
                   ]
legend_elements = [Line2D([0], [0], color='tab:blue', lw=1, label='Red coeff'),
        Line2D([0], [0], color='tab:orange', lw=1, label='Proj coeff'),
                   ]
legend_elements0 = [Line2D([0], [0],color='b', lw=1, label='Loss'),]
legend_elements1 = [Line2D([0], [0], marker='o',color='b', lw=1, label='Train data'),
		   Line2D([0], [0], marker='o', color='r', lw=1, label='Predicted data'),
                   ]


#proj_coeff_u = np.load('/scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/POD/Proj_coeff_U.npy')
#proj_coeff_p = np.load('/scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/POD/Proj_coeff_P.npy')
#print(proj_coeff_u.shape)#[0,:])

#time = red_coeff_mat.red_coeff[:,0,0]
#red_coeff_u1 = red_coeff_mat.red_coeff[:,2,0]
#red_coeff_u2 = red_coeff_mat.red_coeff[:,3,0]
#red_coeff_p1 = red_coeff_mat.red_coeff[:,7,0]
#red_coeff_p2 = red_coeff_mat.red_coeff[:,8,0]

#print(red_coeff_mat.red_coeff.shape)

time = np.linspace(0,1,51)

plt.figure(1)
plt.plot(time, errL2P_mat.errL2P[:])
plt.plot(time,errL2P_proj_mat.errL2P_proj[:])
#plt.plot(np.linspace(0,len(proj_coeff_u[0,:]),len(proj_coeff_u[0,:])),proj_coeff_u[0,:])
plt.legend(handles=legend_elements_P)
plt.xlabel('Time Normalized')
plt.ylabel('Relative Error Pressure')
plt.grid()
plt.savefig('/u/p/psiena/Downloads/rel_err_p.png',format='png', bbox_inches='tight', pad_inches = 0)

plt.figure(2)
plt.plot(time, errL2U_mat.errL2U[:])
plt.plot(time, errL2U_proj_mat.errL2U_proj[:])
#plt.plot(np.linspace(0,len(proj_coeff_u[0,:]),len(proj_coeff_u[0,:])),proj_coeff_u[0,:])
plt.legend(handles=legend_elements_U)
plt.xlabel('Time Normalized')
plt.ylabel('Relative Error Velocity')
plt.grid()
plt.savefig('/u/p/psiena/Downloads/rel_err_u.png',format='png', bbox_inches='tight', pad_inches = 0)


#plt.figure(2)
#plt.plot(red_coeff_u2)
#plt.plot(proj_coeff_u[1,:])
#plt.plot(np.linspace(0,len(proj_coeff_u[0,:]),len(proj_coeff_u[0,:])),proj_coeff_u[0,:])
#plt.legend(handles=legend_elements)
#plt.xlabel('Time')
#plt.ylabel('$C_U$')
#plt.grid()
#plt.savefig('/u/p/psiena/Downloads/Coeff_U_2_W',format='pdf', bbox_inches='tight', pad_inches = 0)

#plt.figure(3)
#plt.plot(red_coeff_p2)
#plt.plot(proj_coeff_p[1,:])
#plt.plot(np.linspace(0,len(proj_coeff_u[0,:]),len(proj_coeff_u[0,:])),proj_coeff_u[0,:])
#plt.legend(handles=legend_elements)
#plt.xlabel('Time')
#plt.ylabel('$C_P$')
#plt.grid()
#plt.savefig('/u/p/psiena/Downloads/Coeff_P_1_W',format='pdf', bbox_inches='tight', pad_inches = 0)

