#import plotly.graph_objects as go
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import red_coeff_mat
#from scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/red_coeff import red_coeff_mat

legend_elements = [Line2D([0], [0], color='blue', lw=1, label='CU_1'),
		   Line2D([0], [0], color='orange', lw=1, label='CU_2'),
                   Line2D([0], [0], color='red', lw=1, label='CU_lift'),
                   Line2D([0], [0], color='m', lw=1, label='BC_U'),
                   ]
legend_elements0 = [Line2D([0], [0],color='tab:blue', lw=1, label='BC pressure'),]
legend_elements1 = [Line2D([0], [0], marker='o',color='b', lw=1, label='Train data'),
		   Line2D([0], [0], marker='o', color='r', lw=1, label='Predicted data'),
                   ]

#coeefs = np.load('/scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/POD/coeefs_U.npy')
BC = np.load('/scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/POD/prova_BC.npy')
BC_P = np.load('/scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/POD/prova_BC_P.npy')
#proj_coef_P = np.load('/scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/POD/Proj_coeff_P.npy')
#proj_coef_U = np.load('/scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/POD/Proj_coeff_U.npy')
#proj_coeff_u = np.load('/scratch/psiena/cardio_rom/17YJunction_cilindrogrande/ITHACAoutput/POD/Proj_coeff_U.npy')
#proj_coeff_p = np.load('/scratch/psiena/cardio_rom/17YJunction_cilindrogrande/ITHACAoutput/POD/Proj_coeff_P.npy')
#prova = np.load('/scratch/psiena/cardio_rom/17YJunction_cilindrogrande/ITHACAoutput/red_coeff/red_coeff.py')
#print(proj_coeff[0,:])

#time = red_coeff_mat.red_coeff[:,0,0]
#red_coeff_ulift = red_coeff_mat.red_coeff[:,1,0]
#red_coeff_u1 = red_coeff_mat.red_coeff[:,2,0]
#red_coeff_u2 = red_coeff_mat.red_coeff[:,3,0]
#red_coeff_plift = red_coeff_mat.red_coeff[:,4,0]
#red_coeff_p1 = red_coeff_mat.red_coeff[:,5,0]
#red_coeff_p2 = red_coeff_mat.red_coeff[:,6,0]
BC1 = np.zeros(51)
BC1_P = np.zeros(51)
j=0
for i in range(0,501,10):
        BC1[j] = BC[0,i]
        BC1_P[j] = BC_P[0,i]
        j = j+1

plt.figure(0)
plt.plot(BC1_P,'tab:blue')
#plt.plot(proj_coef_P[0,:])
#plt.plot(red_coeff_u2, 'orange')
#plt.plot(red_coeff_ulift, 'red')
#33plt.plot(coeefs[0,:])
#plt.plot(BC1[:],'m')
#plt.plot(proj_coeff_u[0,:])
#plt.plot(np.linspace(0,len(proj_coeff_u[0,:]),len(proj_coeff_u[0,:])),proj_coeff_u[0,:])
plt.legend(handles=legend_elements0)
plt.xlabel('Time')
plt.ylabel('BC pressure')
plt.grid()
plt.savefig('/u/p/psiena/Downloads/BC_pressure.pdf',format='pdf', bbox_inches='tight', pad_inches = 0)
#plt.savefig('/u/p/psiena/Downloads/BC_velocity.pdf',format='pdf', bbox_inches='tight', pad_inches = 0)
hh

plt.figure(1)
plt.plot(red_coeff_p1,'blue')
plt.plot(red_coeff_p2, 'orange')
plt.plot(red_coeff_plift, 'red')
plt.plot(BC1_P[:],'m')
#plt.plot(proj_coeff_u[0,:])
#plt.plot(np.linspace(0,len(proj_coeff_u[0,:]),len(proj_coeff_u[0,:])),proj_coeff_u[0,:])
plt.legend(handles=legend_elements)
plt.xlabel('Time')
plt.ylabel('$C_P$')
plt.grid()
plt.savefig('/u/p/psiena/Downloads/Coeff_P_1',format='png', bbox_inches='tight', pad_inches = 0)

plt.figure(1)
plt.plot(red_coeff_u1)
plt.plot(proj_coeff_u[0,:])
#plt.plot(np.linspace(0,len(proj_coeff_u[0,:]),len(proj_coeff_u[0,:])),proj_coeff_u[0,:])
plt.legend(handles=legend_elements)
plt.xlabel('Time')
plt.ylabel('$C_U$')
plt.grid()
plt.savefig('/u/p/psiena/Downloads/Coeff_U_1',format='pdf', bbox_inches='tight', pad_inches = 0)

plt.figure(2)
plt.plot(red_coeff_u2)
plt.plot(proj_coeff_u[1,:])
#plt.plot(np.linspace(0,len(proj_coeff_u[0,:]),len(proj_coeff_u[0,:])),proj_coeff_u[0,:])
plt.legend(handles=legend_elements)
plt.xlabel('Time')
plt.ylabel('$C_U$')
plt.grid()
plt.savefig('/u/p/psiena/Downloads/Coeff_U_2',format='pdf', bbox_inches='tight', pad_inches = 0)

plt.figure(3)
plt.plot(red_coeff_p1)
plt.plot(proj_coeff_p[0,:])
#plt.plot(np.linspace(0,len(proj_coeff_u[0,:]),len(proj_coeff_u[0,:])),proj_coeff_u[0,:])
plt.legend(handles=legend_elements)
plt.xlabel('Time')
plt.ylabel('$C_P$')
plt.grid()
plt.savefig('/u/p/psiena/Downloads/Coeff_P_1',format='pdf', bbox_inches='tight', pad_inches = 0)

plt.figure(4)
plt.plot(red_coeff_p2)
plt.plot(proj_coeff_p[1,:])
#plt.plot(np.linspace(0,len(proj_coeff_u[0,:]),len(proj_coeff_u[0,:])),proj_coeff_u[0,:])
plt.legend(handles=legend_elements)
plt.xlabel('Time')
plt.ylabel('$C_P$')
plt.grid()
plt.savefig('/u/p/psiena/Downloads/Coeff_P_2',format='pdf', bbox_inches='tight', pad_inches = 0)


