import numpy as np
import os
import red_coeff_mat

N_cells = 42000

def read_foam_scal(file):
    '''
    Reads a FOAM file with a scalar field inside
    Returns a vector (numpy array) of the field
    '''
    fileobj = open(file, 'r')
    words = fileobj.read().splitlines()
    words = words[23:N_cells+23]
    for i, w in enumerate(words):
        w2 = float(w)
        words[i] = w2
    words = np.array(words).T
    return words.reshape(-1, 1).flatten()


def read_foam_vec(file):
    '''
    Reads a FOAM file with a veector field inside
    Returns a vector (numpy array) of the field
    '''
    fileobj = open(file, 'r')
    words = fileobj.read().splitlines()
    words = words[23:N_cells+23]
    for i, w in enumerate(words):
        w2 = w[1:-1]
        w3 = [float(elem) for elem in list(w2.split(' '))]
        words[i] = w3
    words = np.array(words)
    return words#.reshape(-1, 1).flatten()



def read_vec(folder, field, N_modes):
    '''
    Read all vectors in POD folders 'POD/1/U', 'POD/2/U',... if 'U' is the field
    Returns an array which groups all the modes
    ONLY FOR VECTOR FIELDS
    '''
    modes = np.zeros((N_cells, 3, N_modes))
    for mode in range(1, N_modes+1, 1):
        file = os.path.join(folder, str(int(mode)), field)
        vec_mode = read_foam_vec(file)
        modes[:, :, int(mode)-1] = vec_mode
    return modes


def read_scal(folder, field, N_modes):
    '''
    Reads all vectors in POD folders 'POD/1/p', 'POD/2/p',... if 'U' is the field
    Returns an array which groups all the modes
    ONLY FOR SCALAR FIELDS
    '''
    modes = np.zeros((N_cells, N_modes))
    for mode in range(1, N_modes+1, 1):
        file = os.path.join(folder, str(int(mode)), field)
        vec_mode = read_foam_scal(file)
        modes[:, int(mode)-1] = vec_mode
    return modes

def read_scal_FOM(folder, field, N_modes):
    '''
    Reads all vectors in POD folders 'POD/1/p', 'POD/2/p',... if 'U' is the field
    Returns an array which groups all the modes
    ONLY FOR SCALAR FIELDS
    '''
    modes = np.zeros((N_cells, N_modes))
    for mode in range(2, N_modes+1, 1):
        file = os.path.join(folder, str(int(mode)), field)
        vec_mode = read_foam_scal(file)
        modes[:, int(mode)-1] = vec_mode
    return modes


folder_modes = '/scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/POD'
Nu = 2
Np = 2
modU = read_vec(folder_modes, 'U', Nu)
modP = read_scal(folder_modes, 'P', Np)
modPlift = read_foam_scal('/scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/Offline/0/Plift0')

modP_conlift = np.zeros([N_cells, 3])
modP_conlift[:,0] = modPlift#/np.linalg.norm(modPlift)
modP_conlift[:,1] = modP[:,0]#/np.linalg.norm(modP[:,0])
modP_conlift[:,2] = modP[:,1]#/np.linalg.norm(modP[:,1])

# leggo snapshots
P_FOM = read_scal_FOM('/scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/Offline', 'p', 51)

# inserisco BC pressione
t_i = 0.0
t_f = 0.5
dt = 0.01
nb_time_steps = (t_f - t_i)/dt + 1
counter = 0
time = t_i
g = np.zeros([int(nb_time_steps)])
while time < t_f + t_f * dt:
    g[counter] = 6*np.sin(6.283185*time)*np.sin(6.283185*time)+2.5*(1-np.exp(-12.566371*time)-np.sin(12.566371*time))
    time = time + dt
    counter = counter + 1

# leggo Plift
Plift = read_foam_scal('/scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/Offline/0/Plift0')

# omogenizzo snap
Pom = np.zeros([N_cells, int(nb_time_steps)])
for i in range(int(nb_time_steps)):
    Pom[:,i] = P_FOM[:,i]-g[i]*Plift[:]

Pom_FOM = read_scal_FOM('/scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/Pomfield', 'Pomfield', 51)

#print('diff',np.linalg.norm(Pom-Pom_FOM)/np.linalg.norm(Pom_FOM))
u, s, vh = np.linalg.svd(Pom, full_matrices=True)
print(u.shape)
hh
#print('modP shape',modP_conlift)

#print('Pom shape',Pom_FOM)

Coeff_proj = modP.T @ Pom_FOM #np.dot(np.dot(modP_conlift,modP_conlift.T),P_FOM)
print(Coeff_proj)
hh
P_proj = modP_conlift @ Coeff_proj

err_proj_P = np.linalg.norm(P_proj - P_FOM,axis = 1)/np.linalg.norm(P_FOM, axis = 1)

print(err_proj_P)

hh

#Err = np.linalg.norm(P_rec-P_FOM.T, axis = 1)/np.linalg.norm(P_FOM.T, axis = 1)
#print('Err',Err)

pRec_ITHACA = read_scal('/scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/Reconstruction', 'pRec', 51)
print('pRec Ithaca shape',pRec_ITHACA.shape)
print('p rec Ithaca',pRec_ITHACA.T[25,:])

names = P_rec[25,:]

# open file in write mode
with open(r'/u/p/psiena/Documents/pRec_dame.txt', 'w') as fp:
    for item in names:
        # write each item on a new line
        fp.write("%f\n" % item)