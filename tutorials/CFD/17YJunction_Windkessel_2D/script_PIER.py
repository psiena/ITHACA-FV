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
#print(modPlift.shape)
#print(modP.shape)
coeffP = red_coeff_mat.red_coeff[:,4:,0]
print('coeff P', coeffP[:,1:].shape)
print('modP', modP.T.shape)
print('Plift',modPlift.shape)

modP_conlift = np.zeros([N_cells, 3])
modP_conlift[:,0] = modPlift
modP_conlift[:,0] = modP[:,0]
modP_conlift[:,0] = modP[:,1]
print(modP_conlift.shape)

P_rec = coeffP[:,:] @ modP_conlift.T #    modP[:, :Np].T
print('p ROM',P_rec[25,:])

P_FOM = read_scal_FOM('/scratch/psiena/cardio_rom/17YJunction_Windkessel/ITHACAoutput/Offline', 'p', 51)
print('P FOM',P_FOM.T[25,:])

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
