import numpy as np
def randPosDiscrete(Nx, Ny, Ngal):
    ''' 
    Given a number of galaxies and the size of the array,
    distribute them in a random, integral way.
    '''
    randx = np.random.random_integers(0, Nx-1, Ngal)
    randy = np.random.random_integers(0, Ny-1, Ngal)
    return (randx, randy)
