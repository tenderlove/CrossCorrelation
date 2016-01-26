def checkXY(xxx, yyy, Nx, Ny):
    '''
    x, y arrays should be between 0 and N?-1. This is broken specifically
    when I use catalogs generated from separate masks for Herschel
    auto-correlation.
    '''
    butteryBiscuit = ((xxx >= 0) & (yyy >= 0) & (xxx < Nx) & (yyy < Ny))
    xxx = xxx[butteryBiscuit]
    yyy = yyy[butteryBiscuit]
    return (xxx, yyy)
