def filterCatalog(xx, yy, Nx, Ny, nanmask=1.0):
    '''
    I'll be using catalogs with sources outside the Herschel area (mask),
    so inject them into a map, apply a mask, and return
    the sources that remain!
    '''
    xymap = np.zeros((Nx, Ny))
    xymap[xx,yy] = 1.0
    xymap *= nanmask
    remainders = np.where(xymap == 1)
    xcat, ycat = remainders[0], remainders[1]
    return (xcat, ycat)
