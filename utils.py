import numpy as np
from scipy.stats import gaussian_kde


# gently taken from https://stackoverflow.com/questions/33793701/pyplot-scatter-to-contour-plot
def density_estimation(m1, m2):
    X, Y = np.mgrid[m1.min():m1.max():100j, m2.min():m2.max():100j]
    
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    
    kernel = gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    
    return X, Y, Z