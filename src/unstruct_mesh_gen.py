import numpy as np
from scipy.spatial import Delaunay

def triangulate_sphere(k=20):
    
    # domain parametrization
    U = np.linspace(0, 2 * np.pi, k)
    V = np.linspace(0, np.pi, k)
    [X, Y] = np.meshgrid(U, V)
    
    #points = np.array([X.flatten(), Y.flatten()]).T
    
    # sphere parametrization
    S1 = np.cos(X) * np.sin(Y)
    S2 = np.sin(X) * np.sin(Y)
    S3 = np.cos(Y)
    
    points = np.array([S1.flatten(), S2.flatten(), S3.flatten()]).T

    # triangulate the points in [0,2pi] x [0,pi]
    tri = Delaunay(np.array([X.flatten(), Y.flatten()]).T)
    
    return points, tri

def triangulate_spheroid(k=20):
    
    # domain parametrization
    U = np.linspace(0, np.pi, k)
    V = np.linspace(0, 2*np.pi, k)
    [X, Y] = np.meshgrid(U, V)
    
    #points = np.array([X.flatten(), Y.flatten()]).T
    
    # sphere parametrization
    S1 = np.cos(X) * np.sin(Y)
    S2 = np.sin(X) * np.sin(Y)
    S3 = 0.5*np.cos(Y)
    
    points = np.array([S1.flatten(), S2.flatten(), S3.flatten()]).T

    # triangulate the points in [0,2pi] x [0,pi]
    tri = Delaunay(np.array([X.flatten(), Y.flatten()]).T)
    
    return points, tri

def triangulate_cylinder(k=20):
    
    # domain parametrization
    U = np.linspace(0, 2 * np.pi, k)
    V = np.linspace(-1, 1, k)
    [X, Y] = np.meshgrid(U, V)
    
    #points = np.array([X.flatten(), Y.flatten()]).T
    
    # sphere parametrization
    S1 = np.cos(X) 
    S2 = np.sin(X)
    S3 = Y
    
    points = np.array([S1.flatten(), S2.flatten(), S3.flatten()]).T

    # triangulate the points in [0,2pi] x [0,pi]
    tri = Delaunay(np.array([X.flatten(), Y.flatten()]).T)

    return points, tri

def triangulate_sphere_cone(k=20):
    
    cone_angle = np.radians(60)
    r = 0.24
    height = 0.2

    # domain parametrization
    U = np.linspace(0, 2 * np.pi, k)
    V = np.linspace(0, np.pi/2 - cone_angle, k)
    W = np.linspace(r - r*np.cos(np.pi/2 - cone_angle), height, k)
    
    S1 = np.zeros((2*k*k,))
    S2 = np.zeros((2*k*k,))
    S3 = np.zeros((2*k*k,))    
    
    # spherical nose
    c = 0
    for i in range(k):
        for j in range(k):
            S1[c] = r*(1-np.cos(V[i]))
            S2[c] = r*np.sin(V[i]) * np.cos(U[j])
            S3[c] = r*np.sin(V[i]) * np.sin(U[j])         
            c += 1

    # cone
    for i in range(k):
        for j in range(k):
            S1[c] = W[i]
            S2[c] = ((S1[c]-r + r*np.cos(np.pi/2 - cone_angle))*np.tan(cone_angle) + np.sin(np.pi/2 - cone_angle) * r) * np.cos(U[j])
            S3[c] = ((S1[c]-r + r*np.cos(np.pi/2 - cone_angle))*np.tan(cone_angle) + np.sin(np.pi/2 - cone_angle) * r) * np.sin(U[j])        
            c += 1        
  
    points = np.array([S1.flatten(), S2.flatten(), S3.flatten()]).T
    
    #fig, ax =  plt.subplots(subplot_kw={"projection": "3d"})
    #ax.scatter(S1, S2, S3)
    #ax.set_box_aspect((np.ptp(S1), np.ptp(S2), np.ptp(S3)))
    #plt.show()
    #triangulate the points in [0,2pi] x [0,pi]
    tri = Delaunay(np.array([S2.flatten(), S3.flatten()]).T)
    
    return points, tri

def triangulate_wedge(k=20):
    
    # domain parametrization
    U = np.linspace(0, 1, k)
    V = np.linspace(-1, 1, k)
    [X, Y] = np.meshgrid(U, V)
    
    #points = np.array([X.flatten(), Y.flatten()]).T
    
    # sphere parametrization
    S1 = X
    S2 = Y
    S3 = -X*np.sin(np.pi/6)
    
    points = np.array([S1.flatten(), S2.flatten(), S3.flatten()]).T

    # triangulate the points in [0,2pi] x [0,pi]
    tri = Delaunay(np.array([X.flatten(), Y.flatten()]).T)

    return points, tri