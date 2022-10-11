import numpy as np

def compute_geometric_properties(simplex_verts):
    
    # compute x and y positions of simplex vertices
    x1, x2, x3 = simplex_verts[0][0], simplex_verts[1][0], simplex_verts[2][0]
    y1, y2, y3 = simplex_verts[0][1], simplex_verts[1][1], simplex_verts[2][1]
    z1, z2, z3 = simplex_verts[0][2], simplex_verts[1][2], simplex_verts[2][2]
    
    # compute length arrays 
    xa, xb = x2 - x1, x3 - x1
    ya, yb = y2 - y1, y3 - y1
    za, zb = z2 - z1, z3 - z1
    
    a = np.array([xa, ya, za])
    b = np.array([xb, yb, zb])
    
    # compute area
    area = 0.5 * np.linalg.norm(np.cross(a,b))
    
    # compute surface unit normal vector
    # depending on surface convention, need to check that normal is
    # is outward by default
    normal = -np.cross(a,b) / (np.linalg.norm(np.cross(a,b)) + 1e-12)    
    
    # compute simplex centroid coordinates
    centroid  = np.zeros((3,))
    
    centroid[0] = (x1 + x2 + x3)/3
    centroid[1] = (y1 + y2 + y3)/3
    centroid[2] = (z1 + z2 + z3)/3    
    
    return area, normal, centroid

def modified_newtonian_driver(points, tri, aoa_vec=[0], Ma_vec=[10], gam=1.4):
    
    print('Running Modified Newtonian Panel Code')
    
    cd_vec, cl_vec = [], []
    
    for aoa in aoa_vec:
        for Ma in Ma_vec:
            print('Parameters: ','Ma =',Ma,' aoa=',np.degrees(aoa))
        
            Cp_max = 2 / (gam * M ** 2) * ((((gam + 1)**2 * M**2) / (4*gam*M**2 - 2*(gam-1)))**(gam/(gam-1))*((1 - gam + 2*gam*M**2)/(gam+1))-1)
            cd, cl, cp = compute_pressure_forces(points, tri, aoa, Cp_max)
            cd_vec.append(cd)
            cl_vec.append(cl)
            cp_vec.append(cp)
    
    print('Complete.')
    return cd_vec, cl_vec, cp_vec

def newtonian_driver(points, tri, aoa_vec=[0]):
    
    print('Running Newtonian Panel Code')
    
    cd_vec, cl_vec, cp_vec = [], [], []
    
    for aoa in aoa_vec:
        print('Parameters: ','aoa=',np.degrees(aoa))
        
        Cp_max = 2   
        cd, cl, cp = compute_pressure_forces(points, tri, aoa, Cp_max)
        cd_vec.append(cd)
        cl_vec.append(cl)
        cp_vec.append(cp)
        
    print('Complete.')
    return cd_vec, cl_vec, cp_vec

def compute_pressure_forces(points, tri, aoa, Cp_max):
    
    cp = []   
    centroids = []
    proj_area = 0
    d = 0
    l = 0     

    for simplex_verts in points[tri.simplices]:

        area, normal, centroid = compute_geometric_properties(simplex_verts)
        centroids.append(centroid)

        if np.dot(np.array([-np.cos(aoa),np.sin(aoa),0]), normal) > 0:

            theta = np.arcsin(np.dot(np.array([-np.cos(aoa),np.sin(aoa),0]), normal))
            cp.append(Cp_max * np.sin(theta)**2)
            p = Cp_max * np.sin(theta)**2 * area
            d += p * np.sin(theta)       
            l += p * np.cos(theta) * np.sign(normal[1])               
            proj_area += area * np.sin(np.arcsin(np.dot(np.array([-np.cos(0),np.sin(0),0]), normal))) 

        else:
            cp.append(0)
            proj_area += area * np.sin(np.arcsin(np.dot(np.array([-np.cos(0),np.sin(0),0]), normal))) 

    cd = d / proj_area
    cl = l / proj_area

    return cd, cl, cp