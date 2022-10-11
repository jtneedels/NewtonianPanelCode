import numpy as np

def compute_geometric_properties(simplex_verts):
    """Computes simplex area, normal vector, and centroid from vertices.

    Args:
        simplex_verts array(float): array of simplex vertices.

    Returns:
        area (float): area of simplex.
        normal array(float): outward oriented surface normal vector.
        centroid array(float): position of centroid of simplex in Cartesian coordinates.

    """
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
    """Computes aerodynamic force coefficients using a Modified Newtonian
    local surface inclination method.

    Args:
        points array(float): array of vertices of triangulated mesh.
        tri array(float): array of triangles from triangulated mesh.
        aoa_vec list(float): list of angle of attack values in radians to be examined.
        Ma_vec list(float): list of Mach number values to be examined.
        gam (float): specific heat ratio.

    Returns:
        cd_vec list(float): list of drag coefficient values.
        cl_vec list(float): list of lift coefficient values.
        cp_vec list(float): list of pressure coefficient values.

    """    
    print('Running Modified Newtonian Panel Code...')
    
    cd_vec, cl_vec = [], []
    
    for aoa in aoa_vec:
        for Ma in Ma_vec:
        
            Cp_max = 2 / (gam * M ** 2) * ((((gam + 1)**2 * M**2) / (4*gam*M**2 - 2*(gam-1)))**(gam/(gam-1))*((1 - gam + 2*gam*M**2)/(gam+1))-1)
            cd, cl, cp = compute_pressure_forces(points, tri, aoa, Cp_max)
            cd_vec.append(cd)
            cl_vec.append(cl)
            cp_vec.append(cp)
            print('Ma =',round(Ma,2),'\t aoa=',round(np.degrees(aoa),2),'\t gamma=',round(gam,2),'\t Cd=',round(cd,4), '\t Cl=',round(cl,4))
    
    print('Complete.')
    return cd_vec, cl_vec, cp_vec

def newtonian_driver(points, tri, aoa_vec=[0]):
    """Computes aerodynamic force coefficients using a Newtonian
    local surface inclination method.

    Args:
        points array(float): array of vertices of triangulated mesh.
        tri array(float): array of triangles from triangulated mesh.
        aoa_vec list(float): list of angle of attack values in radians to be examined.

    Returns:
        cd_vec list(float): list of drag coefficient values.
        cl_vec list(float): list of lift coefficient values.
        cp_vec list(list(float)): list of pressure coefficient values.

    """    
    print('Running Newtonian Panel Code...')
    
    cd_vec, cl_vec, cp_vec = [], [], []
    
    for aoa in aoa_vec:
        
        Cp_max = 2   
        cd, cl, cp = compute_pressure_forces(points, tri, aoa, Cp_max)
        cd_vec.append(cd)
        cl_vec.append(cl)
        cp_vec.append(cp)
        print('aoa=',round(np.degrees(aoa),2),'\t Cd=',round(cd,4), '\t Cl=',round(cl,4))
        
    print('Complete.')
    return cd_vec, cl_vec, cp_vec

def compute_pressure_forces(points, tri, aoa, Cp_max):
    """Compute pressure forces for each simplex using appropriate
    maximum pressure value specified by driver.

    Args:
        points array(float): array of vertices of triangulated mesh.
        tri array(float): array of triangles from triangulated mesh.
        aoa (float): angle of attack (radians).
        Cp_max (float): Cp_max coefficient calculated in driver.

    Returns:
        cd (float): drag coefficient values.
        cl (float): lift coefficient values.
        cp list(float): list of pressure coefficient values.

    """        
    cp, centroids = [], []  
    l, d, proj_area = 0, 0, 0 

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