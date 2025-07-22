import numpy as np 

num_particles = 38
radius = 0.03
volume = 0.000173
# Step 1: Random directions on the sphere
vecs = np.random.normal(size=(num_particles, 3))
vecs /= np.linalg.norm(vecs, axis=1)[:, np.newaxis]

# Step 2: Random radii scaled by cube root for uniform volume distribution
u = np.random.uniform(0, 1, num_particles)
radii = radius * u**(1/3)

# Step 3: Final positions
positions = vecs * radii[:, np.newaxis]
gradients = np.zeros_like(positions)  # Shape: (100, 3)



def cubic_spline_kernel(r, h):
    q = r / h
    sigma = 1 / (np.pi * h**3)
    
    if 0 <= q < 1:
        return sigma * (1 - 1.5*q**2 + 0.75*q**3)
    elif 1 <= q < 2:
        return sigma * 0.25 * (2 - q)**3
    else:
        return 0.0
    
def poly6_kernel(r, h):
    """
    Scalar or vectorized Poly6 kernel for 3D.
    r: scalar or array of distances
    h: smoothing length
    """
    sigma = 315 / (64 * np.pi * h**9)
    result = np.zeros_like(r)
    mask = (r >= 0) & (r <= h)
    result[mask] = sigma * (h**2 - r[mask]**2)**3
    return result

def spiky_kernel_derivative(r, h):
    """
    Scalar derivative dW/dr of Spiky kernel.
    r: scalar or numpy array of distances
    h: smoothing length
    """
    sigma = 15 / (np.pi * h**6)
    dW = np.zeros_like(r)
    mask = (r >= 0) & (r < h)
    dW[mask] = -3 * sigma * (h - r[mask])**2
    return dW
    
def cubic_spline_kernel_gradient(r, h):
    q = r / h
    sigma = 1 / (np.pi * h**3)

    if 0 <= q < 1:
        return sigma * (1 / h) * (-3*q + 2.25*q**2)
    elif 1 <= q < 2:
        return sigma * (1 / h) * (-0.75 * (2 - q)**2)
    else:
        return 0.0
    
sup_rad = 2.0 * radius  # Example support radius 
num_local_iter = 100
rho0 = 1000.0
for _ in range(num_local_iter):

    xi = positions[0]  # First particle position
    rho = 0.0
    grad_c_i = np.zeros(3)  # Gradient of the kernel at xi
    SchurComplement = 0.0 
    for j in range(1, num_particles):
        xj = positions[j]
        r = np.linalg.norm(xi - xj)
        n = (xi - xj) / r
        mass = volume * rho0
        rho += mass* poly6_kernel(r , sup_rad)
        grad_c_j = spiky_kernel_derivative(r, sup_rad) * n 

        gradients[j] = grad_c_j / rho0
        grad_c_i -= grad_c_j 

        SchurComplement += np.dot(gradients[j], gradients[j])

    gradients[0] = grad_c_i / rho0

    c = np.maximum(rho / rho0 - 1.0, 0.0)
    print(c)
    if c < 1e-14:
        break
    
    SchurComplement += np.dot(gradients[0], gradients[0])
    lagrangian = c / SchurComplement

    for j in range(1, num_particles):
        positions[j] += lagrangian * gradients[j]

    positions[0] += lagrangian * gradients[0]

x = 5 
y = 10 
m = 10 




        
