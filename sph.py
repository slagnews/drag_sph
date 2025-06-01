import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from datetime import datetime

np.random.seed(0)

def kernel(r,h):
    """Compute the quintic spline kernel function for SPH
    Arguments:
        r (float or np.ndarray): distance
        h (float): smoothing length.
    Returns:
        (float or np.ndarray): kernel value.
    """
    return 7/(478*np.pi*h**2)*(
        np.maximum(3-r/h, 0)**5
        - 6*np.maximum(2-r/h,0)**5
        + 15*np.maximum(1-r/h,0)**5)

def deriv_kernel(r,h):
    """Compute the derivative of the quintic spline kernel with respect to r
    Arguments:
        r (float or np.ndarray): distance
        h (float): smoothing length.
    Returns:
        (float or np.ndarray): derivative value.
    """
    return 7/(478*np.pi*h**2)*(
        - 5*np.maximum(3-r/h,0)**4
        + 30*np.maximum(2-r/h,0)**4
        - 75*np.maximum(1-r/h,0)**4)

def relative_positions(positions, L):
    """Calculates relative positions and distances between particles.
    Arguments:
        positions (np.ndarray): positions of all particles
        L (float): size of the simulation
    Returns:
        rel_positions (np.ndarray): relative positions of particles
        rel_distances (np.ndarray): distance between particles
    """
    
    # Calculate NxNxD array for N atoms in D dimensions storing relative positions:
    rel_positions = positions[np.newaxis,:,:] - positions[:,np.newaxis,:]

    # Normalize relative positions according to periodic boundary conditions:
    rel_positions = (rel_positions + L/2)%L - L/2

    # Calculate distances from normalized distances:
    rel_distances = np.sqrt(np.sum(rel_positions**2, axis=2))

    return rel_positions, rel_distances

def get_densities(positions, masses, L, h):
    """Calculates the densities at the positions of all particles.
    Arguments:
        positions (np.ndarray): positions of all particles
        masses (np.ndarray):  masses of all particles
        L (float): size of the simulation
        h (float) smoothing length
    Returns:
        (np.ndarray): densities at particle positions
    """
    rel_pos, rel_dist = relative_positions(positions, L)
    return np.sum(masses*kernel(rel_dist, h), axis=1)

def get_density_grid(positions, masses, L, h, grid_size=100):
    """Calculates the densities on a grid
    Arguments:
        positions (np.ndarray): positions of all particles
        masses (np.ndarray): masses of all particles
        h (float): smoothing length
        L (float): size of the simulation
        grid_size (int): size of the grid for density calculation.
    Returns:
        x (np.ndarray): 1D array of x coordinates of the grid
        y (np.ndarray): 1D array of y coordinates of the grid.
        density (np.ndarray): 2D array of shape (grid_size, grid_size) containing the density values at each grid point.
    """
    # Genarate grid
    x = np.linspace(0,L,grid_size)
    y = np.linspace(0,L,grid_size)
    X,Y = np.meshgrid(x,y)

    # Apply periodic boundary conditions to distances
    dx = (X[:,:,None]-positions[:,0])
    dx = (dx+L/2)%L - L/2
    #dy = dx - L*np.round(dx/L)
    dy = (Y[:,:,None]-positions[:,1])#%L
    dy = (dy+L/2)%L - L/2
    #dy = dy - L * np.round(dy / L)

    distances = np.sqrt(dx**2+dy**2)

    # Calculate density
    density = np.sum(masses*kernel(distances, h), axis=2)
    return x,y,density

def cole_pressure(rho, c_s, rho0, gamma):
    """Calculates the pressure from the Cole Equation of State
    Arguments:
        rho (np.ndarray or float): density
        c_s (float): speed of sound
        rho0 (float): reference density
        gamma (float): adiabatic index
    Returns:
        (np.ndarray or float): pressure
    """
    B = c_s**2*rho0/gamma
    return B*((rho/rho0)**gamma-1)

def get_pressures(rho, c_s, rho0):
    """Calculates the pressures at all particle positions
    Arguments:
        rho (np.ndarray): densities at particle positions
        c_s (float): speed of sound
        rho0 (float): reference density
        gamma (float): adiabatic index
    Returns:
        (np.ndarray): pressures at particle positions
    """
    pressures = cole_pressure(rho, c_s, rho0, 7)
    return pressures

def get_pressure_grid(density_grid, c_s, rho0):
    """Calculates the pressures on a grid
    Arguments:
        rho (np.ndarray): density grid
        c_s (float): speed of sound
        rho0 (float): reference density
        gamma (float): adiabatic inde
    Returns:
        (np.ndarray): pressure grid
    """
    pressure_grid = cole_pressure(density_grid, c_s, rho0, 7)
    return pressure_grid

def get_pressure_forces(positions, masses, L, h, c_s, rho0):
    """Calculates the forces from pressure on all particles
    Arguments:
        positions (np.ndarray): positions of all particles
        masses (np.ndarray): masses of all particles
        L (float): size of the simulation
        h (float): smoothing length
        c_s (float): speed of sound
        rho0 (float): reference density
    Returns:
        accelerations (np.ndarray): accelerations of all particles
    """
    rel_pos, distances = relative_positions(positions, L)
    forces = np.zeros_like(positions)
    densities = get_densities(positions, masses, L, h)
    pressures = get_pressures(densities, c_s, rho0)

    distances2 = np.where(distances==0, np.inf, distances)
    nablaW = deriv_kernel(distances, h)[:,:,None]*-rel_pos/distances2[:,:,None]
    for i in range(len(positions)):
        forces[i] = -np.sum((masses*(pressures/densities**2+pressures[i]/densities[i]**2))[:,None]*nablaW[i],axis=0)*masses[i]
    return forces

def central_object_forces(positions, L, central_object):
        """Calculates the force of the central object on all particles
        Arguments:
            positions (np.ndarray): positions of all particles
            L (float): size of the simulation
            central_object (dict): the central object, holding radius, boundary_width and max_force properties
        Returns:
            forces (np.ndarray): the forces of the object on all particles
        """
        dx = positions[:,0] - L/2
        dx = (dx+L/2)%L - L/2
        dy = positions[:,1] - L/2
        dy = (dy+L/2)%L - L/2
        distance = np.sqrt(dx**2+dy**2)
        nx = dx/distance
        ny = dy/distance
        n = np.column_stack((nx,ny))

        strength = 1/(1+np.exp((distance-central_object["radius"])/central_object["boundary_width"]))
        
        forces = central_object["max_force"]*strength[:,None]*n
        return forces


def integrate(initial_positions, initial_velocities, masses, L, h, c_s, rho0, dt, no_steps, central_object=None):
    """Integrate the positions and velocities of particles in a 2D space.
    Arguments:
        initial_positions: 2D numpy array of shape (no_particles, 2) containing the initial x and y coordinates of the particles.
        initial_velocities: 2D numpy array of shape (no_particles, 2) containing the initial x and y velocities of the particles.
        L: Length of the box (for periodic boundary conditions).
        dt: Time step for integration.
        no_steps: Number of time steps to integrate.
        central_object (dict): properties of an optional central object
    Returns:
        positions: 3D numpy array of shape (no_steps, no_particles, 2) containing the x and y coordinates of the particles at each time step.
        velocities: 3D numpy array of shape (no_steps, no_particles, 2) containing the x and y velocities of the particles at each time step.
    """
    start = datetime.now()

    # Initialize data arrays
    no_particles = initial_positions.shape[0]
    positions = np.zeros((no_steps, no_particles, 2))
    velocities = np.zeros((no_steps, no_particles, 2))
    #kinetic_energies = np.zeros((no_steps, no_particles))
    #internal_energies = np.zeros((no_steps, no_particles))

    # Set initial conditions
    positions[0] = initial_positions
    velocities[0] = initial_velocities
    #kinetic_energies[0] = 1/2*np.sum(velocities[0]**2, axis=1)
    
    for i in range(1, no_steps):
        positions[i] = positions[i-1] + velocities[i-1] * dt
        positions[i] %= L

        if central_object:
            velocities[i] = velocities[i-1] + dt*(get_pressure_forces(positions[i], masses, L, h, c_s, rho0)+central_object_forces(positions[i], L, central_object))/masses[:,None]
        else:
            velocities[i] = velocities[i-1] + dt*get_pressure_forces(positions[i], masses, L, h, c_s, rho0)/masses
        """
        # Calculate kinetic energies
        kinetic_energies[i] = 1/2*np.sum(velocities[i]**2, axis=1)

        # Calculate internal energies
        
        densities = get_densities(positions[i], masses, L, h)
        pressures = get_pressures(densities, c_s, rho0)

        distances, rel_pos = relative_positions(positions[i], L)
        distances2 = np.where(distances==0, np.inf, distances) 
        nablaW = deriv_kernel(distances, h)[:,:,None]*-rel_pos/distances2[:,:,None]
        for j in range(no_particles):
            dv = velocities[j] - velocities
            
            dudt = pressures[j]/densities[j]**2+np.sum(masses*(dv[0]*nablaW[j,:,0]+dv[1]*nablaW[j,:,0]))
            internal_energies[i,j] += dudt*dt
        """
        ttg = (datetime.now()-start)/i * (no_steps-i)
        print(f"{(i+1)/no_steps*100:.0f}% done, time to go: {ttg}", end="\r")
    return positions, velocities


def animate_particles(positions, masses, L, h, fps=30, central_object=None, file_name="animation.mp4"):
    """Animates the positions of particles in a 2D space.
    Arguments:
        positions (np.ndarray): positions of all particles
        masses (np.ndarray): masses of all particles
        L (float): size of the simulation
        h (float): smoothing length
        object (dict): optional central object
        file_name (float): file to which to save the animation
    """
    

    no_frames = len(positions)

    fig, ax = plt.subplots()
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_aspect("equal")
    scatter, = ax.plot(positions[0,:,0], positions[0,:,1], linestyle="None", marker="o", markersize=1, color="red")
    x,y,rho = get_density_grid(positions[0], masses, L, h)
    mesh = ax.pcolor(x,y,rho, cmap="viridis")
    if object:
        theta = np.linspace(0,2*np.pi,100)
        ax.plot(L/2+central_object["radius"]*np.cos(theta),L/2+central_object["radius"]*np.sin(theta), color="white")
    start = datetime.now()
    def update(frame):
        """Update the scatter plot for each frame."""
        
        scatter.set_xdata(positions[frame,:,0])
        scatter.set_ydata(positions[frame,:,1])
        x,y,rho = get_density_grid(positions[frame], masses, L, h)
        mesh.set_array(rho.ravel())

        ttg = (datetime.now() - start)/(frame+1) * (no_frames - frame)
        print(f"Frame {frame+1} of {no_frames} done, time to go: {ttg}", end="\r")
        return scatter,mesh

    ani = FuncAnimation(fig, update, frames=no_frames, interval=50)
    plt.show()
    ani.save(file_name, fps=fps, extra_args=['-vcodec', 'libx264'], dpi=300)
    print(f"Animation saved to {file_name}!")