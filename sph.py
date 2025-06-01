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

def get_densities(relative_distances, masses, h):
    """Calculates the densities at the positions of all particles.
    Arguments:
        positions (np.ndarray): positions of all particles
        masses (np.ndarray):  masses of all particles
        L (float): size of the simulation
        h (float) smoothing length
    Returns:
        (np.ndarray): densities at particle positions
    """
    return np.sum(masses*kernel(relative_distances, h), axis=1)

def get_density_grid(positions, masses, L, h, grid_size=100, is_ghost=None):
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
    if not isinstance(is_ghost, np.ndarray):
        is_ghost = np.zeros(len(positions), dtype=bool)

    # Generate grid
    x = np.linspace(0,L,grid_size)
    y = np.linspace(0,L,grid_size)
    X,Y = np.meshgrid(x,y)

    # Apply periodic boundary conditions to distances
    dx = (X[:, :, None] - positions[~is_ghost, 0])
    dx = (dx+L/2)%L - L/2
    dy = (Y[:, :, None] - positions[~is_ghost, 1])
    dy = (dy+L/2)%L - L/2

    distances = np.sqrt(dx**2+dy**2)

    # Calculate density
    density = np.sum(masses[~is_ghost]*kernel(distances, h), axis=2)
    return x,y,density

def ideal_pressure(rho, c_s, rho0, beta=0.07):
    """Calculates the pressure from a modified ideal gas law
    Arguments:
        rho (np.ndarray): density
        c_s (float): speed of sound
        rho0 (float): reference density
        beta (float): coefficient for numerical stability
    Returns:
        (np.ndarray): pressure
    """
    return beta*c_s**2*rho0 + c_s**2*(rho-rho0)


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
    densities = get_densities(distances, masses, h)
    pressures = ideal_pressure(densities, c_s, rho0)

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


def integrate(initial_positions, initial_velocities, masses, L, h, c_s, rho0, dt, no_steps, central_object=None, is_ghost=None):
    """Integrate the positions and velocities of particles in a 2D space.
    Arguments:
        initial_positions: 2D numpy array of shape (no_particles, 2) containing the initial x and y coordinates of the particles.
        initial_velocities: 2D numpy array of shape (no_particles, 2) containing the initial x and y velocities of the particles.
        L: Length of the box (for periodic boundary conditions).
        dt: Time step for integration.
        no_steps: Number of time steps to integrate.
        central_object (dict): properties of an optional central object
        is_ghost (np.ndarray): booleans indicating for each particle if it is a ghost
    Returns:
        positions: 3D numpy array of shape (no_steps, no_particles, 2) containing the x and y coordinates of the particles at each time step.
        velocities: 3D numpy array of shape (no_steps, no_particles, 2) containing the x and y velocities of the particles at each time step.
    """
    if not isinstance(is_ghost, np.ndarray):
        is_ghost = np.array([0]*len(initial_positions))
    
    start = datetime.now()

    # Initialize data arrays
    no_particles = initial_positions.shape[0]
    positions = np.zeros((no_steps, no_particles, 2))
    velocities = np.zeros((no_steps, no_particles, 2))
    forces = np.zeros((no_steps, no_particles, 2))

    # Set initial conditions
    positions[0] = initial_positions
    velocities[0] = initial_velocities
    if central_object:
        forces[0] = get_pressure_forces(positions[0], masses, L, h, c_s, rho0) + central_object_forces(positions[0], L, central_object)
    else:
        forces[0] = get_pressure_forces(positions[0], masses, L, h, c_s, rho0)
    
    for i in range(no_steps-1):
        positions[i+1] = positions[i] + velocities[i] * dt + 1/2*forces[i]/masses[:,None]*dt**2
        positions[i+1,is_ghost] = positions[i,is_ghost]
        positions[i+1] %= L
        if central_object:
            forces[i+1] = get_pressure_forces(positions[i+1], masses, L, h, c_s, rho0) + central_object_forces(positions[i+1], L, central_object)
        else:
            forces[i+1] = get_pressure_forces(positions[i+1], masses, L, h, c_s, rho0)
        velocities[i+1] = velocities[i] + 1/2*(forces[i] + forces[i+1])*dt
        velocities[i+1,is_ghost] = velocities[i,is_ghost]
        

        
        #if central_object:
        #    #forces = get_pressure_forces(positions[i], masses, L, h, c_s, rho0)+central_object_forces(positions[i], L, central_object)
        #    velocities[i+1] = velocities[i] + dt*(get_pressure_forces(positions[i], masses, L, h, c_s, rho0)+central_object_forces(positions[i], L, central_object))/masses[:,None]
        #else:
        #    #forces = get_pressure_forces(positions[i], masses, L, h, c_s, rho0)
        #    velocities[i+1] = velocities[i] + dt*get_pressure_forces(positions[i], masses, L, h, c_s, rho0)/masses
        ttg = (datetime.now()-start)/(i+1) * (no_steps-i-1)
        print(f"{(i+1)/no_steps*100:.0f}% done, time to go: {ttg}", end="\r")
    return positions, velocities


def animate_particles(positions, masses, L, h, fps=30, central_object=None, file_name="animation.mp4", is_ghost=None):
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
    fluid_particles, = ax.plot(positions[0,~is_ghost,0], positions[0,~is_ghost,1], linestyle="None", marker="o", markersize=1, color="red")
    ghost_particles, = ax.plot(positions[0,is_ghost,0], positions[0,is_ghost,1], linestyle="None", marker="o", markersize=1, color="white")
    x,y,rho = get_density_grid(positions[0], masses, L, h, is_ghost=is_ghost)
    mesh = ax.pcolor(x,y,rho, cmap="viridis")
    if central_object:
        theta = np.linspace(0,2*np.pi,100)
        ax.plot(L/2+central_object["radius"]*np.cos(theta),L/2+central_object["radius"]*np.sin(theta), color="white")
    start = datetime.now()
    def update(frame):
        """Update the scatter plot for each frame."""
        
        fluid_particles.set_xdata(positions[frame,~is_ghost,0])
        fluid_particles.set_ydata(positions[frame,~is_ghost,1])
        ghost_particles.set_xdata(positions[frame,is_ghost,0])
        ghost_particles.set_ydata(positions[frame,is_ghost,1])
        x,y,rho = get_density_grid(positions[frame], masses, L, h, is_ghost=is_ghost)
        mesh.set_array(rho.ravel())

        ttg = (datetime.now() - start)/(frame+1) * (no_frames - frame)
        print(f"Frame {frame+1} of {no_frames} done, time to go: {ttg}", end="\r")
        return fluid_particles,ghost_particles,mesh

    ani = FuncAnimation(fig, update, frames=no_frames, interval=50)
    plt.show()
    ani.save(file_name, fps=fps, extra_args=['-vcodec', 'libx264'], dpi=300)
    print(f"Animation saved to {file_name}!")