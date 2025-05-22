import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

np.random.seed(0)

def kernel(r,h):
    """Compute the kernel function for SPH.
    Arguments:
        r: distance
        h: Smoothing length.
    Returns:
        Kernel value.
    """
    return 1/(2*np.pi*h**2)*np.exp(-0.5*(r/h)**2)

def deriv_kernel(r,h):
    return -r/(2*np.pi*h**4)*np.exp(-0.5*(r/h)**2)

def relative_positions(positions, L):
    """Calculates relative positions and distances between particles.
    Arguments:
        positions (np.ndarray): The positions of the particles in cartesian space
        L (float): The size of the simulation space
    Returns:
        rel_positions (np.ndarray): Relative positions of particles
        rel_distances (np.ndarray): The distance between particles
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
        positions (np.ndarray): The particle positions in carthesian space
        masses (np.ndarray):  The particle masses
        L (float): the size of the simulation space
        h (float) the kernel size
    Returns:
        (np.ndarray): densities at particle positions
    """
    rel_pos, rel_dist = relative_positions(positions, L)
    return np.sum(masses*kernel(rel_dist, h), axis=1)

def get_density_grid(positions, masses, L, h, grid_size=100):
    """Calculates the densities on a grid
    Arguments:
        positions: 2D numpy array of shape (no_particles, 2) containing the x and y coordinates of the particles.
        h: Smoothing length.
        L: Length of the box (for periodic boundary conditions).
        grid_size: Size of the grid for density calculation.
    Returns:
        x: 1D numpy array of x coordinates of the grid.
        y: 1D numpy array of y coordinates of the grid.
        density: 2D numpy array of shape (grid_size, grid_size) containing the density values at each grid point.
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

def tait_pressure(rho, c_s, rho0, gamma):
    B = c_s**2*rho0/gamma
    return B*((rho/rho0)**gamma-1)

def get_pressures(densities, c_s, rho0):
    """Calculates the pressures at all particle positions"""
    pressures = tait_pressure(densities, c_s, rho0, 7)
    return pressures

def get_pressure_grid(density_grid, c_s, rho0):
    """Calculates the pressures on a grid"""
    pressure_grid = tait_pressure(density_grid, c_s, rho0, 7)
    return pressure_grid

def accelerations(positions, masses, L, h, c_s, rho0):
    rel_pos, distances = relative_positions(positions, L)
    accelerations = np.zeros_like(positions)
    densities = get_densities(positions, masses, L, h)
    pressures = get_pressures(densities, c_s, rho0)

    distances2 = np.where(distances==0, np.inf, distances)
    nablaW = deriv_kernel(distances, h)[:,:,None]*-rel_pos/distances2[:,:,None]
    for i in range(len(positions)):
        accelerations[i] = -np.sum((masses*(pressures/densities**2+pressures[i]/densities[i]**2))[:,None]*nablaW[i],axis=0)
    return accelerations

def animate_particles(positions, masses, L, h, no_frames, file_name="animation.mp4"):
    """Animate the positions of particles in a 2D space.
    Arguments:
        positions: 3D numpy array of shape (no_steps, no_particles, 2) containing the x and y coordinates of the particles.
        L: Length of the box (for setting limits).
        no_frames: Number of frames to display in the animation.
    """

    fig, ax = plt.subplots()
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_aspect("equal")
    scatter, = ax.plot(positions[0,:,0], positions[0,:,1], linestyle="None", marker="o", markersize=1, color="red")
    x,y,rho = get_density_grid(positions[0], masses, L, h)
    mesh = ax.pcolor(x,y,rho, cmap="viridis")

    def update(frame):
        """Update the scatter plot for each frame."""
        print(f"Frame {frame} of {no_frames} done", end="\r")
        scatter.set_xdata(positions[frame,:,0])
        scatter.set_ydata(positions[frame,:,1])
        x,y,rho = get_density_grid(positions[frame], masses, L, h)
        mesh.set_array(rho.ravel())
        return scatter,mesh

    ani = FuncAnimation(fig, update, frames=no_frames, interval=50)
    plt.show()
    ani.save(file_name, fps=30, extra_args=['-vcodec', 'libx264'], dpi=300)
    print(f"Animation saved to {file_name}!")

def integrate(initial_positions, initial_velocities, masses, L, h, c_s, rho0, dt, no_steps):
    """Integrate the positions and velocities of particles in a 2D space.
    Arguments:
        initial_positions: 2D numpy array of shape (no_particles, 2) containing the initial x and y coordinates of the particles.
        initial_velocities: 2D numpy array of shape (no_particles, 2) containing the initial x and y velocities of the particles.
        L: Length of the box (for periodic boundary conditions).
        dt: Time step for integration.
        no_steps: Number of time steps to integrate.
    Returns:
        positions: 3D numpy array of shape (no_steps, no_particles, 2) containing the x and y coordinates of the particles at each time step.
        velocities: 3D numpy array of shape (no_steps, no_particles, 2) containing the x and y velocities of the particles at each time step.
    """
    positions = np.zeros((no_steps, initial_positions.shape[0], 2))
    velocities = np.zeros((no_steps, initial_positions.shape[0], 2))
    positions[0] = initial_positions
    velocities[0] = initial_velocities

    for i in range(1, no_steps):
        positions[i] = positions[i-1] + velocities[i-1] * dt
        positions[i] %= L  # Apply periodic boundary conditions
        velocities[i] = velocities[i-1] + dt*accelerations(positions[i], masses, L, h, c_s, rho0) # Assuming constant velocity for simplicity
        print(f"{(i+1)/no_steps*100:.0f}% done", end="\r")
    return positions, velocities