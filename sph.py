import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def kernel(dx,h):
    """Compute the kernel function for SPH.
    Arguments:
        x: Position of the particle.
        x0: Position of the reference particle.
        h: Smoothing length.
    Returns:
        Kernel value.
    """
    return 1/np.sqrt(2*np.pi*h**2)*np.exp(-0.5*(dx/h)**2)


def density(positions, h, L, grid_size=100):
    """Calculates the density of particles on grid.
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
    dx = (X[:,:,None]-positions[:,0])#%L
    dx = dx - L * np.round(dx / L)
    dy = (Y[:,:,None]-positions[:,1])#%L
    dy = dy - L * np.round(dy / L)

    # Calculate density
    density = np.sum(kernel(np.sqrt(dx**2 + dy**2), h), axis=2)
    return x,y,density


def animate_particles(positions, L, no_frames):
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
    x,y,rho = density(positions[0], 0.5, L)
    mesh = ax.pcolor(x,y,rho, cmap="viridis")

    def update(frame):
        """Update the scatter plot for each frame."""
        print(f"frame {frame}")
        scatter.set_xdata(positions[frame,:,0])
        scatter.set_ydata(positions[frame,:,1])
        x,y,rho = density(positions[frame], 0.5, L)
        mesh.set_array(rho.ravel())
        return scatter,mesh

    ani = FuncAnimation(fig, update, frames=no_frames, interval=50)
    plt.show()
    ani.save("animation.mp4", fps=30, extra_args=['-vcodec', 'libx264'], dpi=300)

def integrate(initial_positions, initial_velocities, L, dt, no_steps):
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
        velocities[i] = velocities[i-1]  # Assuming constant velocity for simplicity

    return positions, velocities

if __name__ == "__main__":
    no_particles = 1000
    no_steps = 100
    L = 10
    v0 = 0.1
    dt = 0.1
    initial_positions = np.random.random(size=(no_particles,2))*L
    initial_velocities = np.random.normal(size=(no_particles,2))*v0

    positions, velocities = integrate(initial_positions, initial_velocities, L, dt, no_steps)
    fig, ax = plt.subplots()
    x,y,rho = density(positions[0], 0.5, L)
    ax.pcolor(x,y,rho, cmap="viridis")
    ax.scatter(positions[0,:,0], positions[0,:,1], color="red", s=1)
    plt.show()
    animate_particles(positions, L, no_frames=100)
    

