import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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
    scatter, = ax.plot(positions[0,:,0], positions[0,:,1], linestyle="None", marker="o", markersize=5, color="blue")

    def update(frame):
        """Update the scatter plot for each frame."""
        scatter.set_xdata(positions[frame,:,0])
        scatter.set_ydata(positions[frame,:,1])
        return scatter,

    ani = FuncAnimation(fig, update, frames=no_frames, interval=50)
    plt.show()

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
    no_particles = 100
    no_steps = 100
    L = 10
    v0 = 0.1
    dt = 0.1
    initial_positions = np.random.random(size=(no_particles,2))*L
    initial_velocities = np.random.normal(size=(no_particles,2))*v0

    positions, velocities = integrate(initial_positions, initial_velocities, L, dt, no_steps)
    animate_particles(positions, velocities, L, dt=0.1, no_frames=100)
    

