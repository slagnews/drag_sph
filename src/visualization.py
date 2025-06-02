import struct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from src import sph_cpp

import struct
import numpy as np

def read_all_output(filename):
    positions = []
    velocities = []
    with open(filename, "rb") as f:
        while True:
            int_bytes = f.read(4)
            if not int_bytes:
                break  # EOF
            N = struct.unpack("i", int_bytes)[0]

            # Read positions and velocities
            x_bytes = f.read(8 * N)
            y_bytes = f.read(8 * N)
            vx_bytes = f.read(8 * N)
            vy_bytes = f.read(8 * N)

            if any(len(b) < 8 * N for b in [x_bytes, y_bytes, vx_bytes, vy_bytes]):
                print("Warning: Incomplete frame at the end of the file")
                break

            x = np.frombuffer(x_bytes, dtype=np.float64)
            y = np.frombuffer(y_bytes, dtype=np.float64)
            vx = np.frombuffer(vx_bytes, dtype=np.float64)
            vy = np.frombuffer(vy_bytes, dtype=np.float64)

            pos = np.stack((x, y), axis=1)     # shape: (N, 2)
            vel = np.stack((vx, vy), axis=1)   # shape: (N, 2)

            positions.append(pos)
            velocities.append(vel)

    # shape: (timesteps, N, 2)
    return np.array(positions), np.array(velocities)


def kernel(q, sigma):
    w = np.zeros_like(q)
    w[q<=2] = sigma*0.25*(2-q[q<=2])**3
    w[q<=1] = sigma*(1-1.5*q[q<=1]**2+0.75*q[q<=1]**3)
    return w

def get_density_grid(pos_x, pos_y, Lx, Ly, masses, sigma, grid_size=100):
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
    
    # Generate grid
    x = np.linspace(0,Lx,grid_size)
    y = np.linspace(0,Ly,grid_size)
    X,Y = np.meshgrid(x,y)

    # Apply periodic boundary conditions to distances
    dx = (X[:, :, None] - pos_x)
    dx = (dx+Lx/2)%Lx - Lx/2
    dy = (Y[:, :, None] - pos_y)
    dy = (dy+Ly/2)%Ly - Ly/2

    distances = np.sqrt(dx**2+dy**2)

    # Calculate density
    density = np.sum(masses[0]*kernel(distances, sigma), axis=2)
    return x,y,density

def animate(positions, params, masses):
    fig, ax = plt.subplots()
    ax.set_xlim(0,float(params.Lx))
    ax.set_ylim(0,float(params.Ly))
    ax.set_aspect("equal")

    # Plot particles
    particles, = ax.plot(positions[0,:,0], positions[0,:,1], linestyle="None", marker="o", color="red")

    # Plot central object
    phi = np.linspace(0,2*np.pi,100)
    ax.plot(params.Lx/2.0+params.central_radius*np.cos(phi), params.Ly/2.0+params.central_radius*np.sin(phi), color="white")

    # Plot density
    x,y,rho = get_density_grid(positions[0,:,0], positions[0,:,1], params.Lx, params.Ly, masses, params.sigma)
    density = ax.pcolor(x,y,rho)

    def update(frame):
        particles.set_xdata(positions[frame,:,0])
        particles.set_ydata(positions[frame,:,1])
        x,y,rho = get_density_grid(positions[frame,:,0], positions[frame,:,1], params.Lx, params.Ly, masses, params.sigma)
        density.set_array(rho)
        print(f"frame {frame+1} of {len(positions)} done!", end='\r')

    anim = FuncAnimation(fig, update, len(positions))
    anim.save("figures/cpp_anim.mp4")
    plt.show()