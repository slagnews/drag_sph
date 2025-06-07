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
    types_all = []

    with open(filename, "rb") as f:
        while True:
            int_bytes = f.read(4)
            if not int_bytes:
                break
            N = struct.unpack("i", int_bytes)[0]

            x = np.frombuffer(f.read(8 * N), dtype=np.float64)
            y = np.frombuffer(f.read(8 * N), dtype=np.float64)
            vx = np.frombuffer(f.read(8 * N), dtype=np.float64)
            vy = np.frombuffer(f.read(8 * N), dtype=np.float64)
            types = np.frombuffer(f.read(N), dtype=np.uint8)

            if any(len(arr) < N for arr in [x, y, vx, vy, types]):
                print("Warning: Incomplete frame at end of file")
                break

            pos = np.stack((x, y), axis=1)
            vel = np.stack((vx, vy), axis=1)

            positions.append(pos)
            velocities.append(vel)
            types_all.append(types)

    return positions, velocities, types_all



def kernel(q, sigma):
    w = np.zeros_like(q)
    w[q<=2] = sigma*0.25*(2-q[q<=2])**3
    w[q<=1] = sigma*(1-1.5*q[q<=1]**2+0.75*q[q<=1]**3)
    return w

def relative_positions(positions, Lx, Ly):
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
    rel_positions[:,:,0] = (rel_positions[:,:,0] + Lx/2)%Lx - Lx/2
    rel_positions[:,:,1] = (rel_positions[:,:,1] + Ly/2)%Ly - Ly/2

    # Calculate distances from normalized distances:
    rel_distances = np.sqrt(np.sum(rel_positions**2, axis=2))

    return rel_positions, rel_distances

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

def get_velocity_grid(positions, velocities, masses, params, grid_size=100):
    """Interpolated the velocities of particles on a grid
    Arguments:
        positions (np.ndarray): the positions of all particles
        velocities (np.ndarray): velocities of all particles

    """
    _, distances = relative_positions(positions, params.Lx, params.Ly)
    densities = get_densities(distances, masses, params.h)

    x = np.linspace(0, params.Lx, grid_size)
    y = np.linspace(0, params.Ly, grid_size)
    X, Y = np.meshgrid(x, y)

    # Differences between grid and particle positions
    dx = (X[:, :, None] - positions[:, 0])
    dx = (dx + params.Lx / 2) % params.Lx - params.Lx / 2
    dy = (Y[:, :, None] - positions[:, 1])
    dy = (dy + params.Ly / 2) % params.Ly - params.Ly / 2

    distances = np.sqrt(dx**2 + dy**2)
    W = kernel(distances, params.h)  # shape: (grid_size, grid_size, N)

    # Broadcast and weight velocities: (N,2) → (grid, grid, N, 2)
    velocity_contribs = (masses / densities)[:, None] * velocities  # shape (N, 2)
    Vx = np.sum(W * velocity_contribs[:, 0], axis=2)
    Vy = np.sum(W * velocity_contribs[:, 1], axis=2)
    V = np.sqrt(Vx**2+Vy**2)

    return X, Y, Vx, Vy, V

def animate(positions, velocities, params, masses, field_type="density", interval=100):
    fig, ax = plt.subplots()
    ax.set_xlim(0, float(params.Lx))
    ax.set_ylim(0, float(params.Ly))
    ax.set_aspect("equal")

    # Plot particles
    particles, = ax.plot(positions[0][:,0], positions[0][:,1], linestyle="None", marker=".", color="red")

    # Plot in- and outflow boundaries
    ax.axvline(params.inflow_factor*params.Lx, color="black", linestyle="--")
    ax.axvline((1-params.outflow_factor)*params.Lx, color="black", linestyle="--")

    # Plot central object
    phi = np.linspace(0, 2*np.pi, 100)
    ax.plot(params.Lx/2.0 + params.central_radius*np.cos(phi),
            params.Ly/2.0 + params.central_radius*np.sin(phi),
            color="black")

    # Plot field (using imshow for simplicity and speed)
    if field_type == "density":
        _, _, field_data = get_density_grid(positions[0][:,0], positions[0][:,1], params.Lx, params.Ly, masses, params.sigma)
        field_img = ax.imshow(field_data, extent=[0, params.Lx, 0, params.Ly], origin="lower", cmap="viridis")
    elif field_type == "velocity":
        _, _, _, _, field_data = get_velocity_grid(positions[0], velocities[0], masses, params)
        field_img = ax.imshow(field_data, extent=[0, params.Lx, 0, params.Ly], origin="lower", cmap="viridis")

    
    #fig.colorbar(field_img, ax=ax)

    def update(frame):
        particles.set_xdata(positions[frame][:,0])
        particles.set_ydata(positions[frame][:,1])

        if field_type == "density":
            _, _, field_data = get_density_grid(positions[frame][:,0], positions[frame][:,1], params.Lx, params.Ly, masses, params.sigma)
            field_img.set_data(field_data)
        elif field_type == "velocity":
            _, _, _, _, field_data = get_velocity_grid(positions[frame], velocities[frame], masses, params)
            field_img.set_data(field_data)
        print(f"frame {frame+1} of {len(positions)} done, {len(positions[frame])} particles in sim", end='\r')
        return particles

    anim = FuncAnimation(fig, update, frames=len(positions), interval=interval)
    return anim

def animate_types(positions, velocities, types, params, masses, interval=100):
    fig, ax = plt.subplots()
    ax.set_xlim(0, float(params.Lx))
    ax.set_ylim(0, float(params.Ly))
    ax.set_aspect("equal")

    # Plot particles
    inflow_particles, = ax.plot(positions[0][:,0][types[0]==0.0], positions[0][:,1][types[0]==0.0], linestyle="None", marker=".", color="C0")
    mainflow_particles, = ax.plot(positions[0][:,0][types[0]==2.0], positions[0][:,1][types[0]==2.0], linestyle="None", marker=".", color="C1")
    outflow_particles, = ax.plot(positions[0][:,0][types[0]==1.0], positions[0][:,1][types[0]==1.0], linestyle="None", marker=".", color="C2")
    ghost_particles, = ax.plot(positions[0][:,0][types[0]==3.0], positions[0][:,1][types[0]==3.0], linestyle="None", marker=".", color="black")

    # Plot in- and outflow boundaries
    ax.axvline(params.inflow_factor*params.Lx, color="black", linestyle="--")
    ax.axvline((1-params.outflow_factor)*params.Lx, color="black", linestyle="--")

    # Plot central object
    phi = np.linspace(0, 2*np.pi, 100)
    ax.plot(params.Lx/2.0 + params.central_radius*np.cos(phi),
            params.Ly/2.0 + params.central_radius*np.sin(phi),
            color="black")
    
    def update(frame):
        inflow_particles.set_xdata(positions[frame][:,0][types[frame]==0.0])
        inflow_particles.set_ydata(positions[frame][:,1][types[frame]==0.0])

        mainflow_particles.set_xdata(positions[frame][:,0][types[frame]==2.0])
        mainflow_particles.set_ydata(positions[frame][:,1][types[frame]==2.0])

        outflow_particles.set_xdata(positions[frame][:,0][types[frame]==1.0])
        outflow_particles.set_ydata(positions[frame][:,1][types[frame]==1.0])

        ghost_particles.set_xdata(positions[frame][:,0][types[frame]==3.0])
        ghost_particles.set_ydata(positions[frame][:,1][types[frame]==3.0])

        print(f"Frame {frame+1} done", end="\r")

        return inflow_particles, mainflow_particles, outflow_particles, ghost_particles

    anim = FuncAnimation(fig, update, frames=len(positions), interval=interval)
    return anim

def calculate_energies(positions, velocities, params, masses):
    # Total mass
    M = np.sum(masses)

    no_steps = len(positions)
    E_bulk = np.zeros(no_steps)
    E_rand = np.zeros(no_steps)
    E_kinetic = np.zeros(no_steps)
    E_total = np.zeros(no_steps)
    E_object = np.zeros(no_steps)

    for t in range(no_steps):
        v = velocities[t]
        
        momentum = np.sum(masses[:, None] * v, axis=0)
        v_mean = momentum / M

        # Energy from bulk flow
        E_bulk[t] = 0.5 * M * np.dot(v_mean, v_mean)

        # Energy from object potential
        distances = np.sqrt((positions[t][:,0]-params.Lx/2)**2+(positions[t][:,1]-params.Ly/2)**2)
        E_object[t] = np.sum(masses*params.max_force*params.boundary_width*np.log(np.exp((params.central_radius-distances)/params.boundary_width)+1))

        # Total kinetic energy
        kinetic = 0.5 * masses * np.sum(v**2, axis=1)
        E_kinetic[t] = np.sum(kinetic)
        
        # Energy from random motion
        E_rand[t] = E_kinetic[t] - E_bulk[t]

        # Total energy including potential
        E_total[t] = E_kinetic[t] + E_object[t]

    return E_total, E_bulk, E_rand, E_object