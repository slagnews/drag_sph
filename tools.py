import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def read_all_output(filename):
    positions = []
    types_all = []
    drag_forces = []

    with open(filename, "rb") as f:
        while True:
            n_bytes = f.read(4)  # int N
            if not n_bytes:
                break
            if len(n_bytes) < 4:
                print("Warning: incomplete N")
                break
            N = struct.unpack('i', n_bytes)[0]

            x_bytes = f.read(8 * N)
            if len(x_bytes) < 8 * N:
                print("Warning: incomplete pos_x")
                break
            x = np.frombuffer(x_bytes, dtype=np.float64)

            y_bytes = f.read(8 * N)
            if len(y_bytes) < 8 * N:
                print("Warning: incomplete pos_y")
                break
            y = np.frombuffer(y_bytes, dtype=np.float64)

            types_bytes = f.read(N)  # uint8_t types, 1 byte each
            if len(types_bytes) < N:
                print("Warning: incomplete types")
                break
            types = np.frombuffer(types_bytes, dtype=np.uint8)

            drag_bytes = f.read(8)  # single double drag_force
            if len(drag_bytes) < 8:
                print("Warning: incomplete drag_force")
                break
            drag_force = struct.unpack('d', drag_bytes)[0]

            pos = np.stack((x, y), axis=1)

            positions.append(pos)
            types_all.append(types)
            drag_forces.append(drag_force)

    return positions, types_all, drag_forces

def animate_types(positions, types, interval=100):
    fig, ax = plt.subplots()

    # Plot particles
    inflow_particles, = ax.plot(positions[0][:,0][types[0]==0.0], positions[0][:,1][types[0]==0.0], linestyle="None", marker=".", color="C0")
    mainflow_particles, = ax.plot(positions[0][:,0][types[0]==2.0], positions[0][:,1][types[0]==2.0], linestyle="None", marker=".", color="C1")
    outflow_particles, = ax.plot(positions[0][:,0][types[0]==1.0], positions[0][:,1][types[0]==1.0], linestyle="None", marker=".", color="C2")
    ghost_particles, = ax.plot(positions[0][:,0][types[0]==3.0], positions[0][:,1][types[0]==3.0], linestyle="None", marker=".", color="black")

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

    anim = animation.FuncAnimation(fig, update, frames=len(positions), interval=interval)
    return anim

