def read_all_output(filename):
    positions = []
    velocities = []
    types_all = []
    drag_forces = []

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

            drag_force_bytes = f.read(8)
            if len(drag_force_bytes) < 8:
                print("Warning: Missing drag_force at end of file")
                break
            drag_force = struct.unpack("d", drag_force_bytes)[0]

            pos = np.stack((x, y), axis=1)
            vel = np.stack((vx, vy), axis=1)

            positions.append(pos)
            velocities.append(vel)
            types_all.append(types)
            drag_forces.append(drag_force)

    return positions, velocities, types_all, drag_forces