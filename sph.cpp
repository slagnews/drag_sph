#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <numeric>

struct SimulationParams {
	double Lx, Ly, h, rc, v0, dt, c_s, gamma_index, rho0;
	int no_particles, no_cells, no_steps;
	std::array<int, 2> nc;
	
	SimulationParams(int no_steps, double Lx, double Ly, double rho0, double h, double v0, double dt, double c_s, double gamma_index)
		: no_steps(no_steps),
		  Lx(Lx),
		  Ly(Ly),
		  rho0(rho0),
		  h(h),
		  v0(v0),
		  dt(dt),
		  c_s(c_s),
		  gamma_index(gamma_index)
	{
		rc = 2*h;
		nc = {int(Lx/rc), int(Ly/rc)};
		no_cells = nc[0] * nc[1];
	}
};

struct ParticleList {
	// Particle information
	int no_particles;
	const SimulationParams& params;
	
	// Particle information vectors
	std::vector<double> pos_x, pos_y, mass, rho, pressure;
	std::vector<int> cell_idx, indices, cell_counts, cell_start;
	
	// Constructor
	ParticleList(int no_particles,  const SimulationParams& p)
		: params(p),
		  no_particles(no_particles),
		  pos_x(no_particles),
		  pos_y(no_particles),
		  mass(no_particles),
		  rho(no_particles),
		  pressure(no_particles),
		  cell_idx(no_particles),
		  indices(no_particles),
		  cell_counts(p.no_cells),
		  cell_start(p.no_cells+1) {}
		  
	// Helper functions
	std::array<int, 2> get_cell_index(double x, double y) {
		return { int(x / params.rc), int(y / params.rc) };
	}

	int vec_to_scalar_index(int cx, int cy) {
		return cx * params.nc[1] + cy;
	}
	
	void assign_cell() {
	// Assign the corresponding cell to each particle
		for (int i = 0; i < no_particles; ++i) {
			auto cell = get_cell_index(pos_x[i], pos_y[i]);
			int cx = cell[0];
			int cy = cell[1];
			int c = vec_to_scalar_index(cx, cy);
			cell_idx[i] = c;
		}
	}

	void get_sorter() {
	// Get the indices needed to sort particle vectors
		std::iota(indices.begin(), indices.end(), 0);
		
		std::stable_sort(indices.begin(), indices.end(),
		 [&](int a, int b) {
		 	return cell_idx[a] < cell_idx[b];
		 });
	}
	
	void get_cell_start() {
		std::fill(cell_counts.begin(), cell_counts.end(), 0);
		
		for (int i=0; i<no_particles; i++) {
			int cell = cell_idx[indices[i]];
			cell_counts[cell]++;
		}
		
		std::fill(cell_start.begin(), cell_start.end(), 0);
		for (int c=0; c < params.no_cells; c++) {
			cell_start[c+1] = cell_start[c] + cell_counts[c];
		}
	}

	template<typename T>
	std::vector<T> reorder(const std::vector<T>& input, const std::vector<int>& indices) {
		std::vector<T> output(indices.size());
		for (size_t i = 0; i < indices.size(); ++i) {
			output[i] = input[indices[i]];
		}
		return output;
	}


	void compart() {
	// Does all of the above steps, call this after every step to (re)compartmentalize the data.
		assign_cell();
		
		get_sorter();
		
		get_cell_start();
		
		pos_x = reorder(pos_x, indices);
		pos_y = reorder(pos_y, indices);
		mass = reorder(mass, indices);
		rho = reorder(rho, indices);
		pressure = reorder(pressure, indices);
	}
};

double cole_pressure(double rho, double c_s, double rho0, double gamma_index) {
	double B = pow(c_s, 2)*rho0/gamma_index;
	return B*(pow(rho/rho0, gamma_index)-1);
}


int main(int argc, char *argv[]) {
	if (argc < 11) {
		std::cerr << "Not enough arguments given" << std::endl;
		return 1;
	}
	
	SimulationParams params(
		atoi(argv[2]), // no_steps
		atof(argv[3]), // Lx
		atof(argv[4]), // Ly
		atof(argv[5]), // rho0
		atof(argv[6]), // h
		atof(argv[7]), // v0
		atof(argv[8]), // dt
		atof(argv[9]), // c_s
		atof(argv[10]) // gamma_index
	);
	
	ParticleList pl(atoi(argv[1]), params); // Where argv[1] = no_particles
	
	return 0;
}
