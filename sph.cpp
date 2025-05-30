#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <random>
#include <numeric>

struct SimulationParams {
	double res, Lx, Ly, kappa, h, rc, v0, dt, c_s, gamma_index, rho0;
	int no_cells, no_steps, Nx, Ny, init_particles;
	std::array<int, 2> nc;
	
	SimulationParams(double res, int no_steps, double Lx, double Ly, double rho0, double kappa, double v0, double dt, double c_s, double gamma_index)
		: res(res),
		  no_steps(no_steps),
		  Lx(Lx),
		  Ly(Ly),
		  rho0(rho0),
		  kappa(kappa),
		  v0(v0),
		  dt(dt),
		  c_s(c_s),
		  gamma_index(gamma_index)
	{
		h = kappa*res;
		rc = 2*h;
		nc = { int(std::ceil(Lx / rc)), int(std::ceil(Ly / rc)) };
		no_cells = nc[0] * nc[1];
		Nx = int(Lx/res);
		Ny = int(Ly/res);
		init_particles = Nx*Ny;
	}
};

struct ParticleList {
	// Particle information
	int no_particles;
	const SimulationParams& params;
	
	// Particle information vectors
	std::vector<double> pos_x, pos_y, vel_x, vel_y, mass, rho, pressure;
	std::vector<int> cell_idx, indices, cell_counts, cell_start;
	
	// Randomizer generator and distribution
	std::mt19937 engine;
	std::uniform_real_distribution<double> init_x_dist;
	std::uniform_real_distribution<double> init_y_dist;
	
	// Constructor
	ParticleList(const SimulationParams& p)
		: params(p),
		  no_particles(p.init_particles),
		  pos_x(p.init_particles),
		  pos_y(p.init_particles),
		  vel_x(p.init_particles),
		  vel_y(p.init_particles),
		  mass(p.init_particles),
		  rho(p.init_particles),
		  pressure(p.init_particles),
		  cell_idx(p.init_particles),
		  indices(p.init_particles),
		  cell_counts(p.no_cells),
		  cell_start(p.no_cells+1),
		  engine(std::random_device{}()),
		  init_x_dist(0, p.Lx),
		  init_y_dist(0, p.Ly) {}
		  
	// Helper functions
	std::array<int, 2> get_cell_index(double x, double y) {
		return { int(x / params.rc), int(y / params.rc) };
	}

	int vec_to_scalar_index(int cx, int cy) {
		return cx + params.nc[0] * cy;
	}
	
	void init_random() {
	// Randomize the x and y positions over the whole domain
		for (size_t i=0; i<no_particles; ++i) {
			pos_x[i] = init_x_dist(engine);
			pos_y[i] = init_y_dist(engine);
		}
	}
	
	void init_grid() {
		int index = 0;
		for (size_t i=0; i<params.Nx; ++i) {
			for (size_t j=0; j<params.Ny; ++j) {
				pos_x[index] = i*params.res;
				pos_y[index] = j*params.res;
				++index;
			}
		}
	}
	
	void init_mass() {
		for (int i=0; i<no_particles; ++i) {
			mass[i] = params.rho0*params.res*params.res;
		}
	}
	
	void init_vel() {
		for (int i=0; i<no_particles; ++i) {
			vel_x[i] = params.v0;
			vel_y[i] = 0;
		}
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
		atof(argv[1]), // res
		atoi(argv[2]), // no_steps
		atof(argv[3]), // Lx
		atof(argv[4]), // Ly
		atof(argv[5]), // rho0
		atof(argv[6]), // kappa
		atof(argv[7]), // v0
		atof(argv[8]), // dt
		atof(argv[9]), // c_s
		atof(argv[10]) // gamma_index
	);
	
	ParticleList pl(params);
	
	pl.init_grid();
	pl.init_mass();
	pl.init_vel();
	std::cout << pl.vel_x[0] << std::endl;
	
	return 0;
}
