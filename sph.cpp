#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <random>
#include <numeric>

constexpr int neighbor_count = 9;
constexpr int neighbor_offsets[neighbor_count][2] = {
	{-1, -1}, {-1, 0}, {-1, 1},
	{0, -1},  {0, 0},  {0, -1},
	{1, -1},  {1, 0},  {1, 1}
};


struct SimulationParams {
	double res, Lx, Ly, kappa, h, rc, v0, dt, c_s, gamma_index, rho0, sigma, B;
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
		sigma = 10.0/(7.0*M_PI*h*h);
		B = c_s*c_s*rho0/gamma_index;
	}
};

struct ParticleList {
	// Particle information
	int no_particles;
	const SimulationParams& params;
	
	// Particle information vectors
	std::vector<double> pos_x, pos_y, vel_x, vel_y, acc_x, acc_y, mass, rho, pressure;
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
		  acc_x(p.init_particles),
		  acc_y(p.init_particles),
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
	
	std::array<int, 2> scalar_index_to_vec(int cell) const {
		int cx = cell % params.nc[0];
		int cy = cell / params.nc[0];
		return {cx, cy};
	}
	
	inline double kernel(double q) {
		if (q >= 2) return 0.0;
		
		double q2 = q*q;
		double q3 = q2*q;
		
		if (q>=1.0)
			return params.sigma*0.25*(2.0-q)*(2.0-q)*(2.0-q);
		else
			return params.sigma*(1.0 -1.5*q2 +0.75*q3);
	}
	
	inline double cole_pressure(double rho) {
		return params.B*(pow(rho/params.rho0, params.gamma_index)-1);
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
		cell_idx = reorder(cell_idx, indices);
	}
	
	void compute_rho_p_a() {
		std::fill(rho.begin(), rho.end(), 0);
		double h = params.h;
		double rc = params.rc;
		double h2 = h*h;
		double rc2 = rc*rc;
		
		for (int cell=0; cell<params.no_cells; cell++) {
			int start_i = cell_start[cell];
			int end_i = cell_start[cell+1];
			
			for (int i = start_i; i<end_i; ++i) {
			double xi = pos_x[i];
			double yi = pos_y[i];
			double rhoi = 0.0;
			
			auto [cx, cy] = scalar_index_to_vec(cell);
			
				for (auto [dx,dy] : neighbor_offsets) {
					int ncx = cx + dx;
					int ncy = cy + dy;
					
					if (ncx<0 || ncx>=params.nc[0] || ncy<0 || ncy>=params.nc[1]) continue;
					
					int neighbor = vec_to_scalar_index(ncx, ncy);
					int start_j = cell_start[neighbor];
					int end_j = cell_start[neighbor+1];
					
					for (int j = start_j; j<end_j; ++j) {
						double rel_x = xi-pos_x[j];
						double rel_y = yi-pos_y[j];
						double r2 = rel_x*rel_x + rel_y*rel_y;
						
						if (r2 < rc2) {
							double q = std::sqrt(r2 / h2);
							rhoi += mass[j] * kernel(q);
						}
					}
				}
				rho[i] = rhoi;
				pressure[i] = cole_pressure(rhoi);
			}
		}
	}
	void compute_density_pressure() {
		
		std::fill(rho.end(), rho.end(), 0.0);
		
		double h = params.h;
		double rc = params.rc;
		double h2 = h*h;
		
		for (int i=0; i<no_particles; ++i) {
			double xi = pos_x[i];
			double yi = pos_y[i];
			double rhoi = 0.0;
			
			auto cell = get_cell_index(xi, yi);
			int cx = cell[0];
			int cy = cell[1];
			
			for (int dx = -1; dx <= 1; ++dx) {
				for (int dy = -1; dy <= 1; ++dy) {
					int ncx = cx + dx;
					int ncy = cy + dy;
					
					if (ncx<0||ncx>=params.nc[0]||ncy<0||ncy>=params.nc[1])
						continue;
						
					int neighbor_cell = vec_to_scalar_index(ncx, ncy);
					int start = cell_start[neighbor_cell];
					int end = cell_start[neighbor_cell+1];
					
					for (int j = start; j<end; ++j) {
						double rel_x = xi - pos_x[j];
						double rel_y = yi - pos_y[j];
						double r2 = rel_x*rel_x + rel_y*rel_y;
						
						if (r2 < rc*rc) {
							double q = std::sqrt(r2/h2);
							
							double w = kernel(q);
							rhoi += mass[j]*w;
						}
					}
				}
			}
			rho[i] = rhoi;
			pressure[i] = cole_pressure(rhoi);
		}
	}
};



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

	pl.assign_cell();
	pl.get_sorter();
	pl.get_cell_start();
	pl.compart();
	pl.compute_rho_p_a();

	std::cout << "rho[1] = " << pl.rho[1] << std::endl;
	std::cout << "pressure[1] = " << pl.pressure[1] << std::endl;
		
	return 0;
}
