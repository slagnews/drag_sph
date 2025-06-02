#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <random>
#include <numeric>
#include <iomanip>
#include <fstream>

constexpr int neighbor_count = 9;
constexpr int neighbor_offsets[neighbor_count][2] = {
	{-1, -1}, {-1, 0}, {-1, 1},
	{0, -1},  {0, 0},  {0, 1},
	{1, -1},  {1, 0},  {1, 1}
};


struct SimulationParams {
	double res, Lx, Ly, kappa, h, rc, v0, dt, c_s, gamma_index, rho0, sigma, B;
	int no_cells, no_steps, Nx, Ny, init_particles;
	std::array<int, 2> nc;
	double central_radius, boundary_width, max_force;
	SimulationParams(double res, int no_steps, double Lx, double Ly, double rho0, double kappa, double v0, double dt, double c_s, double gamma_index, double central_radius, double boundary_width, double max_force)
		: res(res),
		  no_steps(no_steps),
		  Lx(Lx),
		  Ly(Ly),
		  rho0(rho0),
		  kappa(kappa),
		  v0(v0),
		  dt(dt),
		  c_s(c_s),
		  gamma_index(gamma_index),
		  central_radius(central_radius),
		  boundary_width(boundary_width),
		  max_force(max_force)
	{
		h = kappa*res;
		rc = 3.0*h;
		nc = { int(std::ceil(Lx / rc)), int(std::ceil(Ly / rc)) };
		no_cells = nc[0] * nc[1];
		Nx = int(Lx/res);
		Ny = int(Ly/res);
		init_particles = Nx*Ny;
		sigma = 7.0/(478.0*M_PI*h*h);//10.0/(7.0*M_PI*h*h);
		B = c_s*c_s*rho0/gamma_index;
	}
};

struct ParticleList {
	// Particle information
	int no_particles;
	int current_step;
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
	
	inline double kernel(double q) const {
		if (q >= 3.0) return 0.0;
		else if (q >= 2.0) return params.sigma*pow(3.0-q,5);
		else if (q >= 2.0) return params.sigma*(pow(3.0-q,5) - 6*pow(2.0-q,5));
		else if (q >= 1.0) return params.sigma*(pow(3.0-q,5) - 6*pow(2.0-q,5) + 15*pow(1.0-q,5));
	}

	inline double deriv_kernel(double q) {
		if (q >= 3.0) return 0.0;
		else if (q >= 2.0) return -5*params.sigma*pow(3.0-q,4);
		else if (q >= 2.0) return -5*params.sigma*(pow(3.0-q,4) - 6*pow(2.0-q,4));
		else if (q >= 1.0) return -5*params.sigma*(pow(3.0-q,4) - 6*pow(2.0-q,4) + 15*pow(1.0-q,4));
	}
	
	inline double cole_pressure(double rho) {
		return params.B*(pow(rho/params.rho0, params.gamma_index)-1);
	}
	
	void init_random() {
	// Randomize the x and y positions over the whole domain
		for (int i=0; i<no_particles; ++i) {
			pos_x[i] = init_x_dist(engine);
			pos_y[i] = init_y_dist(engine);
		}
	}
	
	void init_grid() {
		int index = 0;
		for (int i=0; i<params.Nx; ++i) {
			for (int j=0; j<params.Ny; ++j) {
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
		for (int i = 0; i < indices.size(); ++i) {
			output[i] = input[indices[i]];
		}
		return output;
	}


	void compart() {
	// Does all of the above steps, call this after every step to (re)compartmentalize the data.
		assign_cell();
		get_sorter();
		get_cell_start();
		
		auto new_pos_x = reorder(pos_x, indices);
		auto new_pos_y = reorder(pos_y, indices);
		auto new_vel_x = reorder(vel_x, indices);
		auto new_vel_y = reorder(vel_y, indices);
		auto new_acc_x = reorder(acc_x, indices);
		auto new_acc_y = reorder(acc_y, indices);
		auto new_mass = reorder(mass, indices);
		auto new_rho = reorder(rho, indices);
		auto new_pressure = reorder(pressure, indices);
		auto new_cell_idx = reorder(cell_idx, indices);

		pos_x = std::move(new_pos_x);
		pos_y = std::move(new_pos_y);
		vel_x = std::move(new_vel_x);
		vel_y = std::move(new_vel_y);
		acc_x = std::move(new_acc_x);
		acc_y = std::move(new_acc_y);
		mass = std::move(new_mass);
		rho = std::move(new_rho);
		pressure = std::move(new_pressure);
		cell_idx = std::move(new_cell_idx);
	}
	
	void compute_rho_p() {
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

	void compute_a() {
		std::fill(acc_x.begin(), acc_x.end(), 0);
		std::fill(acc_y.begin(), acc_y.end(), 0);
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
			double acc_xi = 0.0;
			double acc_yi = 0.0;
			
			auto [cx, cy] = scalar_index_to_vec(cell);
			
				for (auto [dx,dy] : neighbor_offsets) {
					int ncx = cx + dx;
					int ncy = cy + dy;
					
					if (ncx<0 || ncx>=params.nc[0] || ncy<0 || ncy>=params.nc[1]) continue;
					
					int neighbor = vec_to_scalar_index(ncx, ncy);
					int start_j = cell_start[neighbor];
					int end_j = cell_start[neighbor+1];
					
					for (int j = start_j; j<end_j; ++j) {
						if (i == j) continue;

						double rel_x = xi-pos_x[j];
						double rel_y = yi-pos_y[j];
						double r2 = rel_x*rel_x + rel_y*rel_y;
						
						if (r2 < rc2 && r2 > 1e-12) {
							double r = std::sqrt(r2);
							double q = std::sqrt(r2 / h2);

							double common = -1*mass[j]*
							(pressure[i]/(rho[i]*rho[i]) + 
								pressure[j]/(rho[j]*rho[j]))
							*deriv_kernel(q)/(h*r);

							acc_xi += common * rel_x;
							acc_yi += common * rel_y;
						}
					}
				}
				acc_x[i] = acc_xi;
				acc_y[i] = acc_yi;
			}
		}
	}

	void add_central_object_forces(){
		double Lx = params.Lx;
		double Ly = params.Ly;
		double radius = params.central_radius;
		double width = params.boundary_width;
		double max_force = params.max_force;

		for (int i = 0; i < no_particles; ++i){
			double dx = pos_x[i] - Lx/2.0;
			dx = std::fmod(dx + Lx/2.0, Lx) - Lx/2.0;

			double dy = pos_y[i] - Ly/2.0;
			dy = std::fmod(dy + Ly/2.0, Ly) - Ly/2.0;

			double distance = std::sqrt(dx*dx + dy*dy);

			if (distance > 1e-8){
				double nx = dx/distance;
				double ny = dy/distance;
				
				double strength = 1.0/(1.0+std::exp((distance - radius)/width));
				
				acc_x[i] += max_force * strength * nx;
				acc_y[i] += max_force * strength * ny;
			}
		}
	}

	void compute_v_half() {
		for (int i=0; i<no_particles; ++i) {
			vel_x[i] += 0.5 * params.dt * acc_x[i];
			vel_y[i] += 0.5 * params.dt * acc_y[i];
		}
	}

	void compute_pos() {
		for (int i=0; i<no_particles; ++i) {
			pos_x[i] += params.dt * vel_x[i];
			pos_y[i] += params.dt * vel_y[i];
		}
	}

	void compute_v_full() {
		for (int i=0; i<no_particles; ++i) {
			vel_x[i] += 0.5 * params.dt * acc_x[i];
			vel_y[i] += 0.5 * params.dt * acc_y[i];
		}
	}

	void print_pos() {
		std::cout << std::fixed << std::setprecision(10);  // ← Add this
		for (int i=0; i<no_particles; ++i) {
			std::cout << current_step << "," << pos_x[i] << "," << pos_y[i] << std::endl;
		}
	}

	void write_frame(std::ofstream& out) {
		int N = pos_x.size();
		out.write(reinterpret_cast<const char*>(&N), sizeof(int));
		out.write(reinterpret_cast<const char*>(pos_x.data()), sizeof(double) * N);
		out.write(reinterpret_cast<const char*>(pos_y.data()), sizeof(double) * N);
		out.write(reinterpret_cast<const char*>(vel_x.data()), sizeof(double) * N);
		out.write(reinterpret_cast<const char*>(vel_y.data()), sizeof(double) * N);
	}


	void periodic() {
		for (int i = 0; i < no_particles; ++i) {
			pos_x[i] = std::fmod(pos_x[i], params.Lx);
			if (pos_x[i] < 0) pos_x[i] += params.Lx;

			pos_y[i] = std::fmod(pos_y[i], params.Ly);
			if (pos_y[i] < 0) pos_y[i] += params.Ly;
		}
	}

	void integrate() {
		std::ofstream out("all_output.bin", std::ios::binary);

		for (int step=0; step<params.no_steps; step++) {
			compart();
			compute_rho_p();
			compute_a();
			add_central_object_forces();
			compute_v_half();
			compute_pos();
			periodic();
			compute_v_full();
			write_frame(out);
			current_step += 1;
		}

		out.close();
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
		atof(argv[10]), // gamma_index
		atof(argv[11]), // central_radius
		atof(argv[12]), // boundary_width
		atof(argv[13])  // max_force
	);
	
	ParticleList pl(params);
	
	pl.init_grid();
	pl.init_mass();
	pl.init_vel();

	pl.integrate();
	std::cout << pl.no_particles << std::endl;
		
	return 0;
}
