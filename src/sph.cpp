#include "../include/sph.h"

constexpr int neighbor_offsets[neighbor_count][2] = {
	{-1, -1}, {-1, 0}, {-1, 1},
	{0, -1},  {0, 0},  {0, 1},
	{1, -1},  {1, 0},  {1, 1}
};

SimulationParams::SimulationParams(double res, int no_steps, double Lx, double Ly, double rho0, double kappa, double v0, double dt,
		double c_s, double beta, double central_radius, double inflow_factor, double outflow_factor, double kinematic_viscosity, double max_beta)
	: res(res),
	  no_steps(no_steps),
	  Lx(Lx),
	  Ly(Ly),
	  rho0(rho0),
	  kappa(kappa),
	  v0(v0),
	  dt(dt),
	  c_s(c_s),
	  beta(beta),
	  central_radius(central_radius),
	  inflow_factor(inflow_factor),
	  outflow_factor(outflow_factor),
	  kinematic_viscosity(kinematic_viscosity),
	  max_beta(max_beta)
{
	h = kappa*res;
	rc = 3.0*h;
	nc = { int(std::ceil(Lx / rc)), int(std::ceil(Ly / rc)) };
	no_cells = nc[0] * nc[1];
	Nx = int(Lx/res-1);
	Ny = int(Ly/res);
	init_particles = Nx*Ny;
	sigma = 7.0/(478.0*M_PI*h*h);//10.0/(7.0*M_PI*h*h);
	P0 = beta*c_s*c_s*rho0;
}

ParticleList::ParticleList(const SimulationParams& p)
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
	  type(p.init_particles),
	  cell_counts(p.no_cells),
	  cell_start(p.no_cells+1),
	  max_a(0.0) {}



// Helper functions
std::array<int, 2> ParticleList::get_cell_index(double x, double y) const {
	return { int(x / params.rc), int(y / params.rc) };
}

int ParticleList::vec_to_scalar_index(int cx, int cy) const {
	return cx + params.nc[0] * cy;
}

std::array<int, 2> ParticleList::scalar_index_to_vec(int cell) const {
	int cx = cell % params.nc[0];
	int cy = cell / params.nc[0];
	return {cx, cy};
}

inline double ParticleList::kernel(double q) const {
	if (q >= 3.0) return 0.0;
	else if (q >= 2.0) return params.sigma*pow(3.0-q,5);
	else if (q >= 1.0) return params.sigma*(pow(3.0-q,5) - 6*pow(2.0-q,5));
	else return params.sigma*(pow(3.0-q,5) - 6*pow(2.0-q,5) + 15*pow(1.0-q,5));
}

inline double ParticleList::deriv_kernel(double q) const {
	if (q >= 3.0) return 0.0;
	else if (q >= 2.0) return -5*params.sigma*pow(3.0-q,4);
	else if (q >= 1.0) return -5*params.sigma*(pow(3.0-q,4) - 6*pow(2.0-q,4));
	else return -5*params.sigma*(pow(3.0-q,4) - 6*pow(2.0-q,4) + 15*pow(1.0-q,4));
}

inline double ParticleList::calculate_pressure(double rho) const {
	return params.P0 + params.c_s*params.c_s*(rho-params.rho0);//params.B*(pow(rho/params.rho0, params.gamma_index)-1);
}


void ParticleList::init_grid() {
	int index = 0;
	for (int i=0; i<params.Nx; ++i) {
		for (int j=0; j<params.Ny; ++j) {
			pos_x[index] = (i+1)*params.res;
			pos_y[index] = j*params.res;
			++index;
		}
	}
}

void ParticleList::init_mass() {
	#pragma omp parallel for
	for (int i=0; i<no_particles; ++i) {
		mass[i] = params.rho0*params.res*params.res;
	}
}

void ParticleList::init_vel() {
	#pragma omp parallel for
	for (int i=0; i<no_particles; ++i) {
		vel_x[i] = params.v0;
		vel_y[i] = 0;
	}
}

void ParticleList::init_type() {
	#pragma omp parallel for
	for (int i=0; i< no_particles; ++i) {
		if (pos_x[i] <= params.Lx * params.inflow_factor) {
			type[i] = ParticleType::inflow;
		}

		else if (pos_x[i] >= (params.Lx - params.Lx * params.outflow_factor)) {
			type[i] = ParticleType::outflow;
		}
		else if ((pos_x[i]-params.Lx/2)*(pos_x[i]-params.Lx/2) + (pos_y[i]-params.Ly/2)*(pos_y[i]-params.Ly/2) < params.central_radius*params.central_radius){
			type[i] = ParticleType::ghost;
		}
		else {
			type[i] = ParticleType::mainflow;
		}
	}
}

void ParticleList::assign_cell() {
// Assign the corresponding cell to each particle
	#pragma omp parallel for
	for (int i = 0; i < no_particles; ++i) {
		auto cell = get_cell_index(pos_x[i], pos_y[i]);
		int cx = cell[0];
		int cy = cell[1];
		int c = vec_to_scalar_index(cx, cy);
		cell_idx[i] = c;
	}
}

void ParticleList::manage_inflow_outflow() {
	// First we erase particles that have left the domain
	std::vector<bool> should_erase(no_particles, false);

	for (int i=0; i < no_particles; ++i) {
		if (pos_x[i] > params.Lx) {
			// Particles that leave the simulation on the right side should be removed
			should_erase[i] = true;
		}
		else if (pos_x[i] < 0.0){
			// Particles that leave the simulation on the left side should be removed
			should_erase[i] = true;
		}
	}

	int write_index = 0;
	for (int read_index=0; read_index < no_particles; ++read_index) {
		if (!should_erase[read_index]) {
			if (write_index != read_index) {
				pos_x[write_index] = pos_x[read_index];
				pos_y[write_index] = pos_y[read_index];
				vel_x[write_index] = vel_x[read_index];
				vel_y[write_index] = vel_y[read_index];
				type[write_index] = type[read_index];
			}
			++write_index;
		}
	}

	no_particles = write_index; // Update the number of particles

	// Shorten all of the lists by the new number of particles
	pos_x.resize(no_particles);
	pos_y.resize(no_particles);
	vel_x.resize(no_particles);
	vel_y.resize(no_particles);
	acc_x.resize(no_particles);
	acc_y.resize(no_particles);
	mass.resize(no_particles);
	rho.resize(no_particles);
	pressure.resize(no_particles);
	cell_idx.resize(no_particles);
	indices.resize(no_particles);
	type.resize(no_particles);


	for (int i=0; i < no_particles; ++i) {
		
		if ((pos_x[i] > params.Lx * params.inflow_factor) && (type[i] == ParticleType::inflow)) {
			// If particle went from inflow to mainflow, reassign, and init new particle in inflow
			type[i] = ParticleType::mainflow;

			pos_x.push_back(0.0);
			pos_y.push_back(pos_y[i]); // With y position of previous particle

			vel_x.push_back(params.v0);
			vel_y.push_back(0.0);

			mass.push_back(params.rho0*params.res*params.res);
			type.push_back(ParticleType::inflow);

			// The rest of these push back statements are purely to extend the vectors
			acc_x.push_back(0.0);
			acc_y.push_back(0.0);
			rho.push_back(0.0);
			pressure.push_back(0.0);
			cell_idx.push_back(0.0);
			indices.push_back(0.0);
		} 

		else if ((pos_x[i] >= (params.Lx - params.Lx * params.outflow_factor)) && (type[i] == ParticleType::mainflow)) {
			// If particles went from main flow to outflow, change type
			type[i] = ParticleType::outflow;
		}
	}

	no_particles = pos_x.size(); // Intermediate update of the number of particles
}

void ParticleList::get_sorter() {
// Get the indices needed to sort particle vectors
	std::iota(indices.begin(), indices.end(), 0);
	
	std::stable_sort(indices.begin(), indices.end(),
	 [&](int a, int b) {
	 	return cell_idx[a] < cell_idx[b];
	 });
}

void ParticleList::get_cell_start() {
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
std::vector<T> ParticleList::reorder(const std::vector<T>& input, const std::vector<int>& indices) {
	std::vector<T> output(indices.size());
	for (size_t i = 0; i < indices.size(); ++i) {
		output[i] = input[indices[i]];
	}
	return output;
}


void ParticleList::compart() {
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
	auto new_type = reorder(type, indices);

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
	type = std::move(new_type);
}

void ParticleList::compute_rho() {
	std::fill(rho.begin(), rho.end(), 0);

	#pragma omp parallel for
	for (int cell=0; cell<params.no_cells; cell++) {

		auto [cx, cy] = scalar_index_to_vec(cell);

		for (auto [dx,dy] : neighbor_offsets) {
			int ncx = cx + dx;
			int ncy = cy + dy;

			double y_offset = 0.0;
			int neighbor;
			if (ncx<0 || ncx>=params.nc[0]) continue;

			if (ncy<0) { // Bottom -> Top
				neighbor = vec_to_scalar_index(ncx, ncy+params.nc[1]);
				y_offset = -params.Ly;
			}

			else if (ncy>=params.nc[1]) { // Top -> Bottom
				neighbor = vec_to_scalar_index(ncx, ncy-params.nc[1]);
				y_offset = params.Ly;
			}

			else { // Base
				neighbor = vec_to_scalar_index(ncx, ncy);
			}

			update_rho(cell, neighbor, y_offset);
		}
	}
}

inline void ParticleList::update_rho(int cell, int neighbor, double y_offset) {
	double h = params.h;
	double rc = params.rc;
	double h2 = h*h;
	double rc2 = rc*rc;

	int start_i = cell_start[cell];
	int end_i = cell_start[cell+1];

	int start_j = cell_start[neighbor];
	int end_j = cell_start[neighbor+1];

	for (int i=start_i; i<end_i; ++i) {
		double xi = pos_x[i];
		double yi = pos_y[i];

		for (int j=start_j; j<end_j; ++j) {
			double xj = pos_x[j];
			double yj = pos_y[j]+y_offset;

			double rel_x = xi - xj;
			double rel_y = yi - yj;
			double r2 = rel_x*rel_x + rel_y*rel_y;
			if (r2 < rc2) {
				double q = std::sqrt(r2/h2);
				rho[i] += mass[j] * kernel(q);
			}
		}
	}
}

void ParticleList::compute_p() {
	#pragma omp parallel for
	for (int i=0; i<no_particles; ++i) {
		pressure[i] = calculate_pressure(rho[i]);
	}
}

void ParticleList::compute_forces() {
	std::fill(acc_x.begin(), acc_x.end(), 0);
	std::fill(acc_y.begin(), acc_y.end(), 0);
	
	#pragma omp parallel for
	for (int cell=0; cell<params.no_cells; cell++) {
		auto [cx, cy] = scalar_index_to_vec(cell);

		for (auto [dx,dy] : neighbor_offsets) {
			int ncx = cx + dx;
			int ncy = cy + dy;

			double y_offset = 0.0;
			int neighbor;
			if (ncx<0 || ncx>=params.nc[0]) continue;

			if (ncy<0) { // Bottom -> Top
				neighbor = vec_to_scalar_index(ncx, ncy+params.nc[1]);
				y_offset = -params.Ly;
			}

			else if (ncy>=params.nc[1]) { // Top -> Bottom
				neighbor = vec_to_scalar_index(ncx, ncy-params.nc[1]);
				y_offset = params.Ly;
			}

			else { // Base
				neighbor = vec_to_scalar_index(ncx, ncy);
			}

			update_forces(cell, neighbor, y_offset);
		}
	}
}

inline void ParticleList::update_forces(int cell, int neighbor, double y_offset) {
	double h = params.h;
	double rc = params.rc;
	double rc2 = rc*rc;

	int start_i = cell_start[cell];
	int end_i = cell_start[cell+1];

	int start_j = cell_start[neighbor];
	int end_j = cell_start[neighbor+1];

	for (int i=start_i; i<end_i; ++i) {
		double xi = pos_x[i];
		double yi = pos_y[i];

		for (int j=start_j; j<end_j; ++j) {
			if (type[i] != ParticleType::mainflow) continue; // Only mainflow particles feel forces from other particles
			if (i == j) continue; // Particles do not interact with themselves

			double xj = pos_x[j];
			double yj = pos_y[j];

			double rel_x = xi-xj;
			double rel_y = yi-yj;
			double rel_velx = 0.0;
			double rel_vely = 0.0;
			double r2 = rel_x*rel_x + rel_y*rel_y;
			
			if (r2 < rc2 && r2 > 1e-12) {
				// --------- Pressure Force -----------
				double r = std::sqrt(r2);
				double q = r/h;

				double common = -1*mass[j]*
				(pressure[i]/(rho[i]*rho[i]) + 
					pressure[j]/(rho[j]*rho[j]))
				*deriv_kernel(q)/(h*r);

				acc_x[i] += common * rel_x;
				acc_y[i] += common * rel_y;
				// -------------------------------------

				// --------- Viscosity force -----------
				if (type[j] == ParticleType::ghost) {
					// Ghost particles get a no-slip artificial velocity
					// Ghost particles get a no-slip artificial velocity
					double d_i = std::sqrt(pow(xi-params.Lx/2, 2)+pow(yi-params.Ly/2, 2)) - params.central_radius;
					double d_j = params.central_radius - ((xi-params.Lx/2)*(xj-params.Lx/2) + (yi-params.Ly/2)*(yj-params.Ly/2))/std::sqrt(pow(xi-params.Lx/2,2)+pow(yi-params.Ly/2,2));//std::sqrt(pow(xj-params.Lx/2, 2)+pow(yj-params.Ly/2, 2)) - params.central_radius;
					double beta = 1 + d_j/d_i;
					if ( beta > params.max_beta ) beta = params.max_beta;
					if ( beta < 1 ) beta = 1;
					rel_velx = beta*vel_x[i];
					rel_vely = beta*vel_y[i];
				}
				else {
					// Normal particles just have relative velocities
					rel_velx = vel_x[i] - vel_x[j];
					rel_vely = vel_y[i] - vel_y[j];
				}

				double mu_i = params.kinematic_viscosity*rho[i];
				double mu_j = params.kinematic_viscosity*rho[j];

				double common2 = mass[j]*(mu_i+mu_j)/(rho[i]*rho[j])*(deriv_kernel(q)/r);

				acc_x[i] += common2 * rel_velx;
				acc_y[i] += common2 * rel_vely;

				if (type[j] == ParticleType::ghost){
					acc_x[j] -= common*rel_x + common2*rel_velx;
					acc_y[j] -= common*rel_y + common2*rel_vely;
				}
				// ------------------------------------
			}
		}
	}
}

void ParticleList::compute_v_half() {
	#pragma omp parallel for
	for (int i=0; i<no_particles; ++i) {
		vel_x[i] += 0.5 * params.dt * acc_x[i];
		vel_y[i] += 0.5 * params.dt * acc_y[i];
	}
}

void ParticleList::compute_pos() {
	#pragma omp parallel for
	for (int i=0; i<no_particles; ++i) {
		if (type[i] != ParticleType::ghost){ //Ghost particles should not move
			pos_x[i] += params.dt * vel_x[i];
			pos_y[i] += params.dt * vel_y[i];
		}
	}
}

void ParticleList::compute_v_full() {
	#pragma omp parallel for
	for (int i=0; i<no_particles; ++i) {
		vel_x[i] += 0.5 * params.dt * acc_x[i];
		vel_y[i] += 0.5 * params.dt * acc_y[i];
	}
}


void ParticleList::periodic() {
	#pragma omp parallel for
	for (int i = 0; i < no_particles; ++i) {
		pos_y[i] = std::fmod(pos_y[i], params.Ly);
		if (pos_y[i] < 0) pos_y[i] += params.Ly;
	}
}

void ParticleList::compute_drag_force() {
	drag_force = 0.0;
	for (int i=0; i<no_particles; ++i) {
		if (type[i] == ParticleType::ghost) {
			drag_force += acc_x[i]*mass[i];
		}
	}
}

void ParticleList::compute_max_a() {
	for (int i=0; i<no_particles; ++i) {
		double a = std::sqrt(pow(acc_x[i], 2) + pow(acc_y[i], 2));
		if (a > max_a) {
			max_a = a;
		}
	}
}

void ParticleList::compute_cd() {
	cd = (2*drag_force/(params.v0*params.v0*params.rho0*2*params.central_radius));
}

void ParticleList::integrate() {
	manage_inflow_outflow();
	compart();
	compute_rho();

	compute_p();

	compute_forces();

    compute_drag_force();
    compute_cd();

	compute_v_half();
	compute_pos();
	periodic();
	compute_v_full();

	compute_max_a();
}

// Explicit template instantiations
template std::vector<double> ParticleList::reorder<double>(const std::vector<double>&, const std::vector<int>&);
template std::vector<int> ParticleList::reorder<int>(const std::vector<int>&, const std::vector<int>&);
template std::vector<ParticleType> ParticleList::reorder<ParticleType>(const std::vector<ParticleType>&, const std::vector<int>&);

