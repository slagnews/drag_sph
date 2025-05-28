#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <numeric>

struct ParticleList {
	std::vector<double> pos_x;
	std::vector<double> pos_y;
	std::vector<double> mass;
	std::vector<double> rho;
	std::vector<double> pressure;
	std::vector<int> cell_idx;
};


// Declare simulation parameters
double Lx, Ly, h, v0, dt, c_s, gamma_index, rc, rho0;
int no_particles, no_cells, no_steps;
std::array<int, 2> nc;

// Helper functions
std::array<int, 2> get_cell_index(double x, double y) {
	return { int(x / rc), int(y / rc) };
}

int vec_to_scalar_index(int cx, int cy) {
	return cx * nc[1] + cy;
}

// Assign the corresponding cell to each particle
void assign_cell(ParticleList& pl) {

	for (int i = 0; i < no_particles; ++i) {
		auto cell = get_cell_index(pl.pos_x[i], pl.pos_y[i]);
		int cx = cell[0];
		int cy = cell[1];
		int c = vec_to_scalar_index(cx, cy);
		pl.cell_idx[i] = c;
    }
}

std::vector<int> sorter(int no_particles, const ParticleList& pl) {
	std::vector<int> indices(no_particles);
	std::iota(indices.begin(), indices.end(), 0);
	
	std::stable_sort(indices.begin(), indices.end(),
	 [&](int a, int b) {
	 	return pl.cell_idx[a] < pl.cell_idx[b];
	 });
	 return indices;
}

template<typename T>
std::vector<T> reorder(const std::vector<T>& input, const std::vector<int>& indices) {
	std::vector<T> output(indices.size());
	for (size_t i = 0; i < indices.size(); ++i) {
		output[i] = input[indices[i]];
	}
	return output;
}


/*
double kernel(std::array<double, 2> r_vec, double h) {
	double r = r_vec.norm();
	double q = r/h;
	double sigma = (15/(14*M_PI*pow(h, 2)));

	if (q > 2) {
		return 0;

	} else if ((q > 1) and (q <= 2)) {
		return sigma*pow(2-q, 3);

	} else if ((q >= 0) and (q <= 1)) {
		return sigma*(pow(2-q, 3)-4*pow(1-q, 3));

	} else {
		std::cout << "q smaller than 0" << std::endl;
		return 0;
	}
}

double deriv_kernel(std::array<double, 2> r_vec, double h) {
	double r = r_vec.norm();
	double q = r/h;
	double sigma = (15/(14*M_PI*pow(h, 3)));

	if (q > 2) {
		return 0;

	} else if ((q > 1) and (q <= 2)) {
		return sigma*(-3*pow(2-q, 2));

	} else if ((q >= 0) and (q <= 1)) {
		return sigma*(-3*pow(2-q, 2)+12*pow(1-q, 2));

	} else {
		std::cout << "q smaller than 0" << std::endl;
		return 0;
	}
}
*/

double cole_pressure(double rho, double c_s, double rho0, double gamma_index) {
	double B = pow(c_s, 2)*rho0/gamma_index;
	return B*(pow(rho/rho0, gamma_index)-1);
}


int main(int argc, char *argv[]) {
	if (argc >= 2) {
		no_particles = atoi(argv[1]);
		no_steps = atoi(argv[2]);
		Lx = atof(argv[3]);
		Ly = atof(argv[4]);
		rho0 = atof(argv[5]);
		h = atof(argv[6]);
		v0 = atof(argv[7]);
		dt = atof(argv[8]);
		c_s = atof(argv[9]);
		gamma_index = atof(argv[10]);

		rc = 2*h;
		nc = {int(Lx / rc), int(Ly / rc)};
		no_cells = nc[0] * nc[1];
	}
	
	double particle_mass = Lx*Ly*rho0/2;

	ParticleList pl = {{1.2, 8, 5.5}, {4, 4, 6.2}, {particle_mass, particle_mass, particle_mass}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	
	std::cout << pl.pos_x[0] << std::endl;
	std::cout << pl.pos_x[1] << std::endl;
	std::cout << pl.pos_x[2] << std::endl << std::endl;
	
	assign_cell(pl);
	auto indices = sorter(no_particles, pl);
	pl.pos_x = reorder(pl.pos_x, indices);
	
	std::cout << pl.pos_x[0] << std::endl;
	std::cout << pl.pos_x[1] << std::endl;
	std::cout << pl.pos_x[2] << std::endl;
	
	

	return 0;
}
