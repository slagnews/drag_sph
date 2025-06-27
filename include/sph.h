#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <random>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <chrono>

constexpr int neighbor_count = 9;
extern const int neighbor_offsets[neighbor_count][2];

enum ParticleType {
	inflow,
	outflow,
	mainflow,
	ghost
};

struct SimulationParams {
	double res;
	int no_steps;
	double Lx, Ly, rho0, kappa, v0, dt, c_s, beta, central_radius, inflow_factor, outflow_factor, kinematic_viscosity, max_beta;
	int no_cells, Nx, Ny, init_particles;
	std::array<int, 2> nc;
	double h, rc, P0, sigma;

	SimulationParams(double res, int no_steps, double Lx, double Ly, double rho0, double kappa, double v0, double dt,
		double c_s, double beta, double central_radius, double inflow_factor,
		double outflow_factor, double kinematic_viscosity, double max_beta);
};

struct ParticleList {
	const SimulationParams& params;
	int no_particles;

	std::vector<double> pos_x, pos_y, vel_x, vel_y, acc_x, acc_y, mass, rho, pressure;
	std::vector<int> cell_idx, indices;
	std::vector<ParticleType> type;
	std::vector<int> cell_counts, cell_start;
	double drag_force, cd, max_a;

	ParticleList(const SimulationParams& p);

	// Helper Functions
	std::array<int, 2> get_cell_index(double x, double y) const;
	int vec_to_scalar_index(int cx, int cy) const;
	std::array<int, 2> scalar_index_to_vec(int cell) const;
	inline double kernel(double q) const;
	inline double deriv_kernel(double q) const;
	inline double calculate_pressure(double rho) const;

	// Initialization
	void init_grid();
	void init_mass();
	void init_vel();
	void init_type();

	// Cell and sorting
	void assign_cell();
	void manage_inflow_outflow();
	void get_sorter();
	void get_cell_start();
	void compart();
	template<typename T>
	std::vector<T> reorder(const std::vector<T>& input, const std::vector<int>& indices);

	// Simulation Functions
	void compute_rho();
	void update_rho(int cell, int neighbor, double y_offset);
	void compute_p();
	void compute_forces();
	void update_forces(int cell, int neighbor, double y_offset);
	void compute_v_half();
	void compute_pos();
	void compute_v_full();
	void periodic();

	// Measurables
	void compute_drag_force();
	void compute_max_a();
	void compute_cd();

	// Main Simulation loop
	void integrate();

};

#endif
