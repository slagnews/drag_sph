#include <iostream>
#include <vector>
#include <array>
#include <cmath>

struct Vec2 {
	double x, y;

	Vec2 operator*(double coef) const {
		return {x*coef, y*coef};
	}

	Vec2 operator*(const Vec2& other) const {
		return {x*other.x, y*other.y};
	}

	Vec2 operator+(double coef) const {
		return {x+coef, y+coef};
	}

	Vec2 operator+(const Vec2& other) const {
		return {x+other.x, y+other.y};
	}

	Vec2 operator-(double coef) const {
		return {x-coef, y-coef};
	}

	Vec2 operator-(const Vec2& other) const {
		return {x-other.x, y-other.y};
	}

	Vec2& operator=(const Vec2& other) {
		x = other.x;
		y = other.y;
		return *this;
	}

	Vec2& operator+=(double coef) {
		x += coef;
		y += coef;
		return *this;
	}

	Vec2& operator+=(const Vec2& other) {
		x += other.x;
		y += other.y;
		return *this;
	}

	Vec2& operator-=(double coef) {
		x -= coef;
		y -= coef;
		return *this;
	}

	Vec2& operator-=(const Vec2& other) {
		x -= other.x;
		y -= other.y;
		return *this;
	}

	double norm() const {
		return sqrt(pow(x, 2) + pow(y, 2));
	}

	double dist(const Vec2& other) const {
		return sqrt(pow(x-other.x, 2) + pow(y-other.y, 2));
	}
};


enum class ParticleType {
	flow, in, out, ghost
};

struct Particle {
	Vec2 pos, vel, acc;
	double mass, rho, pressure;
	ParticleType type;
};

class SPHSimulation {
	public:
		double c_s;
		double gamma_index;
		double rho0;
		std::vector<Particle> particles;

	void compute_pressure() {

	}
};


// Declare simulation parameters
double Lx, Ly, h, v0, dt, c_s, gamma_index, rc rho0;
int no_particles, nCells, no_steps;
std::array<int, 2> nc;

// Vector of the particle positions
std::vector<std::array<double, 2>> positions;

// Linked list data
std::vector<int> head;
std::vector<int> lscl;

const int EMPTY = -1;

// Helper functions
std::array<int, 2> get_cell_index(const std::array<double, 2>& pos) {
	return { int(pos[0] / rc), int(pos[1] / rc) };
}

int vec_to_scalar_index(int cx, int cy) {
	return cx * nc[1] + cy;
}

// Constructing the linked list for departmentalization
void construct_linked_list() {
	head.assign(nCells, EMPTY);
	lscl.assign(no_particles, EMPTY);

	for (int i = 0; i < no_particles; ++i) {
		auto cell = get_cell_index(positions[i]);
		int cx = cell[0];
		int cy = cell[1];
		int c = vec_to_scalar_index(cx, cy);

		lscl[i] = head[c];
		head[c] = i;
    }
}

double kernel(Vec2 r_vec, double h) {
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

double deriv_kernel(Vec2 r_vec, double h) {
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
		h = atof(argv[5]);
		v0 = atof(argv[6]);
		dt = atof(argv[7]);
		c_s = atof(argv[8]);
		gamma_index = atof(argv[9]);

		rc = 2*h;
		nc = {int(Lx / rc), int(Ly / rc)};
		nCells = nc[0] * nc[1];
	}

	positions = {{0.5, 0.5}, {0.6, 0.6}, {2.3, 1.1}, {0.4, 0.7}};

	construct_linked_list();

	SPHSimulation sim;
	Particle p1, p2;

	p1.pos = {1, 1};
	p1.vel = {1, 0};
	p1.acc = {0, 0};
	p1.mass = 1.0;
	p1.rho = 0;
	p1.pressure = 0;
	p1.type = ParticleType::flow;
	sim.particles.push_back(p1);

	p2.pos = {2, 1};
	p2.vel = {-1, 0};
	p2.acc = {0, 0};
	p2.mass = 1.0;
	p2.rho = 0;
	p2.pressure = 0;
	p2.type = ParticleType::flow;
	sim.particles.push_back(p2);

	return 0;
}