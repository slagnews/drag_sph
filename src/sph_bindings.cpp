#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "sph.cpp"  // Or include a header if you split into .h/.cpp

namespace py = pybind11;

PYBIND11_MODULE(sph_cpp, m) {
    py::class_<SimulationParams>(m, "SimulationParams")
        .def(py::init<double, int, double, double, double, double, double, double, double, double>())
        .def_readonly("Lx", &SimulationParams::Lx)
        .def_readonly("Ly", &SimulationParams::Ly)
        .def_readonly("Nx", &SimulationParams::Nx)
        .def_readonly("Ny", &SimulationParams::Ny)
        .def_readonly("init_particles", &SimulationParams::init_particles);

    py::class_<ParticleList>(m, "ParticleList")
        .def(py::init<const SimulationParams&>())
        .def("init_grid", &ParticleList::init_grid)
        .def("init_random", &ParticleList::init_random)
        .def("init_mass", &ParticleList::init_mass)
        .def("init_vel", &ParticleList::init_vel)
        .def("integrate", &ParticleList::integrate);
        //.def("set_initial_positions", &ParticleList::set_initial_positions)
        //.def("set_initial_velocities", &ParticleList::set_initial_velocities);
        /**
        .def("get_positions", [](const ParticleList &pl) {
            py::ssize_t N = static_cast<py::ssize_t>(pl.pos_x.size());
            auto result = py::array_t<double>({N, 2});
            auto r = result.mutable_unchecked<2>();
            for (py::ssize_t i = 0; i < N; ++i) {
                r(i, 0) = pl.pos_x[i];
                r(i, 1) = pl.pos_y[i];
            }
            return result;
        });
        */


}