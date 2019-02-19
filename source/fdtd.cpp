//
//  fdtd.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 11.07.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#include "fdtd.hpp"

FDTD::FDTD(Logger* global_log, double accuracy)
{
	/* auto eps = [] (const meep::vec &p) {
		if (p.x() < 2 && p.y() < 3)
			return 12.0;
		return 1.0;
	};

	int argc = 0; char** argv;
	meep::initialize mpi = meep::initialize(argc,argv); // do this even for non-MPI Meep
	double resolution = 20;     // pixels per distance
	meep::grid_volume v = meep::vol2d(5, 10, resolution); // 5x10 2d cell
	meep::structure s(v, eps, meep::pml(1.0));
	meep::fields f(&s);

	f.output_hdf5(meep::Dielectric, v.surroundings());

	double freq = 0.3, fwidth = 0.1;
	meep::gaussian_src_time src(freq, fwidth);
	f.add_point_source(meep::Ey, src, meep::vec(1.1, 2.3));
	
	while (f.time() < f.last_source_time()) {
		f.step();
	}

	f.output_hdf5(meep::Hz, v.surroundings()); */
}

double FDTD::magnetic_rho (double vt, double rho, double phi, double z) const
{
	throw std::logic_error("FDTD::magnetic_rho is not implemented!");
}

double FDTD::electric_rho (double vt, double rho, double phi, double z) const
{
	throw std::logic_error("FDTD::electric_rho is not implemented!");
}

double FDTD::magnetic_phi (double vt, double rho, double phi, double z) const
{
	throw std::logic_error("FDTD::magnetic_phi is not implemented!");
}

double FDTD::electric_phi (double vt, double rho, double phi, double z) const
{
	throw std::logic_error("FDTD::electric_phi is not implemented!");
}

double FDTD::electric_z   (double vt, double rho, double phi, double z) const
{
	throw std::logic_error("FDTD::electric_z is not implemented!");
}

double FDTD::magnetic_z   (double vt, double rho, double phi, double z) const
{
	throw std::logic_error("FDTD::magnetic_z is not implemented!");
}
