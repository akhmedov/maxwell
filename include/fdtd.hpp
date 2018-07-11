//
//  fdtd.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 11.07.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#ifndef fdtd_hpp
#define fdtd_hpp

#include <meep.hpp>

#include "abstract_field.hpp"

// using namespace meep;

struct FDTD : public AbstractField {
	
	FDTD(Logger* global_log = NULL, double accuracy = 1);

	double magnetic_rho (double vt, double rho, double phi, double z) const;
	double electric_rho (double vt, double rho, double phi, double z) const;
	double magnetic_phi (double vt, double rho, double phi, double z) const;
	double electric_phi (double vt, double rho, double phi, double z) const;
	double electric_z   (double vt, double rho, double phi, double z) const;
	double magnetic_z   (double vt, double rho, double phi, double z) const;
	// double electric_x   (double vt, double rho, double phi, double z) const;
	// double magnetic_x   (double vt, double rho, double phi, double z) const;
	// double electric_y   (double vt, double rho, double phi, double z) const; 
	// double magnetic_y   (double vt, double rho, double phi, double z) const;
};

#endif /* fdtd_hpp */
