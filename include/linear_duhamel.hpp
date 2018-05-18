//
//  linear_duhamel.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 16.05.18.
//  Copyright © 2018 Rolan Akhmedov. All rights reserved.
//

#ifndef UNUSED
#define UNUSED(expr) do { (void)(expr); } while (0)
#endif

#ifndef linear_duhamel_hpp
#define linear_duhamel_hpp

#define INIT_NODES 1e5
#define MAX_NODES  1e7

#include "logger.hpp"
#include "config.hpp"
#include "integral.hpp"
#include "phys_math.hpp"
#include "linear_field.hpp"
#include "linear_medium.hpp"
#include "linear_source.hpp"

struct FreeTimeCurrent : public LinearCurrent {
	FreeTimeCurrent (LinearCurrent* on_source, double duration);
	void set_time_depth (const std::function<double(double)>& func);
	double time_shape (double vt) const;
	double rho (double rho, double phi, double z) const;
	double phi (double rho, double phi, double z) const;
	double z (double rho, double phi, double z) const;
private:
	// double duration;
	LinearCurrent* base;
	std::function<double(double)> time_fnc;
};

struct LinearDuramel : public LinearField {
	LinearDuramel (LinearCurrent* source, LinearMedium* medium, LinearField* on, Logger* global_log = NULL);
	void   set_accuracy (double persent);
	double electric_rho (double vt, double rho, double phi, double z) const;
	double electric_phi (double vt, double rho, double phi, double z) const;
	double electric_z   (double vt, double rho, double phi, double z) const;
	double magnetic_rho (double vt, double rho, double phi, double z) const;
	double magnetic_phi (double vt, double rho, double phi, double z) const;
	double magnetic_z   (double vt, double rho, double phi, double z) const;
private:
	Logger* global_log;
	LinearField* on_field;
	double accuracy = 1;
};

#endif /* linear_duhamel_hpp */