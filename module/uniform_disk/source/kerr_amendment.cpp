//
//  kerr_amendment.cpp
//  uniform_disk.module.maxwell
//
//  Created by Rolan Akhmedov on 21.10.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "kerr_amendment.hpp"

#define RHO_EPS 1e-3
#define NU_MAX 10

KerrAmendment::KerrAmendment (double R, double A0, double eps_r, double mu_r, double chi3, Logger* log)
: TransientResponse(R,A0,eps_r,mu_r,log), kerr(chi3) {}

// ===============================================================================================

double KerrAmendment::observed_from (const Point::Cylindrical&) const
{
	throw std::logic_error("KerrAmendment::observed_from is not implemented!");
}

double KerrAmendment::observed_to (const Point::Cylindrical&) const
{
	throw std::logic_error("KerrAmendment::observed_to is not implemented!");
}

// ===============================================================================================

double KerrAmendment::electric_x (const Point::SpaceTime<Point::Cylindrical>& event) const
{
	throw std::logic_error("KerrAmendment::electric_x is not implemented");
}

double KerrAmendment::electric_rho (const Point::SpaceTime<Point::Cylindrical>& event) const
{
	double coeff = std::pow((MU0 * MU) / (EPS0 * EPS), 2);
	coeff *= EPS0 * kerr * std::pow(A0,3) / std::pow(2,7);
	coeff *= std::sqrt(eps_r/(eps_r+kerr));

	auto extantion = [event, this] (double nu) {
		Point::ModalSpaceTime<Point::Cylindrical> mode1(1, nu, event);
		Point::ModalSpaceTime<Point::Cylindrical> mode3(3, nu, event);
		double term1 = this->modal_vmh(mode1) * (jn(0, event.rho() * nu) + jn(2, event.rho() * nu));
		double term3 = this->modal_vmh(mode3) * (jn(2, event.rho() * nu) + jn(4, event.rho() * nu));
		return nu * (std::cos(event.phi()*1) * term1 + std::cos(event.phi()*3) * term3) / 2;
	};

	SimpsonRunge integrate = SimpsonRunge(10 /* init_terms */, 1 /* % */);
	try {
		return coeff * integrate.value(0, NU_MAX, extantion);
	} catch (double val) {
		// TODO: add logging of event
		return 0;
	}
}

double KerrAmendment::electric_phi (const Point::SpaceTime<Point::Cylindrical>&) const
{
	throw std::logic_error("KerrAmendment::electric_phi is not implemented");
}

double KerrAmendment::electric_z (const Point::SpaceTime<Point::Cylindrical>&) const
{
	return 0;
}

double KerrAmendment::magnetic_rho (const Point::SpaceTime<Point::Cylindrical>&) const
{
	throw std::logic_error("KerrAmendment::magnetic_rho is not implemented");
}

double KerrAmendment::magnetic_phi (const Point::SpaceTime<Point::Cylindrical>&) const
{
	throw std::logic_error("KerrAmendment::magnetic_phi is not implemented");
}

double KerrAmendment::magnetic_z (const Point::SpaceTime<Point::Cylindrical>&) const
{
	throw std::logic_error("KerrAmendment::magnetic_z is not implemented");
}

/* */

double KerrAmendment::modal_vmh (const Point::ModalSpaceTime<Point::Cylindrical>& event) const // Vmh * ct
{
	double ct_z = event.ct() * event.ct() - event.z() * event.z();
	if (ct_z < 0) return 0;

	auto unterint_vmh = [event, this] (double z) {
		double j0arg = event.nu() * std::sqrt(2 * event.ct() * (event.z() + z) - 2 * event.z() * z);
		auto tmp = event; tmp.ct() = event.ct()-event.z()+z; tmp.z() = z; // TODO: test copy constructor
		double term1 = j0(j0arg) * (event.ct()-event.z()+z) * this->modal_jm(tmp);

		SimpsonRunge timeIntegr = SimpsonRunge(10 /* init_terms */, 10 /* % */);
		auto term2 = [event, z, this] (double ct) {
			double besselArg = event.nu() * std::sqrt((event.ct() - ct)*(event.ct() - ct) - (event.z() - z)*(event.z() - z));
			auto tmp = event; tmp.ct() = ct; tmp.z() = z; // TODO: test copy constructor
			return (event.ct() - ct) * (j0(besselArg) + jn(2,besselArg)) * ct * this->modal_jm(tmp);
		};

		try{
			return term1 + event.nu() * event.nu() * timeIntegr.value(0, event.ct()-event.z()+z, term2) / 2;
		} catch (double val) {
			// TODO: add logging of event
			return term1;
		}
	};

	SimpsonRunge distIntegr = SimpsonRunge(10 /* init_terms */, 10 /* % */);
	try {
		return distIntegr.value(0, event.z(), unterint_vmh); // TODO: maybe 5R must be used
	} catch (double val) {
		// TODO: add logging of event
		return 0;
	}
}

double KerrAmendment::modal_jm (const Point::ModalSpaceTime<Point::Cylindrical>& event) const // kerr_jm/ct from thesis.pdf
{
	double ct_z = event.ct() * event.ct() - event.z() * event.z();
	double R = this->R;
	if (ct_z < 0) return 0;

	auto underint_j1 = [event, ct_z, R] (double rho) {
		double a = KerrAmendment::alpha(ct_z, rho, R);
		double b = KerrAmendment::beta(ct_z, rho, R);
		double g = KerrAmendment::gamma(ct_z, rho, R);
		double l = KerrAmendment::lambda(ct_z, rho, R);
		double term1 =    j0(event.nu() * rho) * (3*a + b + 3*g + l); 
		double term2 = jn(2, event.nu() * rho) * (3*a + b - 3*g - l);
		return rho * (term1 + term2);
	};

	auto underint_j3 = [event, ct_z, R] (double rho) {
		double a = KerrAmendment::alpha(ct_z, rho, R);
		double b = KerrAmendment::beta(ct_z, rho, R);
		double g = KerrAmendment::gamma(ct_z, rho, R);
		double l = KerrAmendment::lambda(ct_z, rho, R);
		double term1 = jn(2, event.nu() * rho) * (a - b - g + l); 
		double term2 = jn(4, event.nu() * rho) * (a - b + g - l);
		return rho * (term1 + term2);
	};

	double from = std::abs(std::sqrt(ct_z) - R) + RHO_EPS;
	double to = std::sqrt(ct_z) + R - RHO_EPS;
	SimpsonRunge integrate = SimpsonRunge(10 /* init_terms */, 1 /* % */);

	switch (static_cast<int>(event.m())) {
		case 1: case -1: 
			try {
				return integrate.value(from, to, underint_j1);
			} catch (double val) {
				// TODO: add logging of event
				return 0;
			}
		case 3: case -3:
			try {
				return integrate.value(from, to, underint_j3);
			} catch (double val) {
				// TODO: add logging of event
				return 0;
			}
	}

	return 0;
}

/* normed intensity components of linear wave */

double KerrAmendment::alpha (double ct_z, double rho, double R) // alpha / ct
{
	double sqrt_vt_z = std::sqrt(ct_z);
	double i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(ct_z, rho, R);
	return 3 * i1 * i1 * i1_perp;
}

double KerrAmendment::beta (double ct_z, double rho, double R) // beta / ct
{
	double sqrt_vt_z = std::sqrt(ct_z);
	double i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = TransientResponse::int_bessel_001(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(ct_z, rho, R);
	double i2_perp = KerrAmendment::int_bessel_001_perp(ct_z, rho, R);
	return (i2-i1) * ( i1_perp * (i2-i1) + 2 * i1 * (i2_perp - i1_perp) );
}

double KerrAmendment::gamma (double ct_z, double rho, double R)  // gamma / ct
{
	double sqrt_vt_z = std::sqrt(ct_z);
	double i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = TransientResponse::int_bessel_001(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(ct_z, rho, R);
	double i2_perp = KerrAmendment::int_bessel_001_perp(ct_z, rho, R);
	return 3 * (i2-i1) * (i2-i1) * (i2_perp - i1_perp);
}

double KerrAmendment::lambda (double ct_z, double rho, double R)  // lambda / ct
{
	double sqrt_vt_z = std::sqrt(ct_z);
	double i1 = TransientResponse::int_bessel_011(sqrt_vt_z, rho, R);
	double i2 = TransientResponse::int_bessel_001(sqrt_vt_z, rho, R);
	double i1_perp = KerrAmendment::int_bessel_011_perp(ct_z, rho, R);
	double i2_perp = KerrAmendment::int_bessel_001_perp(ct_z, rho, R);
	return i1 * i1 * (i2_perp - i1_perp) + 2 * i1 * i1_perp * (i2-i1);
}

/* partial derivatives of I1 and I2 by time */

double KerrAmendment::int_bessel_011_perp (double ct_z, double rho, double R) // i1_perp / ct
{
	double rho2 = rho * rho;
	double R2 = R * R;

	if (R < std::abs(rho - std::sqrt(ct_z))) return 0.0;
	if (R > rho + std::sqrt(ct_z)) return 0.0;

	double first = (rho2 - R2) * (rho2 - R2) / ct_z;
	first /= std::sqrt( (rho+R)*(rho+R) - ct_z );
	first /= std::sqrt( ct_z - (rho-R)*(rho-R) );

	double second = 2 * (rho2 + R2) - ct_z;
	second /= std::sqrt(4*rho2*R2 - (ct_z - rho2 - R2)*(ct_z - rho2 - R2));

	return (first - second) / (2*M_PI*rho2);
}

double KerrAmendment::int_bessel_001_perp (double ct_z, double rho, double R) // i1_perp / ct
{
	double rho2 = rho * rho;
	double R2 = R * R;

	if (R < std::abs(rho - std::sqrt(ct_z))) return 0.0;
	if (R > rho + std::sqrt(ct_z)) return 0.0;

	double res = ct_z - rho2 + R2;
	res /= std::sqrt(4*rho2*ct_z - (ct_z+rho2-R2) * (ct_z+rho2-R2) );
	return - res / ct_z / M_PI;
}

namespace { extern "C" {

	void tr_module (ModuleManager* core, Logger* log, double R, double A0, double mu, double eps)
	{
		ModuleEntity tr;
		tr.field_cyl_arg = new TransientResponse(R,A0,mu,eps,log);
		const char* name = "UniformDisk.TrancientResponse";
		core->load_module(name,tr);
	}

	void meandr_module (ModuleManager* core, Logger* log, double R, double A0, double tau0, double mu, double eps)
	{
		ModuleEntity squared;
		squared.field_cyl_arg = new SquaredPulse(R,A0,eps,mu,tau0,log);
		const char* name = "UniformDisk.MeandrMonocycle";
		core->load_module(name,squared);
	}

	void kerr_module (ModuleManager* core, Logger* log, double R, double A0, double /*tau0*/, double mu, double eps)
	{
		ModuleEntity kerr;
		kerr.field_cyl_arg = new KerrAmendment(R, A0, eps, mu, 1, log);
		const char* name = "UniformDisk.NonlinearKerrAmendment";
		core->load_module(name,kerr);
	}

	void load_module (ModuleManager* core, Logger* global, double R, double A0, double tau0, double mu, double eps)
	{
		tr_module (core, global, R, A0, mu, eps);
		meandr_module (core, global, R, A0, tau0, mu, eps);
		kerr_module (core, global, R, A0, tau0, mu, eps);
	}

} }
