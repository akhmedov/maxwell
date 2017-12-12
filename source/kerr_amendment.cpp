//
//  kerr_amendment.cpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 21.10.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#include "kerr_amendment.hpp"

KerrMedium::KerrMedium (double realative_mu, double realative_eps, double chi3_electric, double conduct) 
: Homogeneous(realative_mu, realative_eps) {
	if (chi3_electric > 0)
		this->chi3 = chi3_electric;
	else
		std::invalid_argument("Kerr coefficient must be positive.");
	
	if (conduct > 0) 
		this->sigma = conduct;
	else
		std::invalid_argument("Kerr coefficient must be positive.");
}

double KerrMedium::conductivity (double ct, double z) const
{
	UNUSED(ct); UNUSED(z);
	return this->sigma;
}

double KerrMedium::relative_permittivity (double ct, double z, std::size_t term) const
{
	UNUSED(ct); UNUSED(z);
	switch(term) {
		// case 0: return 0;
		case 1: return Homogeneous::relative_permittivity(ct,z);
		case 3: return this->chi3;
		default: return 0;
	}
}

double KerrMedium::relative_permeability (double ct, double z, std::size_t term) const
{
	UNUSED(ct); UNUSED(z);
	switch(term) {
		// case 0: return 0;
		case 1: return Homogeneous::relative_permeability(ct,z);
		case 3: return 0;
		default: return 0;
	}
}

/* */

KerrAmendment::KerrAmendment (MissileField* field, KerrMedium* medium)
: MissileField(*field), NonlinearField(field, medium) { }

double KerrAmendment::electric_rho (double vt, double rho, double phi, double z) const
{

	/* Monte-Carlo Inegration 
	auto field = [this, rho, phi, z] (double vt) {

		auto mode = [this, vt, rho, phi, z] (int m, double varrho, double z_perp, double vt_perp, double nu) {
			double mu_r = nl_medium->relative_permeability(vt,z);
			double sqrt_eps0 = std::sqrt(NonlinearMedium::EPS0);
			double G = this->riemann(nu, vt_perp - vt, z_perp - z);
			double res = std::cos(m * phi);
			res *= m * G * jn(m, nu * rho);
			res /= rho * std::sqrt(nu);
			res *= varrho * this->im_modal_source(m,nu,vt,varrho,z);
			return res * mu_r / sqrt_eps0;
		};

		auto modes_sum = [this, vt, rho, phi, z, mode] (double varrho, double z_perp, double vt_perp, double nu) {
			double sum = 0;
			sum += mode(-1, varrho, z_perp, vt_perp, nu);
			sum += mode(1, varrho, z_perp, vt_perp, nu);
			sum += mode(-3, varrho, z_perp, vt_perp, nu);
			sum += mode(3, varrho, z_perp, vt_perp, nu);
			return sum;
		};

		std::vector<std::pair<double,double>> limits;
		limits.push_back(std::make_pair(0,10e3));
		limits.push_back(std::make_pair(0,10e3));
		limits.push_back(std::make_pair(0,10e3));
		limits.push_back(std::make_pair(0,10e3));
		MonteCarlo integral = MonteCarlo(10e8, limits);
		return integral.value(modes_sum);
	};

	return Math::derivative(field,vt); */

	/* Simpson Inegration
	double R = this->R;
	Simpson I2 = Simpson(100);
	Simpson I3 = Simpson(10e3);
	Simpson I4 = Simpson(10e4);

	auto field = [&I2, &I3, &I4, &R, &mode] (double vt) {
		auto int_zperp = [&I2, &I3, &I4, &R, &mode] (double nu) {
			auto int_vtperp = [&I2, &I4, &R, &mode, &nu] (double vt_perp) {  
				auto int_varrho = [&I2, &I4, &mode, &nu, &vt_perp] (double z_perp) {
					std::cout << nu << ' ' << vt_perp << ' ' << z_perp << std::endl;
					auto modes_sum = [&mode, &z_perp, &vt_perp, &nu] (double varrho) {
						double sum = 0;
						sum += mode(-1, varrho, z_perp, vt_perp, nu);
						sum += mode( 1, varrho, z_perp, vt_perp, nu);
						sum += mode(-3, varrho, z_perp, vt_perp, nu);
						sum += mode( 3, varrho, z_perp, vt_perp, nu);
						return sum;
					};
					return I4.value(0, 10e3, modes_sum);
				};
				return I2.value(vt_perp, vt_perp + 2*R, int_varrho);
			};
			return I3.value(0, 10e2, int_vtperp);
		};
		return I4.value(0, 10e3, int_zperp);
	}; */


	auto field = [this, rho, phi, z] (double vt) {

		auto mode = [this, vt, rho, phi, z] (int m, double rho_perp, double z_perp, double vt_perp, double nu) {
			double mu_r = nl_medium->relative_permeability(vt,z);
			double sqrt_eps0 = std::sqrt(NonlinearMedium::EPS0);
			double G = this->riemann(nu, vt_perp - vt, z_perp - z);
			double res = std::cos(m * phi);
			res *= m * G * jn(m, nu * rho);
			res /= rho * std::sqrt(nu);
			res *= rho_perp * this->im_modal_source(m,nu,vt,rho_perp,z);
			return res * mu_r / sqrt_eps0;
		};

		auto modes_sum = [mode] (double rho_perp, double z_perp, double vt_perp, double nu) {
			double sum = 0;
			sum += mode(-1, rho_perp, z_perp, vt_perp, nu);
			sum += mode( 1, rho_perp, z_perp, vt_perp, nu);
			sum += mode(-3, rho_perp, z_perp, vt_perp, nu);
			sum += mode( 3, rho_perp, z_perp, vt_perp, nu);
			return sum;
		};

		std::vector< std::tuple<double,std::size_t,double> > limits;
		limits.push_back(std::make_tuple(0, 100, 10e3)); // rho_perp
		limits.push_back(std::make_tuple(0, 100, 10e3)); // z_perp
		limits.push_back(std::make_tuple(0, 100, 10e3)); // vt_perp
		limits.push_back(std::make_tuple(0, 100, 10e3)); // nu_perp
		SimpsonMultiDim integral = SimpsonMultiDim(limits);
		return integral.value(modes_sum);
	};


	return Math::derivative(field,vt);
}

double KerrAmendment::electric_phi (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::electric_phi is not implemented");
}

double KerrAmendment::electric_z (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	return 0;
}

double KerrAmendment::magnetic_rho (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_rho is not implemented");
}

double KerrAmendment::magnetic_phi (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_phi is not implemented");
}

double KerrAmendment::magnetic_z (double ct, double rho, double phi, double z) const
{
	UNUSED(ct); UNUSED(phi); UNUSED(rho); UNUSED(z);
	throw std::logic_error("KerrAmendment::magnetic_z is not implemented");
}

/* */

double KerrAmendment::riemann (double nu, double vt_diff, double z_diff) const
{ 
	double distance = vt_diff * vt_diff - z_diff * z_diff;
	if (distance > 0) distance = std::sqrt(distance);
	else std::invalid_argument("Interval is not legal!");

	double eps_r = nl_medium->relative_permittivity(vt_diff,z_diff);
	double mu_r = nl_medium->relative_permeability(vt_diff,z_diff);
	double velocity = NonlinearMedium::C / (eps_r * mu_r);

	return j0(nu * distance) * velocity / 2;
}

double KerrAmendment::im_modal_source (int m, double nu, double ct, double varrho, double z) const
{
	double terms = 0;

	switch (m) {
		case -1: { 
			terms -= 3 * KerrAmendment::N1(-1,nu,ct,varrho,z);
			terms -= 1 * KerrAmendment::N2(-1,nu,ct,varrho,z);
			terms -= 1 * KerrAmendment::N3(-1,nu,ct,varrho,z);
			terms += 3 * KerrAmendment::N4(-1,nu,ct,varrho,z);
			terms += 1 * KerrAmendment::N5(-1,nu,ct,varrho,z);
			terms += 1 * KerrAmendment::N6(-1,nu,ct,varrho,z);
			break;
		}
		case 1: {
			terms -= 3 * KerrAmendment::N1(1,nu,ct,varrho,z);
			terms -= 1 * KerrAmendment::N2(1,nu,ct,varrho,z);
			terms -= 1 * KerrAmendment::N3(1,nu,ct,varrho,z);
			terms -= 3 * KerrAmendment::N4(1,nu,ct,varrho,z);
			terms -= 1 * KerrAmendment::N5(1,nu,ct,varrho,z);
			terms -= 1 * KerrAmendment::N6(1,nu,ct,varrho,z);
			break; 
		}
		case -3: {
			terms -= 1 * KerrAmendment::N1(-3,nu,ct,varrho,z);
			terms += 1 * KerrAmendment::N2(-3,nu,ct,varrho,z);
			terms -= 1 * KerrAmendment::N4(-3,nu,ct,varrho,z);
			terms += 1 * KerrAmendment::N5(-3,nu,ct,varrho,z);
			break; 
		}
		case 3: { 
			terms += 1 * KerrAmendment::N1(3,nu,ct,varrho,z);
			terms += 1 * KerrAmendment::N2(3,nu,ct,varrho,z);
			terms -= 1 * KerrAmendment::N4(3,nu,ct,varrho,z);
			terms += 1 * KerrAmendment::N5(3,nu,ct,varrho,z);
			break; 
		}

		default: terms = 0;
	}
	
	return terms;
}

double KerrAmendment::N1 (int m, double nu, double ct, double varrho, double z) const
{
	// TODO: check for varrho dependency

	double A3 = this->A0 * this->A0 * this->A0;
	double chi = nl_medium->relative_permittivity(ct, z, 3);

	double eps_r = nl_medium->relative_permittivity(ct,z);
	double mu_r = nl_medium->relative_permeability(ct,z);
	double velocity = NonlinearMedium::C / (eps_r * mu_r);
	double em_relation = 1;
	em_relation *= std::sqrt(NonlinearMedium::MU0 * mu_r);
	em_relation /= std::sqrt(NonlinearMedium::EPS0 * eps_r);

	mpf_class vt_z = mpf_class(ct * ct);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );
	double vt = ct / std::sqrt(eps_r * mu_r);

	double i1 = this->int_bessel_011(std::sqrt(vt_z.get_d()), varrho, this->R);
	double i1_perp = velocity * this->vint_bessel_011_perp(vt, z, varrho, this->R);

	double res = 3 * m * A3 * chi * std::sqrt(NonlinearMedium::MU0);
	res /= 64 * std::sqrt(nu);
	res *= em_relation * em_relation * em_relation;
	res *= jn(m,varrho*nu) * i1 * i1 * i1_perp;
	return res;
}

double KerrAmendment::N2 (int m, double nu, double ct, double varrho, double z) const
{
	double A3 = this->A0 * this->A0 * this->A0;
	double chi = nl_medium->relative_permittivity(ct, z, 3);

	double eps_r = nl_medium->relative_permittivity(ct,z);
	double mu_r = nl_medium->relative_permeability(ct,z);
	double velocity = NonlinearMedium::C / (eps_r * mu_r);
	double em_relation = 1;
	em_relation *= std::sqrt(NonlinearMedium::MU0 * mu_r);
	em_relation /= std::sqrt(NonlinearMedium::EPS0 * eps_r);

	mpf_class vt_z = mpf_class(ct * ct);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );

	double vt = ct / std::sqrt(eps_r * mu_r);

	double i1 = this->int_bessel_011(std::sqrt(vt_z.get_d()), varrho, this->R);
	double i2 = this->int_bessel_001(std::sqrt(vt_z.get_d()), varrho, this->R);
	double i1_perp = velocity * this->vint_bessel_011_perp(vt, z, varrho, this->R);
	double i2_perp = velocity * this->vint_bessel_001_perp(vt, z, varrho, this->R);

	double res = m * A3 * chi * std::sqrt(NonlinearMedium::MU0);
	res /= 64 * std::sqrt(nu);
	res *= em_relation * em_relation * em_relation;
	res *= (i2-i1) * (i1_perp*(i2-i1) + 2*i1*(i2_perp-i1_perp));
	return res;
}

double KerrAmendment::N3 (int m, double nu, double ct, double varrho, double z) const
{
	double sigma = nl_medium->conductivity(ct, z);
	double eps_r = nl_medium->relative_permittivity(ct,z);
	double mu_r = nl_medium->relative_permeability(ct,z);
	double em_relation = 1;
	em_relation *= std::sqrt(NonlinearMedium::MU0 * mu_r);
	em_relation /= std::sqrt(NonlinearMedium::EPS0 * eps_r);

	mpf_class vt_z = mpf_class(ct * ct);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );

	double i1 = this->int_bessel_011(std::sqrt(vt_z.get_d()), varrho, this->R);

	double res = m * this->A0 * sigma * std::sqrt(NonlinearMedium::MU0);
	res /= 4 * std::sqrt(nu);
	res *= em_relation * i1 * jn(m,varrho*nu);
	return res;
}

double KerrAmendment::N4 (int m, double nu, double ct, double varrho, double z) const
{
	double A3 = this->A0 * this->A0 * this->A0;
	double chi = nl_medium->relative_permittivity(ct, z, 3);

	double eps_r = nl_medium->relative_permittivity(ct,z);
	double mu_r = nl_medium->relative_permeability(ct,z);
	double velocity = NonlinearMedium::C / (eps_r * mu_r);
	double em_relation = 1;
	em_relation *= std::sqrt(NonlinearMedium::MU0 * mu_r);
	em_relation /= std::sqrt(NonlinearMedium::EPS0 * eps_r);

	mpf_class vt_z = mpf_class(ct * ct);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );

	double vt = ct / std::sqrt(eps_r * mu_r);

	double i1 = this->int_bessel_011(std::sqrt(vt_z.get_d()), varrho, this->R);
	double i2 = this->int_bessel_001(std::sqrt(vt_z.get_d()), varrho, this->R);
	double i1_perp = velocity * this->vint_bessel_011_perp(vt, z, varrho, this->R);
	double i2_perp = velocity * this->vint_bessel_001_perp(vt, z, varrho, this->R);
	double bessel = jn(m-1,varrho*nu) - jn(m+1,varrho*nu);

	double res = A3 * chi * std::sqrt(NonlinearMedium::MU0) * std::sqrt(nu) / 128;
	res *= em_relation * em_relation * em_relation;
	res *= (i2 - i1) * (i2 - i1) * (i2_perp - i1_perp) * bessel;
	return res;
}

double KerrAmendment::N5 (int m, double nu, double ct, double varrho, double z) const
{
	double A3 = this->A0 * this->A0 * this->A0;
	double chi = nl_medium->relative_permittivity(ct, z, 3);

	double eps_r = nl_medium->relative_permittivity(ct,z);
	double mu_r = nl_medium->relative_permeability(ct,z);
	double velocity = NonlinearMedium::C / (eps_r * mu_r);
	double em_relation = 1;
	em_relation *= std::sqrt(NonlinearMedium::MU0 * mu_r);
	em_relation /= std::sqrt(NonlinearMedium::EPS0 * eps_r);

	mpf_class vt_z = mpf_class(ct * ct);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );

	double vt = ct / std::sqrt(eps_r * mu_r);

	double i1 = this->int_bessel_011(std::sqrt(vt_z.get_d()), varrho, this->R);
	double i2 = this->int_bessel_001(std::sqrt(vt_z.get_d()), varrho, this->R);
	double i1_perp = velocity * this->vint_bessel_011_perp(vt, z, varrho, this->R);
	double i2_perp = velocity * this->vint_bessel_001_perp(vt, z, varrho, this->R);
	double bessel = jn(m-1,varrho*nu) - jn(m+1,varrho*nu);
	double mult = i1 * (i2_perp - i1_perp) + 2 * i1_perp * (i2 - i1);

	double res = 3 * A3 * chi * std::sqrt(NonlinearMedium::MU0) * std::sqrt(nu) / 128;
	res *= em_relation * em_relation * em_relation;
	res *= mult * bessel;
	return res;
}

double KerrAmendment::N6 (int m, double nu, double ct, double varrho, double z) const
{
	double sigma = nl_medium->conductivity(ct, z);
	double eps_r = nl_medium->relative_permittivity(ct,z);
	double mu_r = nl_medium->relative_permeability(ct,z);
	double em_relation = 1;
	em_relation *= std::sqrt(NonlinearMedium::MU0 * mu_r);
	em_relation /= std::sqrt(NonlinearMedium::EPS0 * eps_r);

	mpf_class vt_z = mpf_class(ct * ct);
	mpf_div(vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(eps_r * mu_r).get_mpf_t());
	mpf_add( vt_z.get_mpf_t(), vt_z.get_mpf_t(), mpf_class(- z * z).get_mpf_t() );

	double i1 = this->int_bessel_011(std::sqrt(vt_z.get_d()), varrho, this->R);
	double i2 = this->int_bessel_001(std::sqrt(vt_z.get_d()), varrho, this->R);
	double bessel = jn(m-1,varrho*nu) - jn(m+1,varrho*nu);

	double res = this->A0 * sigma * std::sqrt(NonlinearMedium::MU0) * std::sqrt(nu) / 8;
	res *= em_relation * (i2 - i1) * bessel;
	return res;
}

/* */

double KerrAmendment::vint_bessel_011_perp (double vt, double z, double rho, double R) const
{
	double vt_z = vt * vt - z * z;
	double rho2 = rho * rho;
	double R2 = R * R;

	if (vt_z < 0) throw std::invalid_argument("ct-z < 0 is not legal");
	if (vt_z == 0) throw std::invalid_argument("ct-z = 0 is not legal");
	if (rho < 0) throw std::invalid_argument("rho < 0 is not legal");
	if (R <= 0) throw std::invalid_argument("R <= 0 is not legal");

	if (R < std::abs(rho - std::sqrt(vt_z))) return 0.0;
	if (R > rho + std::sqrt(vt_z)) return 0.0;

	double res = vt / rho2;
	res *= (rho2 - R2) * (rho2 - R2) / vt_z;
	res /= std::sqrt( (rho+R)*(rho+R) - vt*vt + z*z );
	res /= std::sqrt( (rho+R)*(rho+R) - vt*vt + z*z );
	return res / (2*M_PI);
}

double KerrAmendment::vint_bessel_001_perp (double vt, double z, double rho, double R) const
{
	double vt_z = vt * vt - z * z;
	double rho2 = rho * rho;
	double R2 = R * R;

	if (vt_z < 0) throw std::invalid_argument("ct-z < 0 is not legal");
	if (vt_z == 0) throw std::invalid_argument("ct-z = 0 is not legal");
	if (rho < 0) throw std::invalid_argument("rho < 0 is not legal");
	if (R <= 0) throw std::invalid_argument("R <= 0 is not legal");

	if (R < std::abs(rho - std::sqrt(vt_z))) return 0.0;
	if (R > rho + std::sqrt(vt_z)) return 0.0;

	double res = - vt / vt_z;
	res *= vt_z - rho2 + R2;
	res /= std::sqrt(4*rho2*vt_z - (vt_z+rho2-R2) * (vt_z+rho2-R2) );
	return res / M_PI;
}
