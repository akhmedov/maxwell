//
//  integral.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 17.08.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef integral_hpp
#define integral_hpp

#include <functional>
#include <cmath>
#include <vector>

struct Integral { };

struct Simpson : public Integral {
	Simpson (std::size_t terms);
	double value (double from, double to, std::function<double(double)> f);
private:
	const std::size_t quadr_terms;
};

struct GaussLaguerre : public Integral {
	GaussLaguerre ();
	GaussLaguerre (std::size_t terms);
	double value (std::function<double(double)> f);
private:
	const std::vector<std::pair<double, double>> polynom_data {
		/*				root				weight				*/
		std::make_pair( 0.0933078120170,	0.218234885940),
		std::make_pair( 0.4926917403020,	0.342210177923),
		std::make_pair( 1.2155954120710,	0.263027577942),
		std::make_pair( 2.2699495262040,	0.126425818106),
		std::make_pair( 3.6676227217510,	0.040206864921),
		std::make_pair( 5.4253366274140,	0.856387780361e-2),
		std::make_pair( 7.5659162266130,	0.121243614721e-2),
		std::make_pair( 10.120228568019,	0.111674392344e-3),
		std::make_pair( 13.130282482176,	0.645992676202e-5),
		std::make_pair( 16.654407708330,	0.222631690710e-6),
		std::make_pair( 20.776478899449,	0.422743038498e-8),
		std::make_pair( 25.623894226729,	0.392189726704e-10),
		std::make_pair( 31.407519169754,	0.145651526407e-12),
		std::make_pair( 38.530683306486,	0.148302705111e-15),
		std::make_pair( 48.026085572686,	0.160059490621e-19)
	};
	const std::size_t quadr_terms;
};

#endif /* integral_hpp */
