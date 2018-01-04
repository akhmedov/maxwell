//
//  plot_model.hpp
//  Maxwell
//
//  Created by Rolan Akhmedov on 27.12.17.
//  Copyright © 2017 Rolan Akhmedov. All rights reserved.
//

#ifndef plot_model_hpp
#define plot_model_hpp

#include "uniform_disk_current.hpp"
#include "gnu_plot.hpp"
#include "manager.hpp"
#include "integral.hpp"
#include "config.hpp"
#include "noise.hpp"

#include <algorithm>
#include <vector>
#include <functional>
#include <tuple> // std::make_tuple
#include <limits> // std::numeric_limits<size_t>::max

struct PlotModel {
	PlotModel (Config* global_conf);
	enum Name { Ex_from_ct = 1, Hy_from_ct, Ex_from_ct_rho };
	void call (const Name& model_name);
private:
	Config* global_conf;
	std::vector<std::function<void(PlotModel&)>> model_pointer;
	void __Ex_from_ct ();
	void __Hy_from_ct ();
	void __Ex_from_ct_rho ();
};

#endif /* plot_model_hpp */
