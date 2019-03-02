//
//  module_manager.hpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#ifndef module_manager_hpp
#define module_manager_hpp

#include "logger.hpp"
#include "linear_medium.hpp"
#include "linear_source.hpp"
#include "abstract_field.hpp"

#include <dlfcn.h>
#include <utility>
#include <string>
#include <iostream>
#include <map>

struct ModuleEntity {
    LinearCurrent* source; // TODO: use union {LinearCurrent, LinearCharge}
    LinearMedium* medium;
    AbstractField* field;
};

class ModuleManager {

public:
    ModuleManager ();
    bool load_module (const std::string& posix_path, Logger* global, double R, double A0, double tau0, double eps, double mu);
    void load_module (const std::string& name, ModuleEntity implementation);

    std::vector<std::string> get_loaded () const;
    ModuleEntity get_module (const std::string& name);

private:
    std::map<std::string, ModuleEntity> module;
};

using LoaderFncPtr = void (*)(ModuleManager*,Logger*,double,double,double,double,double);

#endif /* module_manager */
