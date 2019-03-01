//
//  module_manager.hpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#ifndef module_manager_hpp
#define module_manager_hpp

#include "linear_medium.hpp"
#include "linear_source.hpp"
#include "abstract_field.hpp"

#include <dlfcn.h>
#include <utility>
#include <string>
#include <iostream>
#include <map>

using CStrFncPtr   = const char* (*)();
using SourceFncPtr = LinearCurrent* (*)(); // TODO: pass matrix (castom type)
using MediumFncPtr = LinearMedium* (*)();
using FieldFncPtr  = AbstractField* (*)();

struct ModuleEntity {
    LinearCurrent* source; // TODO: use union {LinearCurrent, LinearCharge}
    LinearMedium* medium;
    AbstractField* field;
};

class ModuleManager {

public:
    ModuleManager ();
    bool load_module (const std::string& posix_path);
    void load_module (const std::string& name, ModuleEntity implementation);

    std::vector<std::string> get_loaded () const;
    ModuleEntity get_module (const std::string& name);

private:
    std::map<std::string, ModuleEntity> module;
};

#endif /* module_manager */
