//
//  module_manager.cpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "module_manager.hpp"

ModuleManager::ModuleManager (Logger* global)
: log(global) {
    // load_module("ZeroFiled", {new HomoMedium(), new NoSource(), new Zero()});
    // load_module("WhiteGaussianNoise", {new HomoMedium(), new NoSource(), new Noise()});
}

bool ModuleManager::load_module (const std::string& posix_path, double R, double A0, double tau0, double eps, double mu)
{
    std::ifstream lib(posix_path, std::ios::binary);
    if (!lib) {
        if (log) log->error("File does not exist: " + posix_path);
        return false;
    }

    void* library = dlopen(posix_path.data(), RTLD_LAZY);
    if (!library) {
        if (log) log->error("File is not a library: " + posix_path);
        return false;
    }

    LoaderFncPtr call = reinterpret_cast<LoaderFncPtr> (dlsym(library, "load_module"));
    if (!call) {
        std::string messeng = "load_module() function have not been found in module library: ";
        if (log) log->error(messeng + posix_path);
        return false;
    }

    if (log) log->info("Module library loaded: " + posix_path);
    call(this, log, R, A0, tau0, eps, mu);
    return true;
}

bool ModuleManager::load_module (const std::string& name, ModuleEntity impl)
{
    auto status = module.insert(std::make_pair(name, impl));
    if (log) {
        log->info("Submodule loaded: " + name);
        if (!status.second) log->warning("Loaded module already existed: " + name);
    }
    return status.second; // false if already existing
}

std::vector<std::string> ModuleManager::get_loaded () const
{
    std::vector<std::string> list;
    for (auto &pair : module)
        list.push_back(pair.first);
    return list;
}

ModuleEntity ModuleManager::get_module (const std::string& name)
{
    auto it = module.find(name);
    if (it == module.end())
        return {nullptr, nullptr, nullptr};
    return it->second;
}
