//
//  module_manager.cpp
//  core.maxwell
//
//  Created by Rolan Akhmedov on 28.02.19
//  Copyright © 2019 Rolan Akhmedov. All rights reserved.
//

#include "module_manager.hpp"

const std::vector<std::string> ModuleManager::format = {"dylib","so","dll"};

ModuleManager::ModuleManager (Logger* global)
: log(global) {
    // load_module("ZeroFiled", {new HomoMedium(), new NoSource(), new Zero()});
    // load_module("WhiteGaussianNoise", {new HomoMedium(), new NoSource(), new Noise()});
}

bool ModuleManager::load_module (const std::string& posix_path, const std::string& libname, double R, double A0, double tau0, double eps, double mu)
{
    std::string relative = this->find_lib (posix_path, libname);
    if (relative.empty()) return false;

    std::ifstream lib(relative, std::ios::binary);
    if (!lib) {
        if (log) log->error("File does not exist: " + relative);
        return false;
    }

    void* library = dlopen(relative.data(), RTLD_LAZY);
    if (!library) {
        if (log) log->error("File is not a library: " + relative);
        return false;
    }

    LoaderFncPtr call = reinterpret_cast<LoaderFncPtr> (dlsym(library, "load_module"));
    if (!call) {
        std::string messeng = "load_module() function have not been found in module library: ";
        if (log) log->error(messeng + relative);
        return false;
    }

    if (log) log->info("Module library loaded: " + relative);
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
        return {nullptr, nullptr, nullptr, nullptr, nullptr};
    return it->second;
}

std::string ModuleManager::find_lib (const std::string& posix_path, const std::string& library)
{
    std::string relative;

    for (auto i : format) {

        std::string tmp = posix_path + "/lib" + library + "." + i;
        std::ifstream lib(tmp, std::ios::binary);

        if (lib) {
            relative = tmp;
            break;
        }
        
        else if (log) log->warning("File does not exist: " + posix_path);
    }

    return relative;
}
