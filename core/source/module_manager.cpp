#include "module_manager.hpp"

ModuleManager::ModuleManager ()
{
    // load_module("ZeroFiled", {new HomoMedium(), new NoSource(), new Zero()});
    // load_module("WhiteGaussianNoise", {new HomoMedium(), new NoSource(), new Noise()});
}

bool ModuleManager::load_module (const std::string& posix_path, Logger* global, double R, double A0, double tau0, double eps, double mu)
{
    void* module = dlopen(posix_path.data(), RTLD_LAZY);
    if (!module) return false; // TODO: log event

    LoaderFncPtr call = reinterpret_cast<LoaderFncPtr> (dlsym(module, "load_module"));
    if (!call) return false; // TODO: log event

    call(this, global, R, A0, tau0, eps, mu);

    return true;
}

void ModuleManager::load_module (const std::string& name, ModuleEntity impl)
{
    auto addition = module.insert(std::make_pair(name, impl)); // TODO: log event
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
