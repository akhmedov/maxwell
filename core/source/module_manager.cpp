#include "module_manager.hpp"

ModuleManager::ModuleManager ()
{
    // load_module("ZeroFiled", {new HomoMedium(), new NoSource(), new Zero()});
    // load_module("WhiteGaussianNoise", {new HomoMedium(), new NoSource(), new Noise()});
}

bool ModuleManager::load_module (const std::string& posix_path)
{
    void* module = dlopen(posix_path.data(), RTLD_LAZY);
    if (!module) return false; // TODO: log event
    CStrFncPtr name = reinterpret_cast<CStrFncPtr> (dlsym(module, "load_module_name"));
    if (!name) return false; // TODO: log event
    SourceFncPtr source = reinterpret_cast<SourceFncPtr> (dlsym(module, "load_module_source"));
    if (!source) return false; // TODO: log event
    MediumFncPtr medium = reinterpret_cast<MediumFncPtr> (dlsym(module, "load_module_medium"));
    if (!medium) return false; // TODO: log event
    FieldFncPtr field = reinterpret_cast<FieldFncPtr> (dlsym(module, "load_module_filed"));
    if (!field) return false; // TODO: log event
    load_module(std::string(name()), {source(), medium(), field()});
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
