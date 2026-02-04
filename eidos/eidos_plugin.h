//
// eidos_plugin.h


// Handle loading libraries
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#include <memory>

// Struct to close the lib automatically - saves a pointer of memory
// struct LibClose
// {
//     void operator()(void* handle) 
//     {
// #ifdef _WIN32
//         FreeLibrary((HINSTANCE)handle);
//         return;
// #endif
//         dlclose(handle);
//     }
// };

// using lib_ptr = std::unique_ptr<void, LibClose>;

// int (*ODE_fn_ptr)(double, const double*, double*, void*);
// int (*jac)(double, const double*, double*, double*, void*);

// using ODE_func_t = std::add_pointer_t<int(double, const double*, double*, void*)>;

// lib_ptr lib = LoadDynLib("lib.so");

// ODE_func_t ODE = LoadSymFromLib<ODE_fn_ptr>(lib, std::string("ODE"));