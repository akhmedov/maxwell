KNOWN BUGS AND NEW FEATURES
======

TODO LIST
------

01. (Cource, Medium, AbstructField) -> CylindricalField<Space> : Field<Space>
02. Magnetic Field calculation
03. GnuPlot -> ScriptManager<GnuPlot>
04. nlohman::json -> boost::property_tree::write_json
05. GMP -> boost::multiprecision
06. Noise class code revive
07. Extend core coverage and refactor of tests
08. SafeManager calculation on shared database
09. Implement 1D and 2D Point::SpaceTime<System>
10. PyPlot -> ScriptManager<PyPlot>
11. GUI of feco, hfss, cst to Maxwell
12. Spectrum calculations in AbstructField<>
13. Add test for MySQL API for evolution coefficient table

NEW MODULES TO IMPLEMENT
------

1. Dipole antenna
2. Bow tie antenna
3. Non-plane bow tie antenna
4. Butt of coaxial cable

BUGS LIST
------

1. fix up_disk precompiler define-vars
2. warnings on GNU GCC 5.4.0 20160607
3. MySQL client select_point execution time
4. Initialization list constractor for Point::SpaceTime<System>
5. override keyword for CartesianField, CylindricalField
6. FFT frequency resolution and magnitude bugs
7. LinearDuhamel<> inegtal exaption returns 0 - it is not correct
8. CalculationManagers divide argument for separate threads and some of them are done faster. Then dead cores are not in use.

UPDATE LIST FOR c++17 STANDART
------

- paralel algorithms
- use std::variant receiver evants
- use std::filsistems::exists() for file check in Dataset:: and GnuPlot::
- use std::cyl_bessel_i to replace usage of depricated functon jn(m,z) from <cmath>
