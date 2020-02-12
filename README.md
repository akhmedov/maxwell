MAXWELL PROJECT
======

ABSTRACT
------

```diff
- WARNING: The project is in the public beta state!
- Some of unit tests fail and there are a lot of uncovered features
- Compilation warnings and "not implemented" exaptions are presented
- The first ready to use release is planning for
- December 31, 2020
```

Maxwell is a library for high performance numerical simulation of
ultra-wideband radar and telecommunication systems. It implements FDTD,
evolutionary approach and other methods of time domain electrodynamics.
Build-in functionality can by used to obtain properties of EM filed
radiated by arbitrary user defined source. Alternatively, the library
provides an interface AbstractFiled to use build-in functions to compute
user defined field components (\vect{E} and \vect{H}).

Maxwell project provides universal and simple to use interface for high
performance parallel calculations on user defined problems that can
continue unexpectedly interrupted calls from previous state. Also, it
allows to use previously calculated points on new call.

Visualization of the simulations and computations is provided throw
saving plot data in GNUPlot or Python script that is totally independent
from current project.

On the other side the library allows to save computed data as
machine learning dataset in JSON or CSV format. The output dataset
is a time series with menadta to each element. It is desined to help
to embed morden datascience methods into time domain electronics.

SOURCE CODE INFORMAITION
------

Source code is wrote within ISO c++14 standard and POSIX compatibility.
Build system consist of the hierarchy of cmake files that provides
automatic building, testing and installing of the library. It is
expected that your system supports c++14 toolchain and has MySQL
client with c++ API.

Core functionality of the library is located in `maxwell/core` dir and
`libmaxwell.so` in build. Solution for a spesific problems like certan
type of antenna presents in dirictry `maxwell/module` as chiled projects
that are dependent on `libmaxwell.so`. Any module is built as a separat
library that can be loaded in runtime by `ModuleManager` (`libmaxwell.so`)
with using of [dlfcn](https://pubs.opengroup.org/onlinepubs/7908799/xsh/dlfcn.h.html).

All modules and core library have `test` dirictory of unit-test cpp-files.
The tests are writn in `ctest` format.

Source code is OS free and build systems supports GNU/Linux, MacOS and
BSD operation systems. The code can be compiled for Windows-like systems
with the help of Linux environment for Windows like MinGW or Windows
Subsystem for Linux. Following OS list is recommended for development
process, the system is tested on the list:

- GNU/Linux Ubuntu 16.04 LTE x86_64 (GNU GCC version 5.4.0 20160607)
- GNU/Linux Ubuntu 18.04 LTE x86_64 (GNU GCC version 7.3.0 20180607)
- MacOS Sierra 10.12.5 x86_64 (Apple LLVM version 8.1.0)
- MacOS Sierra 10.14.3 x86_64 (Apple LLVM version 10.0.0)

BUILDING INSTRUCTIONS FOR UBUNTU 18.04 x64
------

If you are building the library on in other environment follow the wiki,
nevertheless it is strongly recommended to read this tutorial firstly.

At first making sure that your c++ toolchain is ready. A list you need
`gcc`, `g++`, `make`, `cmake`, `ctest` and `git` tools. Also, your system
must suport POSIX standart for dynamic library linking in runtime
[dlfcn](https://pubs.opengroup.org/onlinepubs/7908799/xsh/dlfcn.h.html).
These and some more debuging tools can be installed by apt-get:

```bash
~ $ sudo apt-get install build-essential cmake git
```

On this step you are able to clone or download sources:

```bash
~ $ git clone --depth=1 https://server/path/maxwell.git
```

If you are not able to login into MySQL server with `mysql -u root -p`
then install and configure MySQL server. Also, maxwell project
is using c++ API of MySQL that is provided in `libmysqlclient-dev` package

```bash
~ $ sudo apt-get install mysql-server
~ $ mysql_secure_installation
~ $ sudo apt-get install libmysqlclient-dev
```

If your MySQL server is ready create user and database for maxwell project:

```bash
~ $ mysql -t -u root -p < maxwell/setup.sql
```

Also, you can clean the server from maxwell data, user and tables by

```bash
~ $ mysql -t -u root -p < maxwell/purge.sql
```

The library maxwell is depended on POSIX library [GMP](https://gmplib.org).
It is recommended to build its by yourself, so you will need to install its
dependency and run the make script (or learn more by ./configure --help):

```bash
~ $ sudo apt-get install m4
~ $ make -C maxwell gnump
~ $ ls maxwell/gnump/include/ maxwell/gnump/lib/
```

Maxwell projects does not have plot functions. It saves plot data in the `*.gnp`
script, that is supported separately by [gnuplot project](http://www.gnuplot.info)
so you need to install `gnuplot` package (learn more ./configure --help):

```bash
~ $ sudo apt-get install gnuplot
~ $ gnuplot --version
```

Alternatively, you can build it by yourself:

```bash
~ $ sudo apt-get install libgd-dev libx11-dev
~ $ make -C maxwell gnuplot
~ $ maxwell/gnuplot/bin/gnuplot --version
```

Finally, you are ready to build the library. For this run the
following commands in Maxwell directory:

```bash
~ $ mkdir maxwell/build && cd maxwell/build
build $ cmake ..
build $ make -j5
```

After secsessful build you will find shared object `build/core/libmaxwell.so`
and dynamicly linking libraries `build/MODNAME/libMODNAME.so` for all modules.
Note, library extantion depends from your system: `.so` for GNU/Linux,
`.dylib` for Mac, and `.dll` for Windows.

It is recommended to run tests for successful build. Brouse build dirictory
for `test` subdirictoris and run `ctest` for all of them. Example:

```bash
~ $ cd build/core/test && ctest
```

MAXWELL PROJECT CONTRIBUTORS
------

Name                | Contribution
------------------- | --------------------------------------------
Rolan Akhmedov      | project director <br>
Oleksandr Dumin     | consultation in applied physics methods <br>
Oleh Zahrychanskyi  | build system and multithreading <br>

LIST OF USEFULL COMMANDS
------

Running mysql script in verbose mode

```bash
~ $ mysql -vvv -u root -p < setup.sql
```

Saving database to script

```bash
~ $ mysqldump -u root -p maxwell > maxwell.sql
```

```bash
~ $ stat -x ./build/maxwell
```

```bash
~ $ file ./build/maxwell
```

```bash
~ $ nm -D libgmp.so | grep __gmpf_init
```

```bash
~ $ otools -L example/plot_energy_distribution
```
