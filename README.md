ABSTRACT
======

Maxwell is high performance software unit for numerical simulation of 
ultra-wideband radar and telecommunication systems. It supports FDTD 
and numerical evolutionary approach for simulation of electromagnetic 
processes defined by an arbitrary source of energy. Alternatively, it 
provides an interface to compute analytical obtained electromagnetic 
field defined by field components like Ex, Ey and Ez. These are a set 
of predefined analytical solutions for these kind of problems.
The product is based on Meep and GNUplot open source projects.

Visualization of the simulation is done with GNU Plot.
There are several options to save computed data. It can be saved as 
GNU Plot script, MySQL database script or JSON time series dataset. 
Moreover, the project is using MySQl C++ driver to store computing data
for future call, so the second call of the Maxwell for the same problem 
will use pre-computed data and will be much more faster.

Source code is wrote within ISO c++17 standard and POSIX compatibility. 
Build system is only one Makefile with several available options for 
testing and debugging. Maxwell project provides several interfaces for 
high performance parallel calculations, they are: CPU multi-threading by
std::threads, GPU calculation by CUDA and Message Passing Interface (MIP)
as a legacy of Meep project.

Source code is OS free and build systems supports GNU/Linux, MacOS and 
BSD operation systems. The code can be compiled for Windows-like systems 
with the help of Linux environment for Windows like MinGW or Windows 
Subsystem for Linux. Following OS list is recommended for development 
process, the system is tested on the list:

- GNU/Linux Ubuntu 16.04 LTE x86_64 (GNU GCC version 5.4.0 20160607)
- MacOS Sierra 10.12.5 x86_64 (Apple LLVM version 8.1.0)

BUILD INSTRUCTIONS FOR UBUNTU 16.04 x64
======

m4                 - dependences for gmp <br/>
libx11-dev         - dependences for gnuplot <br/>
libmysqlclient-dev - mysql client connector <br/>

At the first, it is required to install dependences

```bash
~ $ sudo apt-get install mysql-server
~ $ mysql_secure_installation
~ $ sudo apt-get install m4 libx11-dev libmysqlclient-dev
~ $ sudo apt-get install build-essential cmake git
```

Now it is possible to initialize database in MySQL server

```bash
~ $ mysql -t -u root -p < setup.sql
```

The next step is getting the source from VCS

```bash
~ $ git clone https://server/path/maxwell.git
```

Before compiling the protects you need to compile some dependences. Do not use 
`-j` option for this - it will not speed up compilation

```bash
~ $ make gnuplot gnump meep
```

Finally, go to the source directory and compile the project and unit tests for it

```bash
~ $ cd maxwell
maxwell $ make -j4 unit_test maxwell
maxwell $ ./build/unit_test
```

Now, run the example to examine the build

```bash
maxwell $ ./build/maxwell --help
maxwell $ ./build/maxwell --version
maxwell $ ./build/maxwell --noise 20 --plot 3 --cylindric 1.9:0.01:3,0:0.01:3,0,2
```

Also, you can use next commands to clean MySQL and source directory

```bash
maxwell $ make clean
maxwell $ mysql -u root -p < purge.sql
```

USEFULL DEBUG COMMANDS
======

```bash
~ $ make list

~ $ du -h ./build/maxwell
~ $ ls -lh ./build/maxwell
~ $ file ./build/maxwell
~ $ stat -x ./build/maxwell
~ $ nm -D libgmp.so | grep __gmpf_init

~ $ cat <> maxwell.conf | sed '/\(^FLOAT_BITRATE =\).*/ s//\1$(TERMS)/'
~ $ date && nice -n 10 ./build/gnuplot --conf PATH --model NUM && date
```

FAQ FOR MAXWELL'S DEPENDENCIES
======

gnump library
------

https://gmplib.org <br/>
./configure --help

GNU Plot library
------

http://plplot.sourceforge.net/ <br/>
http://www.gnuplot.info <br/>
http://www.gnuplot.info/demo/surface1.html <br/>
http://gnuplot.sourceforge.net/demo/pm3d.html <br/>


```
set term x11 font "times-roman,15,normal"
set ylabel "" font font "Times-New-Roman,15"
set tics font "Times-New-Roman,15"
set key font "Times-New-Roman,15"
set ylabel "f(x) = J_1(ax) J_0(bx) J_0(cx)" offset 0,-5
```

QT support for GNU Plot library
------

```bash
./configure --with-qt --prefix=$(PROJECT_DIR)/gnuplot
sudo apt-get install qtbase5-dev libqt5svg5-dev
```

Working around the MySQL
------

http://www.mysqltutorial.org/understand-mysql-table-types-innodb-myisam.aspx <br/>
https://dev.mysql.com/doc/connector-cpp/en/connector-cpp-installation-source-unix.html <br/>

```bash
~ $ sudo apt-get install mysql-server
~ $ mysql_secure_installation
~ $ mysql -t -u root -p < setup.sql
~ $ mysql -vvv -u root -p < clean.sql
~ $ mysql -u root -p 
~ $ mysqldump -u root -p maxwell > maxwell.sql
```

```sql
mysql> CREATE DATABASE maxwell;
mysql> GRANT ALL PRIVILEGES ON *.* TO 'maxwell'@'localhost' IDENTIFIED BY 'maxwell';
mysql> SELECT User, Host, Password FROM mysql.user;
mysql> DROP USER 'maxwell'@'localhost';
mysql> DROP TABLE IF EXISTS maxwell;
mysql> SELECT user FROM mysql.user GROUP BY user;
mysql> DELETE FROM mysql.user WHERE user = 'maxwell';
mysql> SHOW COLUMNS FROM maxwell.maxwell_header;
mysql> SHOW COLUMNS FROM maxwell.maxwell_data;
```
