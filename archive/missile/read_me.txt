BUILD OPTIONS
make DEBUG=1 -j4
17f77eb rolan.kharkiv@gmail.com Safe event logger

rolan@laptop:maxwell(master)$ time ./build/maxwell --safe --log --plot 1 --cylindric 1.7:0.05:3,0,0,2
Configuration of runtime...
World point (m): 1.7:0.05:3, 0, 0, 2
Source radius (m): 1
Source magnitude (m): 1
Relative epsilon: 1
Relative mu: 1
Kerr medium coefficient: 0
Conductivity (sigma): 0
Noise level: type(WHITE_GAUSS), power(0) mu(0), sigma(0.000000)
Magnetic component terms number: 263
Float bitrae of GMP: 512
Calculation thread number: 4
MySQL client will try to connect to maxwell@localhost...
Model debug informaition will be loged at maxwell.log...
Evaluation progress: 86%
Write comands to maxwell.gnp script... Done.

real	0m0.747s
user	0m0.421s
sys		0m0.012s
