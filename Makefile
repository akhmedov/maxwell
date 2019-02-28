#
#  Makefile
#  Maxwell
#
#  Created by Rolan Akhmedov on 08.07.17.
#  Copyright © 2017 Rolan Akhmedov. All rights reserved.
#

PROJECT_DIR = $(shell pwd)

gnuplot:
	tar xf archive/gnuplot-5.0.6.tar
	mkdir -p gnuplot && cd gnuplot-5.0.6; \
	./configure --with-x --without-qt --prefix=$(PROJECT_DIR)/gnuplot; \
	make && make install; \
	rm -fr $(PROJECT_DIR)/gnuplot-5.0.6

gnump:
	tar xvjf archive/gmp-6.1.2.tar.bz2
	mkdir -p gnump && cd gmp-6.1.2; \
	./configure --enable-cxx --prefix=$(PROJECT_DIR)/gnump; \
	make && make check && make install; \
	rm -fr $(PROJECT_DIR)/gmp-6.1.2
