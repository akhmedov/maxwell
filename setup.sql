--
--  setup.sql
--  Maxwell
--
--  Created by Rolan Akhmedov on 19.01.18.
--  Copyright © 2017 Rolan Akhmedov. All rights reserved.
--

-- user and database init

CREATE DATABASE IF NOT EXISTS maxwell;
CREATE USER IF NOT EXISTS maxwell;
GRANT ALL PRIVILEGES ON maxwell.* TO 'maxwell'@'localhost' IDENTIFIED BY 'maxwell';
USE maxwell;

-- master table headding

CREATE TABLE IF NOT EXISTS maxwell_header (
	id			BIGINT NOT NULL UNIQUE AUTO_INCREMENT,
	radiator	ENUM('uni_disk', 'grad_disk', 'radial_disk') NOT NULL,
	component	ENUM('Ex','Ey','Ez','Hx','Hy','Hz') NOT NULL,
	radius		DOUBLE(20,5) NOT NULL,
	magnitude	DOUBLE(20,5) NOT NULL,
	mu_r		DOUBLE(20,5) NOT NULL,
	eps_r		DOUBLE(20,5) NOT NULL,
	kerr_r		DOUBLE(20,5) NOT NULL,
	duration	DOUBLE(20,5) NOT NULL,
	signal_type	ENUM('turn_on','turn_off','rect','gauss','triangle') NOT NULL,
	PRIMARY KEY (id),
	CONSTRAINT UC_maxwell_header UNIQUE (radiator, component, radius, 
		magnitude, mu_r, eps_r, kerr_r, duration, signal_type)
) ENGINE=InnoDB;

-- master table entity example

/* INSERT INTO maxwell_header SET 
	radiator  = 'uni_disk',
	component = 'Ex',
	radius 	  = 1,
	magnitude = 1,
	mu_r      = 1,
	eps_r     = 1,
	kerr_r    = 1,
	duration  = 0,
	signal_type = 'rect'
; */

-- slave table headding

CREATE TABLE IF NOT EXISTS maxwell_data (
	id		BIGINT NOT NULL UNIQUE AUTO_INCREMENT,
	head_id	BIGINT NOT NULL,
	ct		DOUBLE(20,5) NOT NULL,
	rho		DOUBLE(20,5) NOT NULL,
	phi		DOUBLE(20,5) NOT NULL,
	z		DOUBLE(20,5) NOT NULL,
	noise	DOUBLE(20,5) DEFAULT NULL,
	lineary	DOUBLE(20,5) DEFAULT NULL,
	square	DOUBLE(20,5) DEFAULT NULL,
	kerr	DOUBLE(20,5) DEFAULT NULL,
	PRIMARY KEY (id),
	FOREIGN KEY (head_id) REFERENCES maxwell_header (id) ON DELETE CASCADE,
	CONSTRAINT UC_maxwell_data UNIQUE (ct, rho, phi, z)
) ENGINE=InnoDB;

-- slave table entity example

/* SELECT @problem := id FROM maxwell_header WHERE 
radiator = 'uni_disk' AND component = 'Ex' AND radius = 1 AND 
magnitude = 1 AND mu_r = 1 AND eps_r = 1 AND kerr_r = 1 AND 
duration = 0 AND signal_type = 'rect'; */

/* INSERT INTO maxwell_data SET
	-- head_id = LAST_INSERT_ID(),
	head_id	= @problem,
	ct		= 0,
	rho		= 0,
	phi		= 0,
	z		= 0,
	noise	= NULL,
	lineary	= NULL,
	square	= NULL,
	kerr	= NULL
; */

-- get value set from database

/* select ct,lineary,kerr from maxwell_data where z=0, rho=0; */
