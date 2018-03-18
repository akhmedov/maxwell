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
	component	ENUM('Ex', 'Ey', 'Ez', 'Ephi', 'Erho', 'Hx', 'Hy', 'Hz', 'Hphi', 'Hrho') NOT NULL,
	radius		DOUBLE(20,7) NOT NULL,
	magnitude	DOUBLE(20,7) NOT NULL,
	mu_r		DOUBLE(20,7) NOT NULL,
	eps_r		DOUBLE(20,7) NOT NULL,
	kerr_r		DOUBLE(20,7) NOT NULL,
	duration	DOUBLE(20,7) NOT NULL,
	signal_type	ENUM('turn_on','turn_off','rect','gauss','triangle') NOT NULL,
	PRIMARY KEY (id),
	CONSTRAINT UC_maxwell_header UNIQUE (radiator, component, radius, 
		magnitude, mu_r, eps_r, kerr_r, duration, signal_type)
) ENGINE=InnoDB;

-- slave table headding

CREATE TABLE IF NOT EXISTS maxwell_data (
	id		BIGINT NOT NULL UNIQUE AUTO_INCREMENT,
	head_id	BIGINT NOT NULL,
	ct		DOUBLE(20,7) NOT NULL,
	rho		DOUBLE(20,7) NOT NULL,
	phi		DOUBLE(20,7) NOT NULL,
	z		DOUBLE(20,7) NOT NULL,
	lineary	DOUBLE(20,7) DEFAULT NULL,
	square	DOUBLE(20,7) DEFAULT NULL,
	kerr	DOUBLE(20,7) DEFAULT NULL,
	PRIMARY KEY (id),
	FOREIGN KEY (head_id) REFERENCES maxwell_header (id) ON DELETE CASCADE,
	CONSTRAINT UC_maxwell_data UNIQUE (head_id, ct, rho, phi, z)
) ENGINE=InnoDB;
