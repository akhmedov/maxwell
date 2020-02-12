--
--  setup.sql
--  Maxwell
--
--  Created by Rolan Akhmedov on 19.01.18.
--  Copyright © 2018 Rolan Akhmedov. All rights reserved.
--

-- user and database init

CREATE DATABASE IF NOT EXISTS maxwell;
CREATE USER IF NOT EXISTS maxwell;
GRANT ALL PRIVILEGES ON maxwell.* TO 'maxwell'@'localhost' IDENTIFIED BY 'maxwell';
USE maxwell;

-- master table headding

CREATE TABLE IF NOT EXISTS problem (
	id          BIGINT NOT NULL UNIQUE AUTO_INCREMENT,
	comment     CHAR(128) NOT NULL,
	PRIMARY KEY (id)
) ENGINE=InnoDB;

-- slave table headding

CREATE TABLE IF NOT EXISTS probe (
	id       	BIGINT NOT NULL UNIQUE AUTO_INCREMENT,
	problem_id  BIGINT NOT NULL,
	ct       	DOUBLE(20,7) NOT NULL,
	rho      	DOUBLE(20,7) NOT NULL,
	phi      	DOUBLE(20,7) NOT NULL,
	z        	DOUBLE(20,7) NOT NULL,
	result  	DOUBLE(20,7) DEFAULT NULL,
	PRIMARY KEY (id),
	FOREIGN KEY (problem_id) REFERENCES problem (id) ON DELETE CASCADE,
	CONSTRAINT unique_probe UNIQUE (problem_id, ct, rho, phi, z)
) ENGINE=InnoDB;

-- sub-slave table headding

CREATE TABLE IF NOT EXISTS evolution (
	id       	BIGINT NOT NULL UNIQUE AUTO_INCREMENT,
	problem_id  BIGINT NOT NULL,
	probe_id  	BIGINT NOT NULL,
	m 			INTEGER(3) NOT NULL,
	nu       	DOUBLE(25,7) NOT NULL,
	vmh 	  	DOUBLE(25,7) DEFAULT NULL,
	PRIMARY KEY (id),
	FOREIGN KEY (problem_id) REFERENCES problem  (id) ON DELETE CASCADE,
	FOREIGN KEY (probe_id)   REFERENCES probe (id) ON DELETE CASCADE,
	CONSTRAINT unique_coeff UNIQUE (problem_id, probe_id, m, nu)
) ENGINE=InnoDB;
