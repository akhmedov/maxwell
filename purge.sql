--
--  purge.sql
--  Maxwell
--
--  Created by Rolan Akhmedov on 19.01.18.
--  Copyright © 2017 Rolan Akhmedov. All rights reserved.
--

/* USE maxwell;
ALTER TABLE maxwell_header DROP PRIMARY KEY;
DROP TABLE IF EXISTS maxwell_data;
DROP TABLE IF EXISTS maxwell_header; */

DROP DATABASE maxwell;
DROP USER 'maxwell'@'localhost';
DELETE FROM mysql.user WHERE user = 'maxwell';
