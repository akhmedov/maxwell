--
--  clean.sql
--  Maxwell
--
--  Created by Rolan Akhmedov on 19.01.18.
--  Copyright © 2017 Rolan Akhmedov. All rights reserved.
--

-- delete single value 
-- DELETE FROM maxwel_data WHERE ct_from <= 0.1 AND z_from <= 0.1;

-- delete all values related to problem
-- DELETE FROM maxwell_header WHERE radius = 1;

USE maxwell;
DELETE FROM maxwell_header;
