#
#  archive.sh
#  Maxwell
#
#  Created by Rolan Akhmedov on 23.04.18.
#  Copyright © 2018 Rolan Akhmedov. All rights reserved.
#

#!/bin/bash

mkdir -p archive/$1
mysqldump --databases maxwell -umaxwell -pmaxwell > archive/$1/maxwell.sql
mysql -umaxwell -pmaxwell < clean.sql
mv maxwell-*.log archive/$1/
mv maxwell.gnp archive/$1/
cp maxwell.conf archive/$1/
touch archive/$1/read_me.txt
