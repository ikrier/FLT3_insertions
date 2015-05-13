#!/usr/bin/python
# =======
#   License
# =======
#   This code is released under the GNU General Public License 3.0. A copy
# of this license is in the LICENSE.txt file.
# copyright Irina Krier 2015
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.



import re
import numpy
import csv
from StringIO import StringIO
from itertools import chain

fixedindices=[2, 5, 7, 12, 13, 26]
indices_samples=range(31,7*84-1+31,7)
indices_coverage=range(32,7*84-1+32,7)
indices_upstream_all=[x+2 for x in indices_coverage]
indices_upstream_uniq=[x+3 for x in indices_coverage]
indices_downstream_all=[x+4 for x in indices_coverage]
indices_downstream_uniq=[x+5 for x in indices_coverage]

libnames= ["Lib_"+`i` for i in range(1,85)]
libnames.sort()

covnames=["cov_"+i for i in libnames]
upallnames=["upall_"+i for i in libnames]
upunames=["upuniq_"+i for i in libnames]
downallnames=["downall_"+i for i in libnames]
downunames=["downuniq_"+i for i in libnames]


fileout=open("parsed_flt3_SIs_all_libs.csv","w")

print >>fileout, "length\tsequence\tchromosome\tstartrange\tendrange\tsum_scores\t"+"\t".join(list(chain(*zip(covnames,upallnames,upunames,downallnames,downunames))))

with open("check_insertion_flt3_all_libs_SI","r") as f:
	io = StringIO(f.read().replace('\t', ' '))
        reader = csv.reader(io, delimiter=" ")	
	for row in reader :
		if re.match('.*ChrID.*',".".join(row)) :
			if int(row[26])>0 :
				row=numpy.array(row)
				print >>fileout, "\t".join(row[fixedindices])+"\t",
				print >>fileout, "\t".join(row[list(chain(*zip(indices_coverage,indices_upstream_all,indices_upstream_uniq,indices_downstream_all,indices_downstream_uniq)))])
fileout.close()
