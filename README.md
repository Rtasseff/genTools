genTools
========

A bunch of python2.7 scripts and modules to do some common tasks, very messy

modules
-------------------
The most valuable:
statsUtil, contains some of my commonly used statistic tools.
stand alone other than:
import numpy as np
import scipy.stats as stats
import scipy.linalg
import operator as op
import itertools

dataUtil, some occasionally valuable data operations 
stand alone other than:
import numpy as np 
from operator import itemgetter

manageData, useful for microarray data, but its customized more my
own purpose and not well documented or commented.

cvTools, was valuable for me, but I then found that most
of this can be implemented easier by scipy and numpy built 
in functions.


license
-------------------

Copyright (C) 2003-2013 Institute for Systems Biology
		     Seattle, Washington, USA.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

20120904 RAT
