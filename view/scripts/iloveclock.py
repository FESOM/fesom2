#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Prints out all fesom.clock values for provided restart file.

Usage:

    python iloveclock.py /path/to/the/fesom.year.oce.restart.nc

by default "noleap" or "365_day" calendar will be used (COREII standard).

If you want to change calendar, provide it as second comand line argument:

    python iloveclock.py /path/to/the/fesom.year.oce.restart.nc proleptic_gregorian
    
Valid calendars 'standard', 'gregorian', 'proleptic_gregorian'
'noleap', '365_day', '360_day', 'julian', 'all_leap', '366_day'

Copyright (c) 2018, 2022 FESOM Development Team.

"""
from netCDF4 import Dataset, num2date
import sys

filename = sys.argv[1]
if len(sys.argv)>2:
    calendar = sys.argv[2]
else:
    calendar = '365_day'
    
f = Dataset(filename)

print(20*'*')
print('CALENDAR: ' + calendar)
print(20*'*')
      
for nstamp in range(f.variables['time'].shape[0]):
    estamp = num2date(f.variables['time'][:][nstamp], f.variables['time'].units, calendar)
    sstamp = int(f.variables['time'][nstamp])
    day = (sstamp//86400)+1
    seconds = sstamp%86400
    print(sstamp)
    print("{:5d} {:10d} {:10d}".format(seconds, day, estamp.year))
    print("{:5d} {:10d} {:10d}".format(86400, day, estamp.year))
    print(20*'*')

f.close()