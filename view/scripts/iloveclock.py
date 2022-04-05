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

Copyright (c) 2018, FESOM Development Team.

"""
from netCDF4 import Dataset, num2date
from datetime import timedelta
from datetime import datetime
import sys

filename = sys.argv[1]
if len(sys.argv)>2:
    calendar = sys.argv[2]
else:
    calendar = '365_day'
    
f = Dataset(filename)
# a = num2date(f.variables['time'][:], f.variables['time'].units, '365_day')

print(20*'*')
print('CALENDAR: ' + calendar)
print(20*'*')
      
for nstamp in range(f.variables['time'].shape[0]):
    sstamp = num2date(f.variables['time'][nstamp], f.variables['time'].units, calendar)
    delta = (60 - sstamp.minute)*60
    estamp = num2date(f.variables['time'][:][nstamp] + delta, f.variables['time'].units, calendar)
    seconds_in_day_s = sstamp.hour*3600+sstamp.minute*60
    seconds_in_day_e = estamp.hour*3600+estamp.minute*60

    if calendar in ['noleap', '365_day', '360_day', '366_day']:
        print(sstamp)
        print("{:5d} {:10d} {:10d}".format(seconds_in_day_s, sstamp.dayofyr, sstamp.year))
        print("{:5d} {:10d} {:10d}".format(seconds_in_day_e, estamp.dayofyr, estamp.year))
        print(20*'*')
    else:
        print(sstamp)
        print("{:5d} {:10d} {:10d}".format(seconds_in_day_s, sstamp.timetuple().tm_yday, sstamp.year))
        print("{:5d} {:10d} {:10d}".format(seconds_in_day_e, estamp.timetuple().tm_yday, estamp.year))
        print(20*'*')

print(20*'*')
print('CALENDAR: ' + calendar)
print(20*'*')

f.close()