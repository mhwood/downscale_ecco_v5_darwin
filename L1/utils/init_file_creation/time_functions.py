
import os
import numpy as np
from datetime import datetime, timedelta

def dv_file_name_iter_to_iter_midpoints(file_name_iter, iters_per_output, iters_per_file):
    file_endpoints = np.arange(file_name_iter-1,file_name_iter-1+(iters_per_file+1)*iters_per_output,iters_per_output)
    file_midpoints = file_endpoints[:-1] + np.diff(file_endpoints)/2
    return(file_midpoints)

def iters_to_elapsed_seconds(start_seconds, seconds_per_iter, iter_numbers):
    elapsed_seconds = start_seconds + iter_numbers*seconds_per_iter
    return(elapsed_seconds)

def elapsed_seconds_to_iters(start_seconds, seconds_per_iter, elapsed_seconds):
    iter_numbers = int((elapsed_seconds-start_seconds)/seconds_per_iter)
    return(iter_numbers)

def elapsed_seconds_to_YMDHMS_str(start_YMD,start_HMS,elapsed_seconds):
    year = int(start_YMD[:4])
    month = int(start_YMD[4:6])
    day = int(start_YMD[6:8])
    hour = int(start_HMS[:2])
    minute = int(start_HMS[2:4])
    second = int(start_HMS[4:6])
    ref_datetime = datetime(year,month,day,hour,minute,second)
    elapsed_datetimes = []
    for seconds in elapsed_seconds.tolist():
        deltat = timedelta(seconds = seconds)
        elapsed_datetimes.append(ref_datetime+deltat)
    elapsed_datetime_strs = []
    for dt in elapsed_datetimes:
        YMD_Str = '{:02d}'.format(dt.year)+'{:02d}'.format(dt.month)+'{:02d}'.format(dt.day)
        HMS_Str = '{:02d}'.format(dt.hour) + '{:02d}'.format(dt.minute) + '{:02d}'.format(dt.second)
        elapsed_datetime_strs.append(YMD_Str+'-'+HMS_Str)
    return(elapsed_datetime_strs)


