#
#MODULE ephem_analysis.py.py
#
#*****************************************************************************
"""
    **PURPOSE** --
    This tool that analyzes an ephemeris for HST.

    **DEVELOPER** --
    Gary Bower

    **MODIFICATION HISTORY** --

    o initial implementation - gab 10/14/25
    
"""
    
#*****************************************************************************

#import alignment_util
import ephem_util
import math
import numpy as np
import os
import re
import spss_sys_util
import spst_getopt
import sys
#import stpydb
#import su_gaps
#import su_util
import time_util

#from astropy.coordinates import SkyCoord
#from astropy.coordinates import Angle
#from astroquery.irsa import Irsa
#import astropy.units as u


__version__ = "25.10.14"

# Define constants.
HEADER_STRING = 'HEADER'
DIVIDER_STRING = '---'
DATA_HEADER = 'time                r           theta        phi'
UNITS = 'degrees'
THRESHOLD = 0.2


def run(*args):
    """
    **USAGE** --
        python ephem_analysis.py [ascii_ephem_file]
       [-output=<output file>] 
             <output file> is any valid filename.  If not given, output is sent
                 to the screen.
             
             Outputs:
      ASCII file with e.g., xxx.


    """
        
    # Parse the command line and determine cclist time window.
    allowed_options = ['output=']
    options, parms = spst_getopt.spst_getopt(args, allowed_options)
    print ('options: ', options, 'parms: ', parms)

    if '-output' in options:
        output_file = options['-output']
            
    # Define path to $SPSSPEF.
    spsspef = spss_sys_util.get_environ_variable("SPSSPEF")[0]
    if len(parms) == 1:
        ephem_asc_file = os.path.join(spsspef, parms[0].strip())
    else:
        files = spss_sys_util.eglob(os.path.join(spsspef, '*.asc'))
        ephem_asc_file = files[-1]
        
    header_dict = ephem_util.ascii_ephem(ephem_asc_file).get_header_dict()
    start_time = time_util.spss_time(header_dict['START_TIME'])
        
    # Print the ephemeris header.
    if '-output' in options:
        with open (output_file, 'w') as g:
            print (f'{HEADER_STRING}', file=g)
            print (f'{header_dict}', file=g)
            print (f'{DIVIDER_STRING*35}', file=g)
    else:    
        print (HEADER_STRING)
        print (header_dict)
        print (DIVIDER_STRING*35)
        
    # Read the file.
    data_all = open(ephem_asc_file, 'r').readlines()
    
    # Seperate out the 'top', i.e., the data before HEADER_STRING.
    top, header, data_dict_list = get_data_dict_list(data_all, start_time)
    
    data_header = list(data_dict_list[0].keys())

    if '-output' in options:
        with open (output_file, 'a+') as g:
            write_report (g, data_dict_list, header_dict)
    else:
        write_report (sys.stdout, data_dict_list, header_dict)  
   
    return spss_sys_util.SUCCESS
        
# end def run

def get_data_dict_list(data, start_time):
    """Seperate out the four of the seven data values out of an ASCII ephemeris."""
    
    for i in range(len(data)):

        if data[i].find(HEADER_STRING) != -1:
            top_index = i
            break
            
    top = data[0:top_index]
    
    for i in range(top_index+1, len(data)):
        if data[i].find(DIVIDER_STRING) != -1:
            end_header_index = i
            break
            
    # Separate out header            
    header = data[top_index+1:end_header_index-1]

    # Separate out the data.

    posn_dict_list = []
    posn_string = data[end_header_index+4:len(data)][0] 
    # posn_string split
    posn_list = posn_string.split()   
    for i in range(0, len(posn_list)-6,7):
        # Intialize dictionary for recording data.
        d = {}
        # Convert delta time in seconds to absolute time.
        d['t'] = start_time + time_util.delta_time(float(posn_list[i]))
        # Get EME2000 (Cartesian) coordinates.
        x = float(posn_list[i+1])
        y = float(posn_list[i+2])
        z = float(posn_list[i+3])
        # Convert to spherical coordinates (Physics standard). Then update dictionary.
        conv = eme2000_to_spherical (x, y, z, units=UNITS)
        d.update(conv)
        '''
        # Comment out differentials.
        dx = float(posn_list[i+4])
        dy = float(posn_list[i+5])
        dz = float(posn_list[i+6])
        '''
        posn_dict_list.append(d)
                
    return top, header, posn_dict_list
    
#end def get_data_dict_list

def eme2000_to_spherical (xVal, yVal, zVal, units='degrees'):
    """
    **DEVELOPER** --
    David Bower
    **PURPOSE**
    This tool is used to convert EME2000 (three-dimensional cartesian) to spherical coordinates 
    (physics standard). 

    **MODIFICATION HISTORY**
    - Initial implementation - DAB 10/4/25
    - Added radians to degrees conversion - DAB 10/5/25
    - Allowed user to pass in floats - DAB 10/6/25
    """

    # Perform the translation 
        # Calculate the radical distance 
        # Calculate the polar angle (theta) 
        # Calculate the azimuthal angle (phi)

    # Radical distance 
    radVal = getRad(xVal, yVal, zVal)

    # Polar Angle - just arccos(z/r)
    polAng = math.acos(zVal / radVal)

    # Azimuthal angle - just arctan(y/x)
    # Correcting this: 
    #azAng = math.atan(yVal / xVal)
 
    # Compute Azimuthal Angle. 
    # Use distance projected onto the xy-plane.
    azAng = np.sign(yVal) * math.acos(zVal/getRad(xVal, yVal, 0.))

    # If units=degrees, convert from radians to degrees
    if units == 'degrees': 
        polAng = math.degrees(polAng)
        azAng = math.degrees(azAng)
    
    data_dict = {'r': radVal, 'theta': polAng, 'phi': azAng}

    # Return the data to the user
    return data_dict

# Function for calculating radical distance
def getRad(xVal, yVal, zVal):

    # Initialize r-value and cumulative value 
    rVal = 0
    cumVal = 0

    """
    This segment of code places the x,y,z values in a list, loops over them, takes each value to the power of 2,
    then adds them together. The final cumulative number is stored in a variable previously initialized. 
    """

    vals = [xVal, yVal, zVal] 
    for i in range(len(vals)):
        vals[i] = math.pow(vals[i], 2)
        cumVal += vals[i]

    # Square root the cumulative value for the final answer
    rVal = math.sqrt(cumVal)

    return rVal
    
def write_report (device, data_dict_list, header_dict):
    """Write to the output."""
    print (f'{DATA_HEADER}', file=device)
    print (f'                   km          {UNITS}      {UNITS}', file=device)
    print (f'{DIVIDER_STRING*35}', file=device)
    for d in data_dict_list:
        r_error = abs(1. - d['r']/header_dict['KEPLERIAN_SEMIM_AXIS'])
        if r_error <= THRESHOLD:
            print (f"{d['t']}  {d['r']:>8.2f}  {d['theta']:>8.2f}  {d['phi']:>8.2f}", file=device)
        else:
            pass
            
    # Print trend.
    first = data_dict_list[0]
    last = data_dict_list[-2]
    print (f'{DIVIDER_STRING*35}', file=device)
    print ('   TREND', file=device)
    print (f'{first}', file=device)
    print (f'{last}', file=device)

    return
    
# end def write_report


##### MAIN #####
################
        
if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        print(sys.argv)
        run(*tuple(sys.argv[1:]))
    else:
        run()