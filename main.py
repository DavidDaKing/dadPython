# init main

# *********************************************************************************************
"""
    **PURPOSE**
    This tool is used to assit with three-dimensional cartesian to spherical coordinate system calculations. 

    **MODIFICATION HISTORY**
    - Initial implementation - DAB 10/4/25

"""
# *********************************************************************************************

import sys
import math

# This is the run function 
def run(*args):
    """
        ** USAGE** --
        python main.py {x-val} {y-val} {z-val}
    """
    if not args:
        print("spew out the usage and quit when no args given")

    # Parse the command line for the given three-dimensional cartesian coordinates
        # Acceptable format: python main.py x-val y-val z-val
        # May be positive or negative values (Not implemented yet)
    if len(sys.argv) != 4:
        print("Usage: python3 main.py {x-val} {y-val} {z-val}")
        sys.exit(1)

    # User given coordinate values 
    xVal = sys.argv[1]
    yVal = sys.argv[2]
    zVal = sys.argv[3]


    # Printed 
    print(f"x: {xVal}, y: {yVal}, z: {zVal}")

    # Perform the translation 
        # Calculate the radical distance 
        # Calculate the polar angle 
        # Calculate the azimuthal angle 

    # Radical distance 
    radVal = getRad(xVal, yVal, zVal)

    # Polar Angle - just arccos(z/r)
    # TypeError: unsupported operand type(s) for /: 'str' and 'float' - cast zVal as int
    polAng = math.acos(int(zVal) / radVal)

    # Azimuthal angle - just arctan(y/x) - cast both as ints
    azAng = math.atan(int(yVal) / int(xVal))

    # Radians to degree conversions 
    polDeg = polAng * (180 / math.pi)
    azDeg = azAng * (180 / math.pi)

    # Spew out the data to the user
    print(f"This is the radical distance: {radVal}")
    print(f"This is the polar angle: {polAng} radians | {polDeg} degrees")
    print(f"This is the azimuthal angle: {azAng} radians | {azDeg} degrees")

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
        # TypeError: must be real number, not string - to avoid: cast as int
        vals[i] = math.pow(int(vals[i]), 2)
        cumVal += vals[i]

    # Square root the cumulative value for the final answer
    rVal = math.sqrt(cumVal)

    return rVal


##### MAIN #####
################

if __name__ == '__main__':
    if len(sys.argv) > 1:
        #print(sys.argv)
        run(*tuple(sys.argv[1:]))
    else:
        run()
