# Map-Based-Colocation-Mining
This github includes the code for our map-based regional colocation mining framework.  The following files are included:
 * Code: this folder holds the code for the map-based regional framework
     * functions.cpp: c++ functions used in the regional colcoation framework
     * regional_colocation.ipynb: python code that calls the c++ code in functions.cpp to perform the colocation mining
     * makefile: used to compile the c++ code
 * Data: this folder holds the real-world data sets
 * IntermediateData: this folder holds border region data passed from the python code to the c++ code


### How to Configure
1. Open functions.cpp
2. On lines 109 and 110, set the distance threshold and prevalence threshold variables
3. On line 525, set the input file path variable to the correct file path holding the csv files for the region of choice
4. Open regional_colocation.ipynb
5. In the 5th cell, set the distance threshold and prevalence threshold variables
6. In the 5th cell, set the shapefile path and directory path variables to the correct path of the region of choice
7. Open the makefile and set the BOOST_PATH variable on line 4 to the appropiate boost library path

### How to Compile and Run
1. Change your current directory to the Code directory
2. Run the following command in the terimal to compile the c++ code: **make**
3. Open regional_colocation.ipynb and run all the cells

