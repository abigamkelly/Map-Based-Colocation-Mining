# Map-Based-Colocation-Mining
This github includes the code for our map-based regional colocation mining framework.  The following files are included:
* regional: this folder holds the map-based regional colocation framework code and files
    * Code: this folder holds the code for the map-based regional framework
        * functions.cpp: c++ functions used in the regional colcoation framework
        * regional_colocation.ipynb: python code that calls the c++ code in functions.cpp to perform the colocation mining
        * makefile: used to compile the c++ code
    * Data: this folder holds the real-world data sets
    * IntermediateData: this folder holds border region data passed from the python code to the c++ code
* regular: this folder holds the map-based colocation code and files used on synthetic testing
    * Code: this folder holds the map-based colocation mining code
        * functions.cpp: c++ code for the colocation mining
        * makefile: used to compile and run the code
    * Data: this folder holds the synthetic data sets

### How to Configure (for regional)
1. Open functions.cpp
2. On lines 85 and 86, set the distance threshold and prevalence threshold variables
3. On line 509, set the input file path variable to the correct file path holding the csv files for the region of choice
4. Open regional_colocation.ipynb
5. In the 5th cell, set the distance threshold and prevalence threshold variables
6. In the 5th cell, set the shapefile path and directory path variables to the correct path of the region of choice

### How to Compile and Run (for regional)
1. Change your current directory to the Code directory
2. Run the following command in the terimal to compile the c++ code: **make**
3. Open regional_colocation.ipynb and run all the cells

### How to Configure (for regular)
1. Open functions.cpp 
2. On lines 56 and 60, set the distance threshold and prevalence threshold variables
3. On line 466, set the input file variable to the correct file path of the synthetic data set of choice

### How to Compile and Run (for regular)
1. Change your current directory to the Code directory
2. Run the following command in the terminal to compile the code: **make**
3. Run the following command in the terminal to run the code: **make run**

