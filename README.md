# Map-Based-Colocation-Mining
This is the read me for the regional colocation framework code.
There are three files needed to run the framework: distancce_threshold.ipynb, regional_colocation.ipynb, and c_functions.cpp.

There is a data set included: data.  This folder includes 2 sub-regions with a shapefile.  Currently, the code is set to run with this data set.

First, open the distance_threshold.ipynb.  The directory of the data files needs to be set in the second cell called 'directory'.  Once that variable has been set, run all of the cells.  The distance threshold will be produced and saved in the required_files directory as distance_threshold_parameter.txt.

Second, open a terminal in the current directory and type the following command: g++ -O3 -shared -o c_functions.so -fPIC c_functions.cpp
This compiles the cpp code and saves it as a shared library, so the other code has access to it.

Third, open the regional_colocation.ipynb file.  Adjust the file paths and prevalence_threshold as necessary in cell 2.  Run all the cells.  The prevalent patterns will be printed in the .ipynb file as well as saved in .txt files in the current directory.  Note: when running another data set, move the produced pattern files to another directory or delete them; they will be overwritten if left in the current directory.
