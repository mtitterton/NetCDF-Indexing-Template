Project Documentation: NetCDF File Indexing Template

Overview:  
This program expands NETCDF files by indexing by perfoming the following operations:
1. Finding the surrounding files for each file and indexing data from each surrounding file into a larger NETCDF
2. Removing the repeats in latitude, longitude, and bathymetry data points and saving as a final NETCDF file

This program is designed to be dynamic by adjusting the requested indexing size and is suitable for applications in climate modeling.

Features  
Input handling: Accepts a list of NETCDF file names as input  
String Manipulation: Maps each file name to its associated coordinate (if file name is related to a physical coordinate)  
Data Expansion: Merges data from neighboring files based on if the surrounding coordinate is in the mapping of the file names  
Duplicate Removal: Finds the overlaps between latitude and longitude for each indexing, and eliminates repeats in bathymetry data.  
Latitude and Longitude Extrapolation: Takes each latitude and longitude data for the final NETCDF file and extrapolates an array so each coordinate is individually spaced.  

Workflow  
The program follows these steps:
1. Input loading: 
- Reads in all NetCDF files and converts string names to coordinates
- Identifies surrounding coordinates and creates lists of existing coordinates in file names

2. Data Expansion:
- Extracts data from surrounding files and adds them to larger NetCDf file for each given file name

3. Duplicate Removal:
- Finds repeats in data points and removes them with indexing while creating a more cohesive NetCDF file

4. Output Generation:
- Saves the processed dataset as a final new NetCDF file


Usage Instructions
Prerequisites:
- Python 3.11.9
- Required libraries: 
    - xarray for opening and reading datasets
    - numpy for numerical operations
    - netCDF4 (Dataset) for creating new NetCDF files
    - os for string manipulation

Install dependencies using pip:

pip install netCDF numpy os xarray

Code Structure
Classes
1. File_list:
- Parent class that is used to adjust the desired indexing size from surrounding files (val)
- Establish the size of each file within this class (Nx, Ny)
- Pass through the file list for the rest of the program (self.file_list)

2. string_parse_function:
- Takes all file names, converts them to physical meanings in the form of coordinates, finds the surrounding coordinates for each coordinate, and creates new file by indexing from surrounding files.
- Method "create_netcdf" is for mapping the file name to a new file of the desired size (its original size plus val length on either end), making sure to include all variables and dimensions.
- Method "cut_string" creates the dictionaries for mapping the original file names to their physical coordinate representations.
- Method "check_boxes" creates hypothetical surrounding coordinates for each file and checks if they exist in the coordinate list, then creates a list for each coordinate of surrounding points (later used in a dataframe with pandas)
- Method "add_data" takes the new NetCDF file associated with the original file name and indexes lat, lon, and bathymetry (Band1) data from the existing surrounding files.

3. overlapping_data:
- The lat and lon data slightly overlap between each surrounding file, so this class finds where those repeats occur and then re-index accordingly.
- Method "eliminate_repeats" creates a final NetCDF for each file, evaluates the latitude and longitude values at each boundary (which occur at the original border of the file as this is where the data is added on), find the repeats, then indexes the data accordingly to eliminate the repeats using the same process as in class string_parse_function.
- Method "interpolate_lat_lon" is used later on in method "fill_lat_lon", taking in an array (using lat or lon), finding the spacing between values, and then interpolating the values on each end so each point is evenly spaced.
- Optional: Method "fill_lat_lon" is used to adjust the latitude and longitude arrays from the final dataset, as indexing can sometimes leads to slightly mismatches in the spacing between values with rounding errors. This method doesn't have to be done, but it is helpful for more consistent NetCDF files.

Algorithm Details:
Data Expansion  
This program identifies surrounding files based on converting string names to coordinates for comparison


Output  
The output is a list of NetCDF files with:
- Expanded data from surrounding files
- No duplicate entries

Error Handling:
This program is heavily dependent on the fact that all NetCDF file sizes are the same. If they are not, errors will be thrown, and indexing should be done manually.
This program does include error handling for:
- If latitude and longitude repeats cannot be found to the rounding specifications within the "eliminate_repeats" method in the "overlapping_data" class, then program will return a message saying so. However, this will cause further issues later with indexing for repeats.
Files not of the same sizes or with issues in overlapping data should be handled manually uses similar techniques to the program

Future Improvements
- Create a more dynamic system so each file does not have to be the same size
- Create an alternate method for dealing with individual files that don't fit with the sizes and repeat definitions of the program
