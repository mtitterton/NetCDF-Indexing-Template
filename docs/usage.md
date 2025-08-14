Project Documentation: NetCDF File Indexing Template

Overview:
This program expands NETCDF files by indexing by perfoming the following operations:
1. Finding the surrounding files for each file and indexing data from each surrounding file into a larger NETCDF
2. Removing the repeats in latitude and longitude points (with their associated bathymetry data) and saving as a final NETCDF file

This program is designed to be dynamic by adjusting the requested indexing size and is suitable for applications in climate modeling

Features
Input handling: Accepts a list of NETCDF file names as input
String Manipulation: Maps each file name to its associated coordinate (if file name is related to a physical coordinate)
Data Expansion: Merges data from neighboring files based on if the surrounding coordinate is in the mapping of the file names
Duplicate Removal: Finds the overlaps between latitude and longitude for each indexing, and eliminates repeats in bathymetry data.
Latitude and Longitude Extrapolation: Takes each latitude and longitude data for the final NETCDF file and extrapolates an array so each coordinate is individually spaced.

Workflow
1. Input Loading


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
class File_List: 