import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
from netCDF4 import Dataset
import os
import numpy as np

''' Steps for program:
1. Find surrounding files
2. Concatenate larger file from surrounding files
3. Use indexing to eliminate overlapping data
'''

class File_List:

    def __init__(self, file_list = None):
        self.file_list = file_list
        self.Nx = 8112
        self.Ny = 8112
        self.val = 1622
    

class string_parse_function(File_List):

    def __init__(self, file_list):
        super().__init__(file_list)
        self.row_val = 0
        self.col_val = 0
        self.n_columns = 8
        self.n_rows = len(self.file_list)
        self.surrounding_boxes_list = [[0 for _ in range(self.n_columns)] for _ in range(self.n_rows)]
        self.surrounding_file_names = [[0 for _ in range(self.n_columns)] for _ in range(self.n_rows)]
        self.coordinate_list = []
        self.checking_list = [[0 for _ in range(self.n_columns)] for _ in range(self.n_rows)]
        self.datasets = {}
        self.new_file_names = []


    def create_netcdf(self):
        # method to create a new empty netcdf file for each given file with a larger size
        # value to extend each file by
        for i, item in enumerate(self.file_list):
            # lat and lon lengths
            new_length = self.Nx + (2*self.val)
            # create netcdf file with variables: lat, lon, Band1 (bathymetry) and dimensions: lat, lon
            base, ext = os.path.splitext(item)
            new_name = base + "_new" + ext
            ncfile = Dataset(new_name, mode='w', format='NETCDF4_CLASSIC')
            ncfile.createDimension('lat', new_length)
            ncfile.createDimension('lon', new_length)
            ncfile.createVariable('Band1',np.float32,('lat','lon'))
            ncfile.createVariable('lon', np.float32, ('lon'))
            ncfile.createVariable('lat', np.float32, ('lat'))

            # map the dataset to its original item name
            self.datasets[item] = new_name

            # take lat, lon, bathy data from original set and transfer over to new dataset
            original_dataset = xr.open_dataset(item)
            band_var = original_dataset.variables['Band1']
            band_values = band_var.values
            lon_variable = original_dataset.variables['lon']
            lon_values = lon_variable.values
            lat_variable = original_dataset.variables['lat']
            lat_values = lat_variable.values
            original_dataset.close()
            
            # initialize every value in new datset as np.nan
            ncfile.variables['Band1'][:, :] = np.nan
            ncfile.variables['lon'][:] = np.nan
            ncfile.variables['lat'][:] = np.nan
            
            # fill in original dataset within larger dataset
            ncfile.variables['Band1'][self.val:new_length-self.val, self.val:new_length-self.val] = band_values
            ncfile.variables['lon'][self.val:new_length-self.val] = lon_values
            ncfile.variables['lat'][self.val:new_length-self.val] = lat_values

            self.new_file_names.append(new_name)

            print("iniitalized: " + str(i))


    def cut_string(self):
    # method to manipulate all file names into comparable data points
    # NOTE: negatives for longitude aren't converted, the "0" in the file name indicates negative but is ignored for comparisons

    # dictionaries for going from data point to filename, and reverse
        self.forward = {}
        self.reverse = {}

        # loop through each file name to split at the underscores
        for item in self.file_list:
            # split into latitude and longitude parts for consistency
            longitude_name = item.split('_')[1][1:] 
            latitude_name = item.split('_')[2][2:] 

            # create new instances for each item in the list
            longitude = "" 
            latitude = ""

            # create data points by cycling through each piece of the string, specifically eliminating the x in each point here
            longitude = "".join(["." if char == "x" else char for char in longitude_name])
            latitude = "".join(["." if char == "x" else char for char in latitude_name])
            
            # turn strings into floats, round to 4 places for consistency in comparisons later
            longitude = round(float(longitude), 4)
            latitude = round(float(latitude),4)

            # add to final list of data points
            self.coordinate_list.append([longitude, latitude])

            # append ponts and filenames to dictionaries
            self.forward[(longitude, latitude)] = item
            self.reverse[item] = [longitude,latitude]


    def check_boxes(self):
        # method to check if surrounding files are present within the file list
        # creates 8 hypothetical coordinates (n, e, s, w, nw, ne, sw, se), checks if they are present in the list, and if so, appends them to the surrounding_files item
        # NOTE: item[0] = longitude, item[1] = latitude

        for i, item in enumerate(self.coordinate_list):
            # up 0.25
            # n_coord
            self.checking_list[i][0] =  [round(float(item[0] + 0.25), 4), round(float(item[1]), 4)]
            # right 0.25
            # e_coord
            self.checking_list[i][1] = [round(float(item[0]), 4), round(float(item[1] - 0.25), 4)]
            # down 0.25
            # s_coord 
            self.checking_list[i][2] = [round(float(item[0] - 0.25), 4), round(float(item[1]), 4)]
            # left 0.25
            # w_coord 
            self.checking_list[i][3] = [round(float(item[0]), 4), round(float(item[1] + 0.25), 4)]
            # up and left 0.25
            # nw_coord 
            self.checking_list[i][4] = [round(float(item[0] + 0.25), 4), round(float(item[1] + 0.25), 4)]
            # up and right 0.25
            # ne_coord 
            self.checking_list[i][5] = [round(float(item[0] + 0.25), 4), round(float(item[1] - 0.25), 4)]
            # down and left 0.25
            # sw_coord 
            self.checking_list[i][6] = [round(float(item[0]- 0.25), 4), round(float(item[1] + 0.25), 4)]
            # down and right 0.25
            # se_coord 
            self.checking_list[i][7] = [round(float(item[0]- 0.25), 4), round(float(item[1] - 0.25), 4)]

           # loop through each item to check if surrounding coordinates are present
            for index, item2 in enumerate(self.checking_list[i]):
                if item2 in self.coordinate_list:
                    # if present, add to surrounding files for that file
                    coord_tuple = tuple(item2)
                    file_to_add = self.forward[coord_tuple]
                    self.surrounding_file_names[i][index] = file_to_add
                else:
                    # if not present, fill place with 0
                    self.surrounding_file_names[i][index] = 0


    def add_data(self):
        # method for appending data onto the original file
        # takes the bigger file corresponding to the original file and indexes from the surrounding files

        for index1, item1 in enumerate(self.surrounding_file_names):
            for index, item in enumerate(self.surrounding_file_names[index1]):
                appending_to = self.file_list[index]
                self.val = 1622
                self.Nx = 8112
                self.Ny = 8112
                length = self.Nx + 2*self.val

                print(item)
                if item != 0:
                    if index == 0:
                        # north
                        # grab data from original dataset
                        print(item)
                        ds = xr.open_dataset(item)

                        # assign variables, multiple steps to ensure they're converted to numpy arrays
                        var = ds.variables['Band1']
                        new_var = var.values
                        lon_var = ds.variables['lon']
                        lon_values = lon_var.values
                        lat_var = ds.variables['lat']
                        lat_values = lat_var.values
                        slice_taken_lat = lat_values[0:self.val]
                        slice_taken = new_var[0:self.val, :]
                        ds.close()

                        # grab bigger file
                        file_of_interest = self.datasets[appending_to]
                        with Dataset(file_of_interest, 'r+') as ds:
                            # assign slices from original data slice to new bigger file, making sure to keep the slice sizes the same
                            ds.variables['Band1'][length-self.val:length, self.val:length-self.val] = slice_taken
                            ds.variables['lat'][length-self.val: length] = slice_taken_lat

                    elif index == 1:
                        # east
                        ds = xr.open_dataset(item)
                        var = ds.variables['Band1']
                        new_var = var.values
                        slice_taken = new_var[:, 0:self.val]
                        lon_var = ds.variables['lon']
                        lon_values = lon_var.values
                        lat_var = ds.variables['lat']
                        lat_values = lat_var.values
                        slice_taken_lon = lon_values[0:self.val]
                        ds.close()
                        file_of_interest = self.datasets[appending_to]
                        with Dataset(file_of_interest, 'r+') as ds:
                            ds.variables['Band1'][self.val:length-self.val, length-self.val:length] = slice_taken
                            ds.variables['lon'][length-self.val:length] = slice_taken_lon

                    elif index == 2:
                        # south
                        ds = xr.open_dataset(item)
                        var = ds.variables['Band1']
                        new_var = var.values
                        lon_var = ds.variables['lon']
                        lon_values = lon_var.values
                        lat_var = ds.variables['lat']
                        lat_values = lat_var.values
                        slice_taken = new_var[self.Nx-self.val-1:self.Nx-1, :]
                        slice_taken_lat = lat_values[self.Nx-self.val-1:self.Nx-1]
                        ds.close()
                        file_of_interest = self.datasets[appending_to]
                        with Dataset(file_of_interest, 'r+') as ds:
                            ds.variables['Band1'][0:self.val, self.val:length-self.val] = slice_taken
                            ds.variables['lat'][0:self.val] = slice_taken_lat

                    elif index == 3:
                        # west
                        ds = xr.open_dataset(item)
                        var = ds.variables['Band1']
                        new_var = var.values
                        lon_var = ds.variables['lon']
                        lon_values = lon_var.values
                        lat_var = ds.variables['lat']
                        lat_values = lat_var.values
                        slice_taken = new_var[:, self.Nx-self.val-1:self.Nx-1]
                        slice_taken_lon = lon_values[self.Nx-self.val-1:self.Nx-1]
                        ds.close()
                        file_of_interest = self.datasets[appending_to]
                        with Dataset(file_of_interest, 'r+') as ds:
                            ds.variables['Band1'][self.val:length-self.val, 0:self.val] = slice_taken
                            ds.variables['lon'][0:self.val] = slice_taken_lon

                    elif index == 4:
                        # NW
                        ds = xr.open_dataset(item)
                        var = ds.variables['Band1']
                        new_var = var.values
                        lon_var = ds.variables['lon']
                        lon_values = lon_var.values
                        lat_var = ds.variables['lat']
                        lat_values = lat_var.values
                        slice_taken_lon = lon_values[self.Ny-self.val-1: self.Nx-1]
                        slice_taken_lat = lat_values[0:self.val]
                        slice_taken = new_var[0:self.val, self.Nx-self.val-1:self.Nx-1]
                        ds.close()
                        file_of_interest = self.datasets[appending_to]
                        with Dataset(file_of_interest, 'r+') as ds:
                            ds.variables['Band1'][length-self.val:length, 0:self.val] = slice_taken
                            ds.variables['lon'][0:self.val] = slice_taken_lon
                            ds.variables['lat'][length-self.val:length] = slice_taken_lat


                    elif index == 5:
                        # NE
                        ds = xr.open_dataset(item)
                        var = ds.variables['Band1']
                        new_var = var.values
                        lon_var = ds.variables['lon']
                        lon_values = lon_var.values
                        lat_var = ds.variables['lat']
                        lat_values = lat_var.values
                        slice_taken_lon = lon_values[0:self.val]
                        slice_taken_lat = lat_values[0:self.val]
                        slice_taken = new_var[0:self.val, 0:self.val]
                        ds.close()
                        file_of_interest = self.datasets[appending_to]
                        with Dataset(file_of_interest, 'r+') as ds:
                            ds.variables['Band1'][length-self.val:length, length-self.val:length] = slice_taken
                            ds.variables['lon'][length-self.val:length] = slice_taken_lon
                            ds.variables['lat'][length-self.val:length] = slice_taken_lat

                    elif index == 6:
                        # SW
                        ds = xr.open_dataset(item)
                        var = ds.variables['Band1']
                        new_var = var.values
                        lon_var = ds.variables['lon']
                        lon_values = lon_var.values
                        lat_var = ds.variables['lat']
                        lat_values = lat_var.values
                        slice_taken = new_var[self.Nx-self.val-1:self.Nx-1, self.Nx-self.val-1:self.Nx-1]
                        slice_taken_lon = lon_values[self.Nx- self.val - 1: self.Nx -1]
                        slice_taken_lat = lat_values[self.Nx- self.val - 1: self.Nx -1]
                        ds.close()
                        file_of_interest = self.datasets[appending_to]
                        with Dataset(file_of_interest, 'r+') as ds:
                            ds.variables['Band1'][0:self.val, 0:self.val] = slice_taken
                            ds.variables['lon'][0:self.val] = slice_taken_lon
                            ds.variables['lat'][0:self.val] = slice_taken_lat

                    elif index == 7:
                        # SE
                        ds = xr.open_dataset(item)
                        var = ds.variables['Band1']
                        new_var = var.values
                        lon_var = ds.variables['lon']
                        lon_values = lon_var.values
                        lat_var = ds.variables['lat']
                        lat_values = lat_var.values
                        slice_taken = new_var[self.Nx-self.val-1:self.Nx-1, 0:self.val]
                        slice_taken_lon = lon_values[0:self.val]
                        slice_taken_lat = lat_values[self.Nx- self.val - 1: self.Nx -1]
                        ds.close()
                        file_of_interest = self.datasets[appending_to]
                        with Dataset(file_of_interest, 'r+') as ds:
                            ds.variables['Band1'][0:self.val, length-self.val:length] = slice_taken
                            ds.variables['lon'][length-self.val: length] = slice_taken_lon
                            ds.variables['lat'][0:self.val] = slice_taken_lat
                    else:
                        print("error")



class overlapping_data(File_List):
    def __init__(self, surrounding_files, new_file_names):
        super().__init__()
        self.surrounding_files = surrounding_files
        self.new_file_names = new_file_names

    def eliminate_repeats(self):
        for new_index, item_name in enumerate(self.new_file_names):
            # method for finding where lat and lon overlap, then doing same indexing technique from concatenating larger files with indexes of overlaps

            # initialize all lists as empty and variables as nans
            same_value_index_lat1 = np.nan
            same_value_index_lat2 = np.nan
            same_value_index_lon1 = np.nan
            same_value_index_lon2 = np.nan
            same_value_lat1 = np.nan
            same_value_lat2 = np.nan
            same_value_lon1 = np.nan
            same_value_lon2 = np.nan
            skip_vals_lat1 = np.nan
            skip_vals_lat2 = np.nan
            skip_vals_lon1 = np.nan
            skip_vals_lon2 = np.nan
            lat_first_slice = []
            lat_second_slice = []
            lon_first_slice = []
            lon_second_slice = []
            lat_middle_slice = []
            lon_middle_slice = []
            lat_final = []
            lon_final = []
            length = self.Nx + 2*self.val

            # length around end of original datafile and into new indexing for evaluating where repeat happens
            buffer = 20

            # open dataset and grab values
            evaluator = self.new_file_names[new_index]
            ds = xr.open_dataset(evaluator)
            lat_var = ds.variables['lat']
            lat_values = lat_var.values
            lon_var = ds.variables['lon']
            lon_values = lon_var.values
            first_bound_lat = lat_values[self.val-buffer:self.val+buffer]
            first_bound_lon = lon_values[self.val-buffer:self.val+buffer]
            second_bound_lat = lat_values[length-self.val-buffer:length-self.val+buffer]
            second_bound_lon = lon_values[length-self.val-buffer:length-self.val+buffer]

            # for loops to find where the lat and lon overlap by finding where difference is negative or positive (lat versus lon)
            for i, value in enumerate(first_bound_lat[:-1]):
                difference = first_bound_lat[i+1] - value
                if difference < 0:
                    same_value_index_lat1 = i
                    same_value_lat1 = first_bound_lat[i+1]
                    break
                
            for i, value in enumerate(second_bound_lat[:-1]):
                difference = second_bound_lat[i+1] - value
                if difference < 0:
                    same_value_index_lat2 = i 
                    same_value_lat2 = second_bound_lat[i+1]
                    break

            for i, value in enumerate(first_bound_lon[:-1]):
                difference = first_bound_lon[i+1] - value
                if difference < 0:
                    same_value_index_lon1 = i 
                    same_value_lon1 = first_bound_lon[i+1]
                    break

            for i, value in enumerate(second_bound_lon[:-1]):
                difference = second_bound_lon[i+1] - value
                if difference < 0:
                    same_value_index_lon2 = i 
                    same_value_lon2 = second_bound_lon[i+1]
                    break

            # find the number of values to skip for indexing
            # if the skip_vals variable still return as np.nan, that means something went wrong in the comparison, do it manually for that set
            if not np.isnan(same_value_index_lat1):
                for i in range(len(first_bound_lat)):
                    if first_bound_lat[i] == same_value_lat1:
                        skip_vals_lat1 = i
                        lat_first_slice = lat_values[:self.val-buffer+skip_vals_lat1]
                        break
                if np.isnan(skip_vals_lat1):
                    print("values not close enough for comparison")

            if not np.isnan(same_value_index_lat2):
                for i in range(same_value_index_lat2, -1, -1):
                    if np.isclose(second_bound_lat[i], same_value_lat2, rtol=0, atol=1e-5):
                        skip_vals_lat2 = i
                        lat_second_slice = lat_values[length - self.val - buffer + same_value_index_lat2 + skip_vals_lat2 + 3:]
                        break
                if np.isnan(skip_vals_lat2):
                    print("values not close enough for comparison")

            if not np.isnan(same_value_index_lon1):
                for i in range(len(first_bound_lon)):
                    if first_bound_lon[i] == same_value_lon1:
                        skip_vals_lon1 = i
                        lon_first_slice = lon_values[:self.val-buffer+skip_vals_lon1]
                        break
                if np.isnan(skip_vals_lon1):
                    print("values not close enough for comparison")


            if not np.isnan(same_value_index_lon2):
                for i in range(same_value_index_lon2,-1,-1):
                    if np.isclose(second_bound_lon[i], same_value_lon2, rtol=0, atol=1e-5):
                        skip_vals_lon2 = i
                        lon_second_slice = lon_values[length-self.val-buffer+same_value_index_lon2+skip_vals_lon2+3:]
                        break
                if np.isnan(skip_vals_lon2):
                    print("values not close enough for comparison")

            ds.close()    

            # open final dataset and assign values to it
            ds2 = xr.open_dataset(item_name)
            lat_original = ds2.variables['lat']
            lat_original_values = lat_original.values
            lon_original = ds2.variables['lon']
            lon_original_values = lon_original.values
            band_original = ds2.variables['Band1']
            band_original_values = band_original.values
            lat_middle_slice = lat_original_values
            lon_middle_slice = lon_original_values
            middle_slice_band = band_original_values
            var = ds2.variables['Band1']
            new_var = var.values
            slice_taken = new_var[::]
            ds2.close()

            # concatenate lat left, right, and middle
            lat_final = np.concatenate((lat_first_slice, lat_middle_slice, lat_second_slice))
            # concatenate lon left, right, and middle
            lon_final = np.concatenate((lon_first_slice, lon_middle_slice, lon_second_slice))


            # find difference between original length and new length for lat/lon
            bottom_lat_diff = self.val - len(lat_first_slice)
            top_lat_diff = self.val - len(lat_second_slice)
            left_lon_diff = self.val - len(lon_first_slice)
            right_lon_diff = self.val - len(lon_second_slice)

            # Save to new NetCDF
            base, ext = os.path.splitext(item_name)
            new_name = base + "_band" + ext
            ncfile1 = Dataset(new_name, mode='w', format='NETCDF4_CLASSIC')
            ncfile1.createDimension('lat', 11356)
            ncfile1.createDimension('lon', 11356)
            ncfile1.createVariable('lat', np.float32, ('lat'))
            ncfile1.createVariable('lon', np.float32, ('lon'))
            ncfile1.createVariable('Band1', np.float32, ('lat', 'lon'))
            ncfile1.variables['Band1'][:,:] = np.nan
            ncfile1.variables['Band1'][self.val+1:length-self.val+1, self.val+1:length-self.val+1] = slice_taken
            ncfile1.variables['lat'][:] = np.nan
            ncfile1.variables['lon'][:] = np.nan
            ncfile1.variables['lat'][bottom_lat_diff:length-top_lat_diff] = lat_final
            ncfile1.variables['lon'][left_lon_diff:length-right_lon_diff] = lon_final
            ncfile1.close()

            # use same indexing technique as before to take data from surrounding datasets and use the skip lat and lon values
            for index_val, item2  in enumerate(self.surrounding_file_names[new_index]):

                if item2 != 0:
                    length = self.Nx + 2*self.val
                    if index_val == 0:
                        ds = xr.open_dataset(item2)
                        var = ds.variables['Band1']
                        new_var = var.values
                        slice_taken = new_var[same_value_index_lat2+skip_vals_lat2-20+3:self.val, :]
                        gap = same_value_index_lat2+skip_vals_lat2-20+3
                        ds.close()
                        with Dataset(new_name, 'r+') as ds:
                            ds.variables['Band1'][length-self.val: length-gap, self.val:length-self.val] = slice_taken
                    # case north

                    elif index_val == 1:
                        ds = xr.open_dataset(item2)
                        var = ds.variables['Band1']
                        new_var = var.values
                        slice_taken = new_var[:, same_value_index_lon2+skip_vals_lon2-20+3:self.val]
                        gap  = same_value_index_lon2+skip_vals_lon2-20+3
                        ds.close()
                        with Dataset(new_name, 'r+') as ds:
                            ds.variables['Band1'][self.val:length-self.val, length-self.val: length-gap] = slice_taken
                    # case east

                    elif index_val == 2:
                        ds = xr.open_dataset(item2)
                        var = ds.variables['Band1']
                        new_var = var.values
                        slice_taken = new_var[self.Nx-self.val:self.Nx-self.val+self.val-buffer+skip_vals_lat1, :]
                        gap = self.val - (self.val - buffer+skip_vals_lat1)
                        ds.close()
                        with Dataset(new_name, 'r+') as ds:
                            ds.variables['Band1'][gap:self.val, self.val:length-self.val] = slice_taken
                        # case south

                    elif index_val == 3:
                        ds = xr.open_dataset(item2)
                        var = ds.variables['Band1']
                        new_var = var.values
                        slice_taken = new_var[:, self.Nx-self.val:self.Nx-self.val+self.val-buffer+skip_vals_lon1]
                        gap = self.val - (self.val - buffer+skip_vals_lon1)
                        ds.close()
                        with Dataset(new_name, 'r+') as ds:
                            ds.variables['Band1'][self.val:length-self.val, gap:self.val] = slice_taken
                        # case west

                    elif index_val == 4:
                        ds = xr.open_dataset(item2)
                        var = ds.variables['Band1']
                        new_var = var.values
                        slice_taken = new_var[same_value_index_lat2+skip_vals_lat2-20+3:self.val, self.Nx-self.val:self.Nx-self.val+self.val-buffer+skip_vals_lon1]
                        gap1 = same_value_index_lat2+skip_vals_lat2-20+3
                        gap2 = self.val - (self.val - buffer+skip_vals_lon1)
                        ds.close()
                        with Dataset(new_name, 'r+') as ds:
                            ds.variables['Band1'][length-self.val: length-gap1, gap2:self.val] = slice_taken
                        # case NW 

                    elif index_val == 5:
                        ds = xr.open_dataset(item2)
                        var = ds.variables['Band1']
                        new_var = var.values
                        slice_taken = new_var[same_value_index_lat2+skip_vals_lat2-20+3:self.val, same_value_index_lon2+skip_vals_lon2-20+3:val]
                        gap1 = same_value_index_lat2+skip_vals_lat2-20+3
                        gap2  = same_value_index_lon2+skip_vals_lon2-20+3
                        ds.close()
                        with Dataset(new_name, 'r+') as ds:
                            ds.variables['Band1'][length-self.val: length-gap1, length-self.val: length-gap2] = slice_taken

                        # case NE
                    elif index_val == 6:
                        ds = xr.open_dataset(item2)
                        var = ds.variables['Band1']
                        new_var = var.values
                        slice_taken = new_var[Nx-self.val:self.Nx-self.val+self.val-buffer+skip_vals_lat1, self.Nx-self.val:self.Nx-self.val+self.val-buffer+skip_vals_lon1]
                        gap1 = self.val - (self.valval - buffer+skip_vals_lat1)
                        gap2 = self.val - (self.val - buffer+skip_vals_lon1)
                        ds.close()
                        with Dataset(new_name, 'r+') as ds:
                            ds.variables['Band1'][gap1:self.val, gap2:self.val] = slice_taken

                        # case SW

                    elif index_val == 7:
                        ds = xr.open_dataset(item2)
                        var = ds.variables['Band1']
                        new_var = var.values
                        slice_taken = new_var[self.Nx-self.val:self.Nx-self.val+self.val-buffer+skip_vals_lat1, same_value_index_lon2+skip_vals_lon2-20+3:self.val]
                        gap1 = self.val - (self.val - buffer+skip_vals_lat1)
                        gap2 = same_value_index_lon2+skip_vals_lon2-20+3
                        ds.close()
                        with Dataset(new_name, 'r+') as ds:
                            ds.variables['Band1'][gap1:self.val, length-self.val: length-gap2] = slice_taken
                        # case SE

    def interpolate_lat_lon(self, arr):
        arr = np.array(arr, dtype=float)
        n = len(arr)
        # Find indices of non-NaN values
        valid_idx = np.where(~np.isnan(arr))[0]
        if len(valid_idx) < 2:
            raise ValueError("Need at least two valid values for extrapolation.")
        result = arr.copy()
        # Extrapolate at the start
        start = valid_idx[0]
        if start > 0:
            slope = (arr[valid_idx[1]] - arr[valid_idx[0]]) / (valid_idx[1] - valid_idx[0])
            for i in range(start - 1, -1, -1):
                result[i] = result[i + 1] - slope
        # Extrapolate at the end
        end = valid_idx[-1]
        if end < n - 1:
            slope = (arr[valid_idx[-1]] - arr[valid_idx[-2]]) / (valid_idx[-1] - valid_idx[-2])
            for i in range(end + 1, n):
                result[i] = result[i - 1] + slope

        return result
    
    def fill_lat_lon(self):
        for item in enumerate(self.new_file):

            ds = xr.open_dataset(item)
            lat_var = ds.variables['lat']
            lon_var = ds.variables['lon']

            lat_values = lat_var.values
            lon_values = lon_var.values
            lat_filled = self.interpolate_lat_lon(lat_values)
            lon_filled = self.interpolate_lat_lon(lon_values)

            ds.variables['lat'] = lat_filled
            ds.variables['lon'] = lon_filled
            ds.close()


def main():
    working_file_list = [

"ncei19_n40x25_w074x00_2018v2_WGS84.nc",
"ncei19_n41x00_w074x00_2015v1_WGS84.nc",
"ncei19_n39x00_w075x00_2018v2_WGS84.nc",  
"ncei19_n40x25_w074x25_2018v2_WGS84.nc",  
"ncei19_n41x00_w074x25_2015v1_WGS84.nc",
"ncei19_n39x00_w075x25_2014v1_WGS84.nc",  
"ncei19_n40x25_w074x75_2014v1_WGS84.nc",  
"ncei19_n41x25_w072x00_2015v1_WGS84.nc",
"ncei19_n39x00_w075x50_2014v1_WGS84.nc",  
"ncei19_n40x25_w075x00_2014v1_WGS84.nc", 
"ncei19_n41x25_w072x25_2015v1_WGS84.nc",
"ncei19_n39x25_w074x75_2018v2_WGS84.nc", 
"ncei19_n40x25_w075x25_2014v1_WGS84.nc",  
"ncei19_n41x25_w072x50_2015v1_WGS84.nc",
"ncei19_n39x25_w075x00_2018v2_WGS84.nc",  
"ncei19_n40x50_w074x00_2018v2_WGS84.nc", 
"ncei19_n41x25_w072x75_2015v1_WGS84.nc",
"ncei19_n39x25_w075x25_2018v2_WGS84.nc",  
"ncei19_n40x50_w074x25_2018v2_WGS84.nc",  
"ncei19_n41x25_w073x00_2016v1_WGS84.nc",
"ncei19_n39x25_w075x50_2014v1_WGS84.nc",  
"ncei19_n40x75_w073x00_2015v1_WGS84.nc",  
"ncei19_n41x25_w073x25_2016v1_WGS84.nc",
"ncei19_n39x50_w074x50_2018v2_WGS84.nc",  
"ncei19_n40x75_w073x25_2015v1_WGS84.nc",  
"ncei19_n41x25_w073x50_2015v1_WGS84.nc",
"ncei19_n39x50_w074x75_2018v2_WGS84.nc",  
"ncei19_n40x75_w073x50_2015v1_WGS84.nc",  
"ncei19_n41x25_w073x75_2015v1_WGS84.nc",
"ncei19_n39x50_w075x25_2018v2_WGS84.nc",  
"ncei19_n40x75_w073x75_2015v1_WGS84.nc",  
"ncei19_n41x25_w074x00_2015v1_WGS84.nc",
"ncei19_n39x50_w075x50_2018v2_WGS84.nc", 
"ncei19_n40x75_w074x00_2015v1_WGS84.nc",  
"ncei19_n41x50_w072x00_2016v1_WGS84.nc",
"ncei19_n39x50_w075x75_2014v1_WGS84.nc",  
"ncei19_n40x75_w074x25_2015v1_WGS84.nc", 
"ncei19_n41x50_w072x25_2016v1_WGS84.nc",
"ncei19_n39x75_w074x25_2018v2_WGS84.nc",  
"ncei19_n41x00_w072x25_2015v1_WGS84.nc",  
"ncei19_n41x50_w072x50_2016v1_WGS84.nc",
"ncei19_n39x75_w074x50_2018v2_WGS84.nc",  
"ncei19_n41x00_w072x50_2015v1_WGS84.nc",  
"ncei19_n41x50_w072x75_2016v1_WGS84.nc",
"ncei19_n39x75_w075x50_2014v1_WGS84.nc",  
"ncei19_n41x00_w072x75_2015v1_WGS84.nc",  
"ncei19_n41x50_w073x00_2016v1_WGS84.nc",
"ncei19_n39x75_w075x75_2014v1_WGS84.nc",  
"ncei19_n41x00_w073x00_2015v1_WGS84.nc", 
"ncei19_n41x50_w074x00_2015v1_WGS84.nc",
"ncei19_n40x00_w074x25_2018v2_WGS84.nc",  
"ncei19_n41x00_w073x25_2015v1_WGS84.nc",  
"ncei19_n41x50_w074x25_2015v1_WGS84.nc",
"ncei19_n40x00_w075x25_2014v1_WGS84.nc",  
"ncei19_n41x00_w073x50_2015v1_WGS84.nc", 
"ncei19_n40x00_w075x50_2014v1_WGS84.nc",  
"ncei19_n41x00_w073x75_2015v1_WGS84.nc"
    ]

    working_new_file_names = [
"ncei19_n40x25_w074x00_2018v2_WGS84_new3.nc",
"ncei19_n41x00_w074x00_2015v1_WGS84_new3.nc",
"ncei19_n39x00_w075x00_2018v2_WGS84_new3.nc",  
"ncei19_n40x25_w074x25_2018v2_WGS84_new3.nc",  
"ncei19_n41x00_w074x25_2015v1_WGS84_new3.nc",
"ncei19_n39x00_w075x25_2014v1_WGS84_new3.nc",  
"ncei19_n40x25_w074x75_2014v1_WGS84_new3.nc",  
"ncei19_n41x25_w072x00_2015v1_WGS84_new3.nc",
"ncei19_n39x00_w075x50_2014v1_WGS84_new3.nc",  
"ncei19_n40x25_w075x00_2014v1_WGS84_new3.nc", 
"ncei19_n41x25_w072x25_2015v1_WGS84_new3.nc",
"ncei19_n39x25_w074x75_2018v2_WGS84_new3.nc", 
"ncei19_n40x25_w075x25_2014v1_WGS84_new3.nc",  
"ncei19_n41x25_w072x50_2015v1_WGS84_new3.nc",
"ncei19_n39x25_w075x00_2018v2_WGS84_new3.nc",  
"ncei19_n40x50_w074x00_2018v2_WGS84_new3.nc", 
"ncei19_n41x25_w072x75_2015v1_WGS84_new3.nc",
"ncei19_n39x25_w075x25_2018v2_WGS84_new3.nc",  
"ncei19_n40x50_w074x25_2018v2_WGS84_new3.nc",  
"ncei19_n41x25_w073x00_2016v1_WGS84_new3.nc",
"ncei19_n39x25_w075x50_2014v1_WGS84_new3.nc",  
"ncei19_n40x75_w073x00_2015v1_WGS84_new3.nc",  
"ncei19_n41x25_w073x25_2016v1_WGS84_new3.nc",
"ncei19_n39x50_w074x50_2018v2_WGS84_new3.nc",  
"ncei19_n40x75_w073x25_2015v1_WGS84_new3.nc",  
"ncei19_n41x25_w073x50_2015v1_WGS84_new3.nc",
"ncei19_n39x50_w074x75_2018v2_WGS84_new3.nc",  
"ncei19_n40x75_w073x50_2015v1_WGS84_new3.nc",  
"ncei19_n41x25_w073x75_2015v1_WGS84_new3.nc",
"ncei19_n39x50_w075x25_2018v2_WGS84_new3.nc",  
"ncei19_n40x75_w073x75_2015v1_WGS84_new3.nc",  
"ncei19_n41x25_w074x00_2015v1_WGS84_new3.nc",
"ncei19_n39x50_w075x50_2018v2_WGS84_new3.nc", 
"ncei19_n40x75_w074x00_2015v1_WGS84_new3.nc",  
"ncei19_n41x50_w072x00_2016v1_WGS84_new3.nc",
"ncei19_n39x50_w075x75_2014v1_WGS84_new3.nc",  
"ncei19_n40x75_w074x25_2015v1_WGS84_new3.nc", 
"ncei19_n41x50_w072x25_2016v1_WGS84_new3.nc",
"ncei19_n39x75_w074x25_2018v2_WGS84_new3.nc",  
"ncei19_n41x00_w072x25_2015v1_WGS84_new3.nc",  
"ncei19_n41x50_w072x50_2016v1_WGS84_new3.nc",
"ncei19_n39x75_w074x50_2018v2_WGS84_new3.nc",  
"ncei19_n41x00_w072x50_2015v1_WGS84_new3.nc",  
"ncei19_n41x50_w072x75_2016v1_WGS84_new3.nc",
"ncei19_n39x75_w075x50_2014v1_WGS84_new3.nc",  
"ncei19_n41x00_w072x75_2015v1_WGS84_new3.nc",  
"ncei19_n41x50_w073x00_2016v1_WGS84_new3.nc",
"ncei19_n39x75_w075x75_2014v1_WGS84_new3.nc",  
"ncei19_n41x00_w073x00_2015v1_WGS84_new3.nc", 
"ncei19_n41x50_w074x00_2015v1_WGS84_new3.nc",
"ncei19_n40x00_w074x25_2018v2_WGS84_new3.nc",  
"ncei19_n41x00_w073x25_2015v1_WGS84_new3.nc",  
"ncei19_n41x50_w074x25_2015v1_WGS84_new3.nc",
"ncei19_n40x00_w075x25_2014v1_WGS84_new3.nc",  
"ncei19_n41x00_w073x50_2015v1_WGS84_new3.nc", 
"ncei19_n40x00_w075x50_2014v1_WGS84_new3.nc",  
"ncei19_n41x00_w073x75_2015v1_WGS84_new3.nc"
    ]
   
    working_surrounding_files = [['ncei19_n40x50_w074x00_2018v2_WGS84.nc', 0, 0, 'ncei19_n40x25_w074x25_2018v2_WGS84.nc', 'ncei19_n40x50_w074x25_2018v2_WGS84.nc', 0, 'ncei19_n40x00_w074x25_2018v2_WGS84.nc', 0], ['ncei19_n41x25_w074x00_2015v1_WGS84.nc', 'ncei19_n41x00_w073x75_2015v1_WGS84.nc', 'ncei19_n40x75_w074x00_2015v1_WGS84.nc', 'ncei19_n41x00_w074x25_2015v1_WGS84.nc', 0, 'ncei19_n41x25_w073x75_2015v1_WGS84.nc', 'ncei19_n40x75_w074x25_2015v1_WGS84.nc', 'ncei19_n40x75_w073x75_2015v1_WGS84.nc'], ['ncei19_n39x25_w075x00_2018v2_WGS84.nc', 0, 0, 'ncei19_n39x00_w075x25_2014v1_WGS84.nc', 'ncei19_n39x25_w075x25_2018v2_WGS84.nc', 'ncei19_n39x25_w074x75_2018v2_WGS84.nc', 0, 0], ['ncei19_n40x50_w074x25_2018v2_WGS84.nc', 'ncei19_n40x25_w074x00_2018v2_WGS84.nc', 'ncei19_n40x00_w074x25_2018v2_WGS84.nc', 0, 0, 'ncei19_n40x50_w074x00_2018v2_WGS84.nc', 0, 0], [0, 'ncei19_n41x00_w074x00_2015v1_WGS84.nc', 'ncei19_n40x75_w074x25_2015v1_WGS84.nc', 0, 0, 'ncei19_n41x25_w074x00_2015v1_WGS84.nc', 0, 'ncei19_n40x75_w074x00_2015v1_WGS84.nc'], ['ncei19_n39x25_w075x25_2018v2_WGS84.nc', 'ncei19_n39x00_w075x00_2018v2_WGS84.nc', 0, 'ncei19_n39x00_w075x50_2014v1_WGS84.nc', 'ncei19_n39x25_w075x50_2014v1_WGS84.nc', 'ncei19_n39x25_w075x00_2018v2_WGS84.nc', 0, 0], [0, 0, 0, 'ncei19_n40x25_w075x00_2014v1_WGS84.nc', 0, 0, 0, 0], ['ncei19_n41x50_w072x00_2016v1_WGS84.nc', 0, 0, 'ncei19_n41x25_w072x25_2015v1_WGS84.nc', 'ncei19_n41x50_w072x25_2016v1_WGS84.nc', 0, 'ncei19_n41x00_w072x25_2015v1_WGS84.nc', 0], ['ncei19_n39x25_w075x50_2014v1_WGS84.nc', 'ncei19_n39x00_w075x25_2014v1_WGS84.nc', 0, 0, 0, 'ncei19_n39x25_w075x25_2018v2_WGS84.nc', 0, 0], [0, 'ncei19_n40x25_w074x75_2014v1_WGS84.nc', 0, 'ncei19_n40x25_w075x25_2014v1_WGS84.nc', 0, 0, 'ncei19_n40x00_w075x25_2014v1_WGS84.nc', 0], ['ncei19_n41x50_w072x25_2016v1_WGS84.nc', 'ncei19_n41x25_w072x00_2015v1_WGS84.nc', 'ncei19_n41x00_w072x25_2015v1_WGS84.nc', 'ncei19_n41x25_w072x50_2015v1_WGS84.nc', 'ncei19_n41x50_w072x50_2016v1_WGS84.nc', 'ncei19_n41x50_w072x00_2016v1_WGS84.nc', 'ncei19_n41x00_w072x50_2015v1_WGS84.nc', 0], ['ncei19_n39x50_w074x75_2018v2_WGS84.nc', 0, 0, 'ncei19_n39x25_w075x00_2018v2_WGS84.nc', 0, 'ncei19_n39x50_w074x50_2018v2_WGS84.nc', 'ncei19_n39x00_w075x00_2018v2_WGS84.nc', 0], [0, 'ncei19_n40x25_w075x00_2014v1_WGS84.nc', 'ncei19_n40x00_w075x25_2014v1_WGS84.nc', 0, 0, 0, 'ncei19_n40x00_w075x50_2014v1_WGS84.nc', 0], ['ncei19_n41x50_w072x50_2016v1_WGS84.nc', 'ncei19_n41x25_w072x25_2015v1_WGS84.nc', 'ncei19_n41x00_w072x50_2015v1_WGS84.nc', 'ncei19_n41x25_w072x75_2015v1_WGS84.nc', 'ncei19_n41x50_w072x75_2016v1_WGS84.nc', 'ncei19_n41x50_w072x25_2016v1_WGS84.nc', 'ncei19_n41x00_w072x75_2015v1_WGS84.nc', 'ncei19_n41x00_w072x25_2015v1_WGS84.nc'], [0, 'ncei19_n39x25_w074x75_2018v2_WGS84.nc', 'ncei19_n39x00_w075x00_2018v2_WGS84.nc', 'ncei19_n39x25_w075x25_2018v2_WGS84.nc', 'ncei19_n39x50_w075x25_2018v2_WGS84.nc', 'ncei19_n39x50_w074x75_2018v2_WGS84.nc', 'ncei19_n39x00_w075x25_2014v1_WGS84.nc', 0], ['ncei19_n40x75_w074x00_2015v1_WGS84.nc', 0, 'ncei19_n40x25_w074x00_2018v2_WGS84.nc', 'ncei19_n40x50_w074x25_2018v2_WGS84.nc', 'ncei19_n40x75_w074x25_2015v1_WGS84.nc', 'ncei19_n40x75_w073x75_2015v1_WGS84.nc', 'ncei19_n40x25_w074x25_2018v2_WGS84.nc', 0], ['ncei19_n41x50_w072x75_2016v1_WGS84.nc', 'ncei19_n41x25_w072x50_2015v1_WGS84.nc', 'ncei19_n41x00_w072x75_2015v1_WGS84.nc', 'ncei19_n41x25_w073x00_2016v1_WGS84.nc', 'ncei19_n41x50_w073x00_2016v1_WGS84.nc', 'ncei19_n41x50_w072x50_2016v1_WGS84.nc', 'ncei19_n41x00_w073x00_2015v1_WGS84.nc', 'ncei19_n41x00_w072x50_2015v1_WGS84.nc'], ['ncei19_n39x50_w075x25_2018v2_WGS84.nc', 'ncei19_n39x25_w075x00_2018v2_WGS84.nc', 'ncei19_n39x00_w075x25_2014v1_WGS84.nc', 'ncei19_n39x25_w075x50_2014v1_WGS84.nc', 'ncei19_n39x50_w075x50_2018v2_WGS84.nc', 0, 'ncei19_n39x00_w075x50_2014v1_WGS84.nc', 'ncei19_n39x00_w075x00_2018v2_WGS84.nc'], ['ncei19_n40x75_w074x25_2015v1_WGS84.nc', 'ncei19_n40x50_w074x00_2018v2_WGS84.nc', 'ncei19_n40x25_w074x25_2018v2_WGS84.nc', 0, 0, 'ncei19_n40x75_w074x00_2015v1_WGS84.nc', 0, 'ncei19_n40x25_w074x00_2018v2_WGS84.nc'], ['ncei19_n41x50_w073x00_2016v1_WGS84.nc', 'ncei19_n41x25_w072x75_2015v1_WGS84.nc', 'ncei19_n41x00_w073x00_2015v1_WGS84.nc', 'ncei19_n41x25_w073x25_2016v1_WGS84.nc', 0, 'ncei19_n41x50_w072x75_2016v1_WGS84.nc', 'ncei19_n41x00_w073x25_2015v1_WGS84.nc', 'ncei19_n41x00_w072x75_2015v1_WGS84.nc'], ['ncei19_n39x50_w075x50_2018v2_WGS84.nc', 'ncei19_n39x25_w075x25_2018v2_WGS84.nc', 'ncei19_n39x00_w075x50_2014v1_WGS84.nc', 0, 'ncei19_n39x50_w075x75_2014v1_WGS84.nc', 'ncei19_n39x50_w075x25_2018v2_WGS84.nc', 0, 'ncei19_n39x00_w075x25_2014v1_WGS84.nc'], ['ncei19_n41x00_w073x00_2015v1_WGS84.nc', 0, 0, 'ncei19_n40x75_w073x25_2015v1_WGS84.nc', 'ncei19_n41x00_w073x25_2015v1_WGS84.nc', 'ncei19_n41x00_w072x75_2015v1_WGS84.nc', 0, 0], [0, 'ncei19_n41x25_w073x00_2016v1_WGS84.nc', 'ncei19_n41x00_w073x25_2015v1_WGS84.nc', 'ncei19_n41x25_w073x50_2015v1_WGS84.nc', 0, 'ncei19_n41x50_w073x00_2016v1_WGS84.nc', 'ncei19_n41x00_w073x50_2015v1_WGS84.nc', 'ncei19_n41x00_w073x00_2015v1_WGS84.nc'], ['ncei19_n39x75_w074x50_2018v2_WGS84.nc', 0, 0, 'ncei19_n39x50_w074x75_2018v2_WGS84.nc', 0, 'ncei19_n39x75_w074x25_2018v2_WGS84.nc', 'ncei19_n39x25_w074x75_2018v2_WGS84.nc', 0], ['ncei19_n41x00_w073x25_2015v1_WGS84.nc', 'ncei19_n40x75_w073x00_2015v1_WGS84.nc', 0, 'ncei19_n40x75_w073x50_2015v1_WGS84.nc', 'ncei19_n41x00_w073x50_2015v1_WGS84.nc', 'ncei19_n41x00_w073x00_2015v1_WGS84.nc', 0, 0], [0, 'ncei19_n41x25_w073x25_2016v1_WGS84.nc', 'ncei19_n41x00_w073x50_2015v1_WGS84.nc', 'ncei19_n41x25_w073x75_2015v1_WGS84.nc', 0, 0, 'ncei19_n41x00_w073x75_2015v1_WGS84.nc', 'ncei19_n41x00_w073x25_2015v1_WGS84.nc'], [0, 'ncei19_n39x50_w074x50_2018v2_WGS84.nc', 'ncei19_n39x25_w074x75_2018v2_WGS84.nc', 0, 0, 'ncei19_n39x75_w074x50_2018v2_WGS84.nc', 'ncei19_n39x25_w075x00_2018v2_WGS84.nc', 0], ['ncei19_n41x00_w073x50_2015v1_WGS84.nc', 'ncei19_n40x75_w073x25_2015v1_WGS84.nc', 0, 'ncei19_n40x75_w073x75_2015v1_WGS84.nc', 'ncei19_n41x00_w073x75_2015v1_WGS84.nc', 'ncei19_n41x00_w073x25_2015v1_WGS84.nc', 0, 0], [0, 'ncei19_n41x25_w073x50_2015v1_WGS84.nc', 'ncei19_n41x00_w073x75_2015v1_WGS84.nc', 'ncei19_n41x25_w074x00_2015v1_WGS84.nc', 'ncei19_n41x50_w074x00_2015v1_WGS84.nc', 0, 'ncei19_n41x00_w074x00_2015v1_WGS84.nc', 'ncei19_n41x00_w073x50_2015v1_WGS84.nc'], [0, 0, 'ncei19_n39x25_w075x25_2018v2_WGS84.nc', 'ncei19_n39x50_w075x50_2018v2_WGS84.nc', 'ncei19_n39x75_w075x50_2014v1_WGS84.nc', 0, 'ncei19_n39x25_w075x50_2014v1_WGS84.nc', 'ncei19_n39x25_w075x00_2018v2_WGS84.nc'], ['ncei19_n41x00_w073x75_2015v1_WGS84.nc', 'ncei19_n40x75_w073x50_2015v1_WGS84.nc', 0, 'ncei19_n40x75_w074x00_2015v1_WGS84.nc', 'ncei19_n41x00_w074x00_2015v1_WGS84.nc', 'ncei19_n41x00_w073x50_2015v1_WGS84.nc', 'ncei19_n40x50_w074x00_2018v2_WGS84.nc', 0], ['ncei19_n41x50_w074x00_2015v1_WGS84.nc', 'ncei19_n41x25_w073x75_2015v1_WGS84.nc', 'ncei19_n41x00_w074x00_2015v1_WGS84.nc', 0, 'ncei19_n41x50_w074x25_2015v1_WGS84.nc', 0, 'ncei19_n41x00_w074x25_2015v1_WGS84.nc', 'ncei19_n41x00_w073x75_2015v1_WGS84.nc'], ['ncei19_n39x75_w075x50_2014v1_WGS84.nc', 'ncei19_n39x50_w075x25_2018v2_WGS84.nc', 'ncei19_n39x25_w075x50_2014v1_WGS84.nc', 'ncei19_n39x50_w075x75_2014v1_WGS84.nc', 'ncei19_n39x75_w075x75_2014v1_WGS84.nc', 0, 0, 'ncei19_n39x25_w075x25_2018v2_WGS84.nc'], ['ncei19_n41x00_w074x00_2015v1_WGS84.nc', 'ncei19_n40x75_w073x75_2015v1_WGS84.nc', 'ncei19_n40x50_w074x00_2018v2_WGS84.nc', 'ncei19_n40x75_w074x25_2015v1_WGS84.nc', 'ncei19_n41x00_w074x25_2015v1_WGS84.nc', 'ncei19_n41x00_w073x75_2015v1_WGS84.nc', 'ncei19_n40x50_w074x25_2018v2_WGS84.nc', 0], [0, 0, 'ncei19_n41x25_w072x00_2015v1_WGS84.nc', 'ncei19_n41x50_w072x25_2016v1_WGS84.nc', 0, 0, 'ncei19_n41x25_w072x25_2015v1_WGS84.nc', 0], ['ncei19_n39x75_w075x75_2014v1_WGS84.nc', 'ncei19_n39x50_w075x50_2018v2_WGS84.nc', 0, 0, 0, 'ncei19_n39x75_w075x50_2014v1_WGS84.nc', 0, 'ncei19_n39x25_w075x50_2014v1_WGS84.nc'], ['ncei19_n41x00_w074x25_2015v1_WGS84.nc', 'ncei19_n40x75_w074x00_2015v1_WGS84.nc', 'ncei19_n40x50_w074x25_2018v2_WGS84.nc', 0, 0, 'ncei19_n41x00_w074x00_2015v1_WGS84.nc', 0, 'ncei19_n40x50_w074x00_2018v2_WGS84.nc'], [0, 'ncei19_n41x50_w072x00_2016v1_WGS84.nc', 'ncei19_n41x25_w072x25_2015v1_WGS84.nc', 'ncei19_n41x50_w072x50_2016v1_WGS84.nc', 0, 0, 'ncei19_n41x25_w072x50_2015v1_WGS84.nc', 'ncei19_n41x25_w072x00_2015v1_WGS84.nc'], ['ncei19_n40x00_w074x25_2018v2_WGS84.nc', 0, 0, 'ncei19_n39x75_w074x50_2018v2_WGS84.nc', 0, 0, 'ncei19_n39x50_w074x50_2018v2_WGS84.nc', 0], ['ncei19_n41x25_w072x25_2015v1_WGS84.nc', 0, 0, 'ncei19_n41x00_w072x50_2015v1_WGS84.nc', 'ncei19_n41x25_w072x50_2015v1_WGS84.nc', 'ncei19_n41x25_w072x00_2015v1_WGS84.nc', 0, 0], [0, 'ncei19_n41x50_w072x25_2016v1_WGS84.nc', 'ncei19_n41x25_w072x50_2015v1_WGS84.nc', 'ncei19_n41x50_w072x75_2016v1_WGS84.nc', 0, 0, 'ncei19_n41x25_w072x75_2015v1_WGS84.nc', 'ncei19_n41x25_w072x25_2015v1_WGS84.nc'], [0, 'ncei19_n39x75_w074x25_2018v2_WGS84.nc', 'ncei19_n39x50_w074x50_2018v2_WGS84.nc', 0, 0, 'ncei19_n40x00_w074x25_2018v2_WGS84.nc', 'ncei19_n39x50_w074x75_2018v2_WGS84.nc', 0], ['ncei19_n41x25_w072x50_2015v1_WGS84.nc', 'ncei19_n41x00_w072x25_2015v1_WGS84.nc', 0, 'ncei19_n41x00_w072x75_2015v1_WGS84.nc', 'ncei19_n41x25_w072x75_2015v1_WGS84.nc', 'ncei19_n41x25_w072x25_2015v1_WGS84.nc', 0, 0], [0, 'ncei19_n41x50_w072x50_2016v1_WGS84.nc', 'ncei19_n41x25_w072x75_2015v1_WGS84.nc', 'ncei19_n41x50_w073x00_2016v1_WGS84.nc', 0, 0, 'ncei19_n41x25_w073x00_2016v1_WGS84.nc', 'ncei19_n41x25_w072x50_2015v1_WGS84.nc'], ['ncei19_n40x00_w075x50_2014v1_WGS84.nc', 0, 'ncei19_n39x50_w075x50_2018v2_WGS84.nc', 'ncei19_n39x75_w075x75_2014v1_WGS84.nc', 0, 'ncei19_n40x00_w075x25_2014v1_WGS84.nc', 'ncei19_n39x50_w075x75_2014v1_WGS84.nc', 'ncei19_n39x50_w075x25_2018v2_WGS84.nc'], ['ncei19_n41x25_w072x75_2015v1_WGS84.nc', 'ncei19_n41x00_w072x50_2015v1_WGS84.nc', 0, 'ncei19_n41x00_w073x00_2015v1_WGS84.nc', 'ncei19_n41x25_w073x00_2016v1_WGS84.nc', 'ncei19_n41x25_w072x50_2015v1_WGS84.nc', 'ncei19_n40x75_w073x00_2015v1_WGS84.nc', 0], [0, 'ncei19_n41x50_w072x75_2016v1_WGS84.nc', 'ncei19_n41x25_w073x00_2016v1_WGS84.nc', 0, 0, 0, 'ncei19_n41x25_w073x25_2016v1_WGS84.nc', 'ncei19_n41x25_w072x75_2015v1_WGS84.nc'], [0, 'ncei19_n39x75_w075x50_2014v1_WGS84.nc', 'ncei19_n39x50_w075x75_2014v1_WGS84.nc', 0, 0, 'ncei19_n40x00_w075x50_2014v1_WGS84.nc', 0, 'ncei19_n39x50_w075x50_2018v2_WGS84.nc'], ['ncei19_n41x25_w073x00_2016v1_WGS84.nc', 'ncei19_n41x00_w072x75_2015v1_WGS84.nc', 'ncei19_n40x75_w073x00_2015v1_WGS84.nc', 'ncei19_n41x00_w073x25_2015v1_WGS84.nc', 'ncei19_n41x25_w073x25_2016v1_WGS84.nc', 'ncei19_n41x25_w072x75_2015v1_WGS84.nc', 'ncei19_n40x75_w073x25_2015v1_WGS84.nc', 0], [0, 0, 'ncei19_n41x25_w074x00_2015v1_WGS84.nc', 'ncei19_n41x50_w074x25_2015v1_WGS84.nc', 0, 0, 0, 'ncei19_n41x25_w073x75_2015v1_WGS84.nc'], ['ncei19_n40x25_w074x25_2018v2_WGS84.nc', 0, 'ncei19_n39x75_w074x25_2018v2_WGS84.nc', 0, 0, 'ncei19_n40x25_w074x00_2018v2_WGS84.nc', 'ncei19_n39x75_w074x50_2018v2_WGS84.nc', 0], ['ncei19_n41x25_w073x25_2016v1_WGS84.nc', 'ncei19_n41x00_w073x00_2015v1_WGS84.nc', 'ncei19_n40x75_w073x25_2015v1_WGS84.nc', 'ncei19_n41x00_w073x50_2015v1_WGS84.nc', 'ncei19_n41x25_w073x50_2015v1_WGS84.nc', 'ncei19_n41x25_w073x00_2016v1_WGS84.nc', 'ncei19_n40x75_w073x50_2015v1_WGS84.nc', 'ncei19_n40x75_w073x00_2015v1_WGS84.nc'], [0, 'ncei19_n41x50_w074x00_2015v1_WGS84.nc', 0, 0, 0, 0, 0, 'ncei19_n41x25_w074x00_2015v1_WGS84.nc'], ['ncei19_n40x25_w075x25_2014v1_WGS84.nc', 0, 0, 'ncei19_n40x00_w075x50_2014v1_WGS84.nc', 0, 'ncei19_n40x25_w075x00_2014v1_WGS84.nc', 'ncei19_n39x75_w075x50_2014v1_WGS84.nc', 0], ['ncei19_n41x25_w073x50_2015v1_WGS84.nc', 'ncei19_n41x00_w073x25_2015v1_WGS84.nc', 'ncei19_n40x75_w073x50_2015v1_WGS84.nc', 'ncei19_n41x00_w073x75_2015v1_WGS84.nc', 'ncei19_n41x25_w073x75_2015v1_WGS84.nc', 'ncei19_n41x25_w073x25_2016v1_WGS84.nc', 'ncei19_n40x75_w073x75_2015v1_WGS84.nc', 'ncei19_n40x75_w073x25_2015v1_WGS84.nc'], [0, 'ncei19_n40x00_w075x25_2014v1_WGS84.nc', 'ncei19_n39x75_w075x50_2014v1_WGS84.nc', 0, 0, 'ncei19_n40x25_w075x25_2014v1_WGS84.nc', 'ncei19_n39x75_w075x75_2014v1_WGS84.nc', 0], ['ncei19_n41x25_w073x75_2015v1_WGS84.nc', 'ncei19_n41x00_w073x50_2015v1_WGS84.nc', 'ncei19_n40x75_w073x75_2015v1_WGS84.nc', 'ncei19_n41x00_w074x00_2015v1_WGS84.nc', 'ncei19_n41x25_w074x00_2015v1_WGS84.nc', 'ncei19_n41x25_w073x50_2015v1_WGS84.nc', 'ncei19_n40x75_w074x00_2015v1_WGS84.nc', 'ncei19_n40x75_w073x50_2015v1_WGS84.nc']]
    working_class = string_parse_function(working_file_list)
    working_class.create_netcdf()
    working_class.cut_string()
    working_class.check_boxes()
    working_class.add_data()
    # new_file_names = working_class.new_file_names
    # surrounding_file_names = working_class.surrounding_file_names
    working_class2 = overlapping_data(working_surrounding_files, working_new_file_names)
    working_class2.eliminate_repeats()
    working_class2.fill_lat_lon()

if __name__ == main():
    main()