# NetCDF-Indexing-Template
A python program for adding bathymetry data from surrounding NetCDF files without having repeated data

This program:
1. Converts a string name into geographic coordinates
2. Indexes data from surrounding coordinates within the NetCDF file.
3. Removes duplicate entries from the file to ensure clean, consistent datasets.

---


## Features
- **Name-to-coordinate conversion** — Input a list of file names that are converted to coordinates rather than manually inputting latitude and longitude
- **Spatial indexing** — Sorts through list and finds surrounding coordinates for indexing into newly sized NetCDF
- **Duplicate removal** — Removes overlaps after indexing surrounding data
- **Batch support** — Works on list of NetCDF files in a directory

---

## Installation
Clone the repository and install dependencies:

```bash
git clone https://github.com/mtitterton/NetCDF-Indexing-Template.git
cd NetCDF-Indexing-Template
pip install -r requirements.txt
```

## Usage
You can find detailed usage instructions in the docs/usage.md file.

## Requirements
- Python 3.8+
- netCDF4
- numpy
- pandas
- xarray
- os
