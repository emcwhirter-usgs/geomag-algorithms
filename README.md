Geomag Algorithms
=================

Geomag algorithms includes tools to fetch, process, and output geomag data.

## Supported Formats ##

### [Edge](readme_io.md#edge) ###

Use an Edge server for data input.

### [IAGA](readme_io.md#iaga2002)###

Use IAGA2002 formatted files for input and/or output.

### [PCDCP](readme_io.md#) ###

Use PCDCP formatted files for input and/or output.


---
## Supported Algorithms ##

### [DeltaF](./docs/DeltaF_usage.md) ###
Calculate DeltaF from geographic, or observatory coordinates.

### [K-USGS](./docs/KUSGS_usage.md) ###
Calculate approved USGS version of K-Index. Some history behind K-USGS.

### [XYZ](./docs/XYZ_usage.md) ###
Rotate data between coordinate systems. From HEZ or HDZ to XYZ and back.


---
## Getting Started ##

### [Install](readme_dependency_install.md) ###
First time install. Walk through dependencies and other considerations.

### [Use](readme_usage.md) ###
Details and examples for proper usage. Get started quickly.

Basic usage:

  - Use the main script, `geomag.py -h`
  - In python scripts, `import geomagio`

### [Develop](readme_develop_install.md) ###
Development dependencies discussed here. Project is built with Grunt and Node
and is written primarily in Python 2.7.
