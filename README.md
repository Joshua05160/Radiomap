# Radiomap

> Open-source code for Radio Map-Enabled 3D Trajectory and Communication
Optimization for Low-Altitude Air-Ground Cooperation （GC-Workshop-2026）
> 
## Project Overview

Radiomap is the open-source code for the paper "Radio Map-Enabled 3D Trajectory and Communication Optimization for Low-Altitude Air-Ground Cooperation". The primary functionalities include extracting building information from OpenStreetMap (OSM) data, extracting wireless channel information from Ranplan, and performing problem modeling and algorithm solving.

## Features

- 🗺️ OSM Data Processing: Load and process geographic data from OpenStreetMap files
- 🏢 Building Extraction: Automatically extract and visualize building geometry information
- 📡 Radio Modeling: Perform radio signal propagation analysis based on geographical environment

## Project Structure

```
Radiomap/
├── Algorithm/              # Core algorithm implementation
├── AuxiliaryFunction/      # Auxiliary function library
├── extractAndVisualizeBuildingsFromOSM.m  # Main program for OSM building extraction
└── README.md              # Project documentation
```

## Basic Usage

1. **Prepare OSM data file**: Download the required area `.osm` file from [OpenStreetMap](https://www.openstreetmap.org/)

2. **Run building extraction**:
```matlab
% Example code
osmFile = 'your_area.osm';
buildings = extractAndVisualizeBuildingsFromOSM(osmFile);
```
3. **View visualization results**: The program will automatically generate a map visualization window
4. **Prepare Radiomap Data**: Import the `.osm` file to [Ranplan](https://www.ranplanwireless.com/), perform detailed modeling refinement as needed, then export the radiomap data.
## Algorithm Description

The project employs the following core algorithms:

- **SCA-Based Communication Scheduling Optimization**
- **SCA-Based Power Optimization**
- **SCA-Based UAVs' Trajectory Initialization**
- **PSO-CM-Based UAVs' Trajectory Optimization**

## Acknowledgments

- Thanks to the OpenStreetMap community for providing open geographic data
- Thanks to the Ranplan Wireless for their wireless planning and optimization software

### Latest Version (2025-09-13)
- ✅ Implemented OSM file loading and processing functionality
- ✅ Added building extraction and visualization module
- ✅ Established project infrastructure
