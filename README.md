# HydroQC Mini (EMODnet Bathymetry QC)

A small Python QC tool that generates quick quality checks and diagnostic maps for an EMODnet bathymetry/DTM GeoTIFF.

## What this project does
Given an input GeoTIFF, it produces:
- `outputs/report.md` (QC summary report)
- `outputs/summary.json` (machine-readable stats)
- `outputs/figures/` (PNG maps: depth, missing-data mask, outliers, slope, roughness)

## Data source (not included in this repo)
The raw GeoTIFF is not uploaded to GitHub to keep the repository lightweight.
Download the tile(s) from EMODnet GeoViewer:
https://emodnet.ec.europa.eu/geoviewer/#!/

After download, place the GeoTIFF in:
`data_raw/`

Example:
`data_raw/D5_2024.tif`

## Setup
Create/activate your environment, then install dependencies:

```bash
pip install -r requirements.txt

Run
python src/run_qc.py data_raw/D5_2024.tif --out outputs
