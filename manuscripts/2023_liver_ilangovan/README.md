# Translational Radiation Research (TRRaC) - Spaced Out Data No More: Genomic Harmonization Meets Machine Learning in Murine Livers

The purpose of this repository is to provide the necessary code to perform our analysis and generate figures necessary for our the publication *Spaced Out Data No More: Genomic Harmonization Meets Machine Learning in Murine Livers*. 

# Authorship Details
- Repository Owner: Hari Ilangovan, Senior Data Scientist, SAIC, NASA Information, Data and Analytics Services, NASA Langley Research Center


|Version History | Date | 
|----------| ----- |
|v0| 12/10/2022 | 
|v1 | 11/2/2023 | 


Publication Authorship:
- Hari Ilangovan<sup>1</sup>
- Prachi Kothiyal<sup>2</sup>
- Katherine A. Hoadley<sup>3</sup>
- Robin Elgart<sup>4</sup>
- Greg Eley<sup>2</sup>
- Parastou Eslami<sup>5</sup>

<sup>1</sup> Science Applications International Corporation (SAIC), Reston, VA 20190, USA
<sup>2</sup>Scimentis LLC, Statham, GA 30666, USA
<sup>3</sup>Department of Genetics, Computational Medicine Program, Lineberger Comprehensive Cancer Center, University of North Caroline at Chapel Hill, Chapel Hill, NC 27599, USA
<sup>4</sup>University of Houston, Houston, TX 77204, USA
<sup>5</sup>Universal Artificial Intelligence Inc, Boston, MA 02130, USA

# Directories

The description for directories are as follows:
- `environments` directory contains the library versions for the open source Python and R modules used to execute the Jupyter notebooks and . `conda` can be used to manage separate environments for R and Python. Please refer to the [anaconda](https://docs.anaconda.com/free/anaconda/install/) page for further information on environment setup and configuration.
- `notebooks` directory contains the analysis notebooks used to transform source data and execute analysis described in the publication.
- `scripts` directory contains the helper functions that support the function calls from the R & Python Jupyter notebooks.

# Downloading the Source Data
All data used in this analysis can be obtained directly from Open Science Data Repository (OSDR). The following data set pages can be accessed to download the unnormalized RNA-seq counts data and their associated metadata files:
- [OSD-47](https://osdr.nasa.gov/bio/repo/data/studies/OSD-47)
- [OSD-168](https://osdr.nasa.gov/bio/repo/data/studies/OSD-168)
- [OSD-242](https://osdr.nasa.gov/bio/repo/data/studies/OSD-242)
- [OSD-245](https://osdr.nasa.gov/bio/repo/data/studies/OSD-245)
- [OSD-379](https://osdr.nasa.gov/bio/repo/data/studies/OSD-379)

The source data files can be used to overwrite the placeholder counts matrix and metadata files in the `./data/raw_counts/` and `./data/metadata/` directories, respectfully. Metadata files for each study are already populated for convenience. The counts need to be downloaded with the latest from OSDR.

Note that [OSD-168](https://osdr.nasa.gov/bio/repo/data/studies/OSD-168) consists of two rodent research missions (Rodent Research 1 and Rodent Research 3). These should be segmented by mission with the ID `168_rr1` an `168_rr3`. The mission associated with each sample can be found from the metadata. 

The filenames of placeholders are designed to match the calls in the notebooks and scripts. If file names are changed from defaults, the corresponding sections in the source code will also need to be updated.

# Running the Notebooks

After source data are downloaded, the notebooks should be executed in the following order:
1. Pseudogene Filter Gene Symbol Conversion
    - Generates the `glds_no_pseudogenes_wExternal.csv` using a BioMaRt query that consists of protein coding genes with external symbols that are used to subset the complete list of genes assayed from the OSD data sets. 
2. Pre-filtering
    - Generates the `prefiltered_pseudogenes.csv` file that contains the set of gene features that will be used to subset each counts matrix in the `Python Mouse Liver Manuscript - Results` and `R Mouse Liver Manuscript - Results` notebooks.
3. Python Mouse Liver Manuscript - Results
    - The table of contents and in-line comments describe the various files required to execute the analysis code.
4. R Mouse Liver Manuscripts - Results
    - The table of contents and in-line comments describe the various files required to execute the analysis code and generate GSEA results.

