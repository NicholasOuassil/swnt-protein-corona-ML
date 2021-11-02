# swnt-protein-corona-ML
Ensemble classifiers and data scraping tools for the prediction of protein surface adsorption to single-walled carbon nanotubes

This work is published as "Supervised Learning Model Predicts Protein Adsorption to Carbon Nanotubes" by Nicholas Ouassil*, Rebecca L. Pinals*, Jackson Travis Del Bonis-O'Donnell, Jeffery W. Wang, and Markita P. Landry in [*biorXiv*](https://doi.org/10.1101/2021.06.19.449132) 

*Co-Authors 


## Included Files and What They Do

1. **Prep for NetSurfP** (Jupyter Notebook) 

    * This script collects data for processing with **NetSurfP 2.0** it will generate a text file. Copy this text file into the prompt on the NetSurfP 2.0 website. Download data and place into an excel sheet then revisit this script and let it process the excel sheet. 

2. **Power Study Script** (Jupyter Notebook)

   * Reproduces the labeling scheme used in figure 1b of the published work.
   * *Requires* the output of __Prep For NetsurP__ 

3. **Data Workup Script** (Jupyter Notebook)
    
   * This script is used to do the majority of the prepping for running through 
   * *Requires* the output of __Prep For NetsurP__ and Power Value from __Power Study Script__

4. **Classification Script** (Jupyter Notebook)

    * Reproduces the classification experiments from the manuscript
    * *Requires* the output of __Data Workup Script__ 


5. **Hyperparameter Optimization** (Jupyter Notebook)

    * Uses GridSearchCV to find the best hyperparameters. Be Careful here as runtimes can rapidly increase with too many parameters
    * *Requires* the output of __Data Workup Script__ 

6. **Different Classifier Testing** (Jupyter Notebook)

    * Uses Generic Architectures to Identify the Best Architecture for the Data
    * *Requires* the output of __Data Workup Script__ 

6. **Data Prep Functions, Interpro Scraping, UniProt NetSurfP Scraping** (Python Scripts)

    * Used throughout the other notebooks in order to get to results
    * Code can be adjusted in these files in order to enhanced feature mining

