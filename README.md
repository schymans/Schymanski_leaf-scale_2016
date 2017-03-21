# Schymanski_leaf-scale_2016

Code and data to reproduce computational results and plots presented in: 
[**Leaf-scale experiments reveal important omission in the Penman-Monteith equation.**](http://www.hydrol-earth-syst-sci-discuss.net/hess-2016-363)
Hydrol. Earth Syst. Sci. Discuss., doi:10.5194/hess-2016-363, 2017

You can either access the files directly on the cloud.sagemath.com server [here]( https://cloud.sagemath.com/projects/e66470cf-1fa4-48bc-8a49-2513f177cf8a/files/Schymanski_leaf-scale_2016/), or follow the below instructions.

## Instructions on how to use these files
To use these files, you can either:

- Download all files as zip file and extract to a local directory on your computer.
- Install [sage](http://www.sagemath.org) on your local computer and start the ipython notebook to work with the files, 

or:

- Create an account at [sagemath cloud](https://cloud.sagemath.com)
- Create a new project
- Start newly created project
- Click on "Create or upload files..."
- Create folders `data`, `figures` and `temp`
- Either:
    - Drag and drop all the downloaded and extracted files onto the field "Drop files to upload" of the respective folders
    - Click on "Files" at the top left of the browser window to see the files. Click on any of them to execute and edit.
- or:
    - To reproduce all results, you just need to create the above folders, and only copy files in the base folder and in `data`, then execute all cells in Worksheet `Worksheet_update`. 


or:

- Access the published worksheets directly on the cloud.sagemath.com server [here]( https://cloud.sagemath.com/projects/e66470cf-1fa4-48bc-8a49-2513f177cf8a/files/Schymanski_leaf-scale_2016/)

## Description of folders
### Main folder
The main folder contains the various worksheets (.ipynb files), the README.md and licence information. 
### data
The data folder contains data files accessed by the worksheets
### figures
Figures produced by the worksheets are placed in the figures worksheet. They are equivalent to those in the paper.
### temp
The temp folder contains temporary files created by the worksheets for exchange of data and variables between worksheets, such as .sage and .sobj files.



## Description of files
### .sage files
These files consist of input code extracted from the .ipynb worksheets with the same name. They are used internally to re-use variables and equations defined in various worksheets.
### .sobj files
These files contain snapshots of worksheet data for quick import into other worksheets.
### Worksheet_setup.ipynb
Contains code to setup the other worksheets. Only relevant if you are interested in the inner workings of the code.
### Worksheet_update.ipynb
Converts Worksheet_setup.ipynb to a .sage file, so that it can be more easily imported by other worksheets.
### leaf_enbalance_eqs.ipynb
Contains variable definitions and equations related to the numerical leaf energy balance model.
### E_PM_eqs.ipynb
Contains variable definitions and equations related to the analytical solutions of the leaf energy balance. It builds on definitions and equations provided in leaf_enbalance_eqs.ipynb
### leaf_chamber_eqs.ipynb
Equations to compute wind tunnel energy balance etc. It builds on definitions and equations provided in leaf_enbalance_eqs.ipynb
### stomatal_cond_eqs.ipynb
Equations to compute stomatal conductance of perforated foils
### leaf_chamber_data.ipynb
Evaluation and plotting of experimental data from leaf wind tunnel and additional data plots.

