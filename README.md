This repository contains scripts that automate the setup of finite element simulations of stent expansion in patient-specific geometries of total cavopulmonary connection in Abaqus CAE v 2020. Three pairs of files can be used to create a simulation in which the stent is expanded to either straighten the patient geometry ('striaghten' files) or conform to the patient geometry without modelling implantation ('virtual catheter' files), or to confrom to the patient geometry with modelling implantation ('tracking' files). For each pair of files, one facilitates preprocessing of the geometry data to gether information that allows the model to be built ('preprocessing' file), and the second can be launched through the Abaqus scripting interface to generate the finite element model ('ASI Build' file). 

An overview of code useage can be found here: 
[Code Overview.pdf](https://github.com/user-attachments/files/16910242/Code.Overview.pdf)

More information on the background of the project can be found here: 
[Project_Report.pdf](https://github.com/user-attachments/files/16910252/Project_Report.pdf)
