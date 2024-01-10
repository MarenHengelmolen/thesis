
This prototype is developed to help users prepare their CFD simulations for urban areas in the open source CFD software OpenFOAM. 

First, it identifies geometric errors that might affect the meshing process essential for accurate results. Subsequently, it proposes mesh parameters based on recent CFD guidelines for rectangular domains. 

The following geometric validations are performed: 

* Separate building and terrain validation by the [val3dity](https://github.com/tudelft3d/val3dity) tool,
* Topological relationships between buildings and terrain validation,
* Sharp angles identification,
* Short edges identification,
* Sliver triangles identification,
* Overlapping buildings identification.

The results are provided in the form of reports, OBJ files and/or TXT files.

To help users with defining adequate mesh parameters, the prototype returns the following output: 

* The 3D model aligned with the incoming flow simulated in OpenFOAM (OBJ),
* Configuration files containing mesh parameters for CFD simulation in urban areas (blockMeshDict and snappyHexMeshDict),
* Advice and information on the input model aiming to help users find a balance between simulation accuracy and performance (e.g. number of buildings satisfying one
of the guidelines, buildings needed to run a realistic simulation). This is provided in the form of OBJ and TXT files.

## Data formats 
The prototype supports OBJ and STL files. 

## Getting started 

* Open Ubuntu
* Go to the website directory
* Run ```flask --app website --debug run```
* Open the generated link
  
First, you need to insert your input model. Then, the prototype leads you to another page in which you need to insert some parameters. Explanation is given for each parameter by clicking on them. After inserting your parameters, the prototype performs the geometric validations and mesh parameters definition previously addressed and returns the results. 

