### Short description  
*cryo3d* is a matlab-based software that allows for fast 3D protein reconstruction based on cryogenic particle images. The software was developped at [Image Processing and Analysis Group](http://medicine.yale.edu/bioimaging/ipa/), Yale School of Medicine, Dept. of Diagnostic Radiology. 

### The project structure  
*/doc* : contains project description and structure  
*/script* : *.m scripts of the pipeline  
*/src* : *.m functions, arranged by the different pipeline step, e.g., preprocessing, best_match, etc; mrc folder contains function to read-write *.mrc files  
*/test* : temporal functions that are under development  

### Running the whole pipeline

**User must provide configuration file as an input for the `\script\cryo3d.m` function (you can find an example of such file for both Windows and Linux in `script` folder.**    

### Main preprocessing components and their outputs
	1. (...)		->	fit_ctfs.m 		-> fit-ctfs.mat
	2. fit-ctfs.mat 	-> 	proprocess_images.m 	-> pre-imgs.mrcs
	3. (...)		->	init_volume.m 		-> ini-vol.mrc
	4. (...)		->	coordinate_axes.m 	-> coord-ax.mat

### Author information
The main skeleton of the software was implemented by Nicha Dvornek, some functions were done by Hemant Tagare, Victoria Rudakova. 

