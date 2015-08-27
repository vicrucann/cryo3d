### Short description  
*cryo3d* is a matlab-based software that allows for fast 3D protein reconstruction based on cryogenic particle images. The software was developped at [Image Processing and Analysis Group](http://medicine.yale.edu/bioimaging/ipa/), Yale School of Medicine, Dept. of Diagnostic Radiology under supervision and with participation of **Prof. Hemant Tagare**. 

### The project structure  
`/doc` : contains project description-related documents    
`/script` : main `*.m` scripts of the pipeline and main workspace (the Matlab should be run from that folder)  
`/src` : `*.m` functions, arranged by the different pipeline step, e.g., preprocessing, best_match, etc; `mrc` folder contains function to read-write `*.mrc` files  
`/test` : temporal functions that are under development, subject to delete for the final release version   

### Source code  

To download the source code of the pipeline, run the following commands:  
```
git clone https://github.com/vicrucann/cryo3d
``` 
To get the submodules (cacharr and rshell-mat):  
```
cd cryo3d
git submodule update --init --recursive
```
When planning to edit (and then push back to the github) any of the submodules, checkout to the master branch for each of the submodules:  
```
cd src/rshell-mat
git checkout master
cd ..
cd cacharr
git checkout master
```

### Running the whole pipeline

**User must provide configuration file as an input for the `\script\cryo3d.m` function (you can find an example of such file for both Windows and Linux in `script` folder.**    

### Main preprocessing components and their outputs
```
	1. (*.star, *.mrcs)     	->	fit_ctfs.m 	            -> fit-ctfs.mat
	2. (fit-ctfs.mat, *.mrcs) 	-> 	proprocess_images.m 	-> pre-imgs.mrcs
	3. (*.mrc)	            	->	init_volume.m   		-> ini-vol.mrc
	4. (theta)          		->	coordinate_axes.m       -> coord-ax.mat
```

### Author information
All of function for fast best match method were written by Nicha Dvornek. The improvements like main framework organization, distributor package [rshell-mat](https://github.com/vicrucann/rshell-mat) and caching class [cacharr](https://github.com/vicrucann/cacharr) were written by Victoria Rudakova. For inquires on additional packages: vicrucann(at)gmail(dot)com. 
