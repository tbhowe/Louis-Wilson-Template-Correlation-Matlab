# LW_template_correlation
A Matlab implementation of the moving-window template correlation method for detection of hippocampal replay during REM sleep, as used in Louie &amp; Wilson 2001.

=========== INSTRUCTIONS FOR USE =============================
1. copy files from "code" directory into a directory on your machine, and set the directory as your PWD. 
2. Open the example script m-file, "LW_example_script_1.m" and change the load path so that it points to your data. Example data is included in the folder called "test data"
3. Run the script. As the script was designed to run on an HPC using parallel workers, this will be slow on a native machine, so for a "quick" demo, set mode to 1 (debug).


========== DESCRIPTION OF METHOD =============================

A full description of the method implemented here can be found in the following Journal article, which is open access:

Louie K, Wilson MA.
Temporally structured replay of awake hippocampal ensemble activity during rapid eye movement sleep.
Neuron. 2001 Jan;29(1):145-56.
https://www.ncbi.nlm.nih.gov/pubmed/11182087
