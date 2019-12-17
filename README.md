# LW_template_correlation
A Matlab implementation of the moving-window template correlation method for detection of hippocampal replay during REM sleep, as used in Louie &amp; Wilson 2001.

=========== INSTRUCTIONS FOR USE =============================
1. copy files into a directory, and set the directory as your PWD. 
2. Open the example script m-file, "LW_example_script_1.m" and change the load path so that it points to your data
3. Run the script. As the script was designed to run on an HPC using parallel workers, this will be slow on a native machine, so for a "quick" demo, set mode to 1 (debug).
