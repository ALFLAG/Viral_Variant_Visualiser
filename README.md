This is Viral Variant Visualiser - VVV

from single-ends or paired-end NGS data, it allows the recontruction of the most present viral genome and determine the genetic variant composition of the viral population, by outputting a PNG image of the viral population and other complementary results files.


This tool can be launch in two ways:     
    . command line: ```python3 vvv_directory/SCRIPTS_PYTHON/viral_variant_visualiser.py``` # with all the requiered argument, see the help section. In that case, make sure you comment lines 133, 139, 141 and uncomment lines 134, 140, 142.   
    . via GALAXY: the file ```viral_variant_visualiser.xml``` can be added to a GALAXY instance together with the ```vvv_directory/SCRIPTS_PYTHON/viral_variant_visualiser.py``` script. In this case, make sure you update the command section of the xml file, with the correct path to the script, and the correct path to ```vvv_directory/```


To use this tool, clone this repository and modify the PATH variable in the ```vvv_directory/snakefile_ViralVariantVisualiser``` into the path were vvv_directory was placed.  
This tool requieres snakemake 5.10 installed in a conda environment named ```snakes_update``` that can be found in the ```vvv_directory/ENVS/``` directory.  
This tool requieres the configuration of a megablast database to determine the metagenomic composition of the sample (nt database from NCBI for example, or RefSeq). And then, modify the BLASTDB and ANSES_NT_DIR variable to the database you used.
This tool requieres the configuration of a custom database using the files of the VRL section of GenBank (viral sequences, naturally found in nature, no synthetic sequence). To do that, launch the script ```vvv_directory/SCRIPTS_VDB/initiate_vdb.sh``` in the ```vvv_directory/SCRIPTS_VDB/``` after activating the ```snakes_update``` environment.  

When you use this tool, please cite:  
A.Flageul, P.Lucas, E.Hirchaud, F.Touzain, Y.Blanchard, N.Eterradossi, P.Brown, B.Grasland, 2020, Viral Variant Visualizer - A novel bioinformatic tool to quickly and simply visualize viral genetic diversity from next generation sequencing raw data, Virus Research
