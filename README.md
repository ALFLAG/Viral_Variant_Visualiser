This is Viral Variant Visualiser - VVV

from single-ends or paired-end NGS data, it allows the recontruction of the most present viral genome and determine the genetic variant composition of the viral population, by outputting a PNG image of the viral population and other complementary results files.

To use this tool, clone this repository.   
This tool can be launch in two ways:     
    . command line: ```python3 vvv_directory/SCRIPTS_PYTHON/launch_VVV_v2.py``` # with all the requiered argument, see the help section  
    . via GALAXY: the file ```viral_variant_visualiser.xml``` can be added to a GALAXY instance together with the ```vvv_directory/SCRIPTS_PYTHON/launch_VVV_v2.py``` script. In this case, make sure you update the command section of the xml file, with the correct path to the script, and the correct path to ```vvv_directory/```

After downloading all this repository, modify the PATH variable in the ```vvv_directory/snakefile_viral_variant_visualiser``` into the path were vvv_directory was placed.  
This tool requieres snakemake 5.10 installed in a conda environment named ```snakes_update``` that can be found in the ```vvv_directory/ENVS/``` directory.  
This tool requieres the configuration of a megablast database to determine the metagenomic composition of the sample (nt database from NCBI for example, or RefSeq). And then, modify the BLASTDB and ANSES_NT_DIR variable to the database you used.
This tool requieres the configuration of a custom database using the files of the VRL section of GenBank (viral sequences, naturally found in nature, no synthetic sequence). To do that, launch the script ```vvv_directory/SCRIPTS_VDB/initiate_vdb.sh``` in the vvv_directory after activating the ```snakes_update``` environment.  
