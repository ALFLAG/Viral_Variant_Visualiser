This is Viral Variant Visualiser - VVV

from single-ends or paired-end NGS data, it allows the recontruction of the most present viral genome and determine the genetic variant composition of the viral population, by outputting a PNG image of the viral population and oter complementary results files.

This tool can be lauch in two ways:  
    command line: ```python3 vvv_directory/SCRIPTS_PYTHON/launch_VVV_v2.py``` # with all the requiered argument, see the help section  
    via galaxy: ```viral_variant_visualiser.xml``` file can be added to galaxy together with the ```SCRIPTS_PYTHON/launch_VVV_v2.py``` script  

After cloning this tool, modify the PATH variable in the ```vvv_directory/snakefile_viral_variant_visualiser``` into the path were vvv_directory was placed.  
This tool requieres snakemake 5.10 installed in a conda environment named 'snakes_update' that can be found in the ```ENVS/``` directory.  
This tool requiere the configuration of a database to determine the metagenomic composition of the sample (nt database from NCBI for exemple, or RefSeq).  
This tool requiere the implementation of a custom database made from the files of the VRL section of Genbank (viral sequences, naturally found in nature, no synthetic sequence). To do that, launch the script ```SCRIPTS_VDB/initiate_vdb.sh``` in the vvv_directory.  
