#! /usr/bin/env bash
# launch this script in the vvv_directory/SCRIPTS_VDB/ directory

mkdir ../VDB/megablast_index/
mkdir ../VDB/viral_gbk_files/ && cd ../VDB/viral_gbk_files

# update the refseq in genbank format
mkdir refseq && cd refseq
wget http://ftp.ncbi.nih.gov/genomes/Viruses/all.gbk.tar.gz
gunzip all.gbk.tar.gz
tar -xf all.gbk.tar

mv */*.gbk . 

# rename refseq file, in order the name comes with version
for i in `ls | grep -P '\.gbk'`
do
    mv $i `grep VERSION $i | cut -f 6 --delimiter=" "`.gbk
done

mv *.gbk ../ && cd ../ && rm -r refseq/

## update the genbank databse
# download the genbank files from ncbi ftp server
# these files are multi gb files from the GENBANK database

for i in `seq 1 40`
# at some point, the 40 should be modified, but, right now (2019), we have time
do
    wget http://ftp.ncbi.nih.gov/genbank/gbvrl$i'.seq.gz'
done

# unzip the genbank files
gunzip gbvrl*.seq.gz
cd ../viral_gbk_files

# these files are multigenbank files, it is necessary to split them all
for i in `seq 1 40`
do
    python3 ../../SCRIPTS_VDB/split_multi_gbk_file.py gbvrl$i'.seq'
done

# extract all fasta seq from genbank files
ls | grep -P '\.gbk' > all_file_names
python3 ../../SCRIPTS_VDB/extract_fasta_from_gb.py all_file_names all_seq.fasta && rm all_file_names

cd ../

mv viral_gbk_files/all_seq.fasta megablast_index/vdb.fasta && cd megablast_index/

# to index database 
source activate megablast
makeblastdb -in vdb.fasta -dbtype nucl
conda deactivate

grep ">" vdb.fasta | sed 's/>//g' > vdb_ids


echo The Viral DataBase has been successfully updated.
echo for regular megablast, i.e. use the following command
echo blastn -task megablast -db path_to_vdb/megablast_index/vdb.fasta -query input.fasta -dbtype nucl -o output.megablast

