cd ~/vdb/viral_genbank_files.ori/
pwd
for i in `seq 0 9`
do
    for a in `seq 0 9`
    do
        echo $i$a

        chmod 777  *$i$a.*.gbk
        rm *$i$a.*.gbk
    done
done
