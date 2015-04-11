#!/bin/sh

curl --proxy proxy.mgmt.hpc.mssm.edu:8123 --data-urlencode sequence@$1_curated.fa http://www.ncbi.nlm.nih.gov/genomes/FLU/Database/annotation.cgi > $1_temp.txt

cat $1_temp.txt | awk "/id='tbl'/ {p=1}; p; /id='xml'/ {p=0}" > $1_feature_temp.txt

#parse the first line, starting with &gt;Feature and add the rest of the line to the temporary file
#add the rest of the matching lines, from line 2 to the end and then remove the last line (which contains the xml div after the table we're interestedin)
#finally replace the url encoded "&gt;" with ">" and save to the final .anno file
a=$(head -1 $1_feature_temp.txt)
b=$(echo $a | awk -F"<pre>&gt;" '{print ">"$2}')
echo $b > $1.anno
tail -n +2 $1_feature_temp.txt | head -n -1 >> $1.anno 
sed -i 's/&gt;/>/g' $1.anno

cat $1_temp.txt | awk "/ffn_h/ {p=1}; p; /aln_h/ {p=0}" > $1_cds_temp.txt
a=$(head -1 $1_cds_temp.txt)
b=$(echo $a | awk -F"&gt;" '{print ">"$2}')
echo $b > $1_cds.fa
tail -n +2 $1_cds_temp.txt | head -n -1 >> $1_cds.fa
sed -i 's/&gt;/>/g' $1_cds.fa

rm -f $1_*temp.txt
