

Given two GFA, merge

if [ $# -lt 3 ]; 
then
	echo "sh merge_gfa_wkflow.sh <GFA> <GFA>"
else
	fa_file=$2
	fa_file2=$3
	minimap2 -c fa_file fa_file2 | grep "tp:A:p" | sort -k6 > sorted_primary_alignment.paf
	./merge_gfa 



fi


