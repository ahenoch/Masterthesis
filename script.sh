#!/bin/bash

set -e

usage () { 
	echo "Usage: $0 [-u <user>] [-d <database>] [-p <number procsses>]" 1>&2; exit 1; 
}

while getopts ":u:d:p:" option ; do
        case ${option} in
        u) user=${OPTARG} ;;
        d) database=${OPTARG} ;;
        p) proc=${OPTARG} ;;
        *) usage ;;
        esac
done
shift $((OPTIND -1))

if [ -z "$proc" ] || [ -z "$user" ] || [ -z "$database" ]; then
	usage
fi

scriptdir=$(realpath $(dirname $0))

if ! [ -d "$scriptdir/Pipeline" ]; then
	mkdir $scriptdir/Pipeline
	mkdir $scriptdir/Pipeline/{ALN,FASTA,LRIscan}
else
	rm -r $scriptdir/Pipeline
	mkdir $scriptdir/Pipeline
	mkdir $scriptdir/Pipeline/{ALN,FASTA,LRIscan}
fi	

conda activate masterthesis

echo 'Working on creating FASTAs.'

for i in $(seq 1 8); do 
	for j in $(seq $(( $i+1 )) 8); do

		$scriptdir/Scripts/retrieve.sh -u $user -d $database -o $scriptdir/Pipeline/FASTA -g $i -g $j -s B -r 120 -x seg_${i}_${j}_neg.fasta -n
		
python3 << END
from Bio import SeqIO
records = SeqIO.parse("$scriptdir/Pipeline/FASTA/seg_${i}_${j}_neg.fasta", "fasta")
count = SeqIO.write(records, "$scriptdir/Pipeline/ALN/seg_${i}_${j}_neg.aln", "clustal")
END

	done 
done

echo 'Finished creating FASTAs.'

conda deactivate
