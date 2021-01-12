#!/bin/bash
#by Alexander Henoch

set -e

usage () { 
	echo "Usage: $0 [-f <fasta input path>] [-c <country info input path>] [-o <output directory>] [-p <number procsses>] [-a <database with alignments (switch)>]" 1>&2; exit 1; 
}

# -------------------- variable checkup

while getopts ":f:c:o:p:a" option ; do
	case ${option} in
	f) infile=${OPTARG} ;;
	c) ininfo=${OPTARG} ;;
	o) outdir=${OPTARG} ;;
	p) proc=${OPTARG} ;;
	a) align=1 ;;
	*) usage ;;
	esac
done
shift $((OPTIND -1))

if ! [ "${outdir: -1}" = "/" ]; then
	outdir+=/
fi

if ! [ -f "$infile" ] || ! [ -f "$ininfo" ] || [ -z "$proc" ] || ! [ -d "$outdir" ]; then
	usage
fi

if ! [ -d "${outdir}.tmp" ]; then
	mkdir ${outdir}.tmp
else
	rm -r ${outdir}.tmp
	mkdir ${outdir}.tmp
fi

scriptdir=$(dirname $0)

# -------------------- temporary files

tmpcdhit=${outdir}.tmp/reprs
tmpinfo=${outdir}.tmp/info.csv
tmpgnms=${outdir}.tmp/gnms.csv
tmpcorr=${outdir}.tmp/corr.csv
tmpstrain=${outdir}.tmp/strns.tsv
tmpfsta=${outdir}.tmp/in.fasta

# -------------------- output files

outbank=${outdir}genbank/
outclstr=${outdir}clusters.tsv
outduplic=${outdir}duplicates.tsv
outgnms=${outdir}genomes.tsv
outstrain=${outdir}strains.tsv

# -------------------- preparations

echo "Working on preparations."

sed 's/ /_/g' $ininfo > $tmpinfo

if ! [ -d "$outbank" ]; then
	mkdir $outbank
	rsync -avP "ftp.ncbi.nih.gov::genbank/gbvrl*.seq.gz" $outbank
	gunzip ${outbank}gbvrl*.seq.gz
else
	echo "local genbank found, no creation necessary."
fi

awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq $0 } END { print seq }' $infile | sed 's/ /_/g;N;s/\n/|/' > $tmpgnms
python3 ${scriptdir}/wrapper.py $tmpgnms $tmpcorr $outbank 0
#####änder wrapper und readingframe auf only csv input mit comma und nicht | außerdem alles etwas bereiningen und verschönern (sauberes backbone database.sh script)
#####frequency.py ist okay
#####readingframe naja
#####wrapper zur Zeit noch Müll
#####aber alles funktioniert ohne probleme und ergebnisse sehen gut aus
#####frequency anpassen (settings) zu viele unclustered sosnt pipeline sauber
#####wrapper und readingframe als class?

echo "Preparations finished."

echo "Working on genome clustering."

python3 ${scriptdir}/frequency.py $tmpcorr $outclstr

echo "Table with clusters finished."

# -------------------- database genomes

echo "Working on genome components and strain informations."

> $tmpstrain

echo -e "accession\tstrain\tsegment\tgenome" > $outgnms
echo -e "strain\tspecies\tcity\tsubcountry\tcountry\tyear\tsubtype\thost" > $outstrain

python3 ${scriptdir}/readingframe.py $tmpinfo $tmpcorr $proc $outgnms $tmpstrain

sed -i 's/\r//' $outgnms
sed 's/\r//' $tmpstrain | sort -u >> $outstrain

echo "Table with genome components and database with strain informations finished."

# -------------------- database duplicates

sed 's/;/\n/; s/\r//' $tmpcorr > $tmpfsta
cd-hit-est -M 10000 -s 1 -c 1 -T $proc -i $tmpfsta -o $tmpcdhit

echo -e "duplicate\toriginal" > $outduplic

while read line; do
if ! [[ "$line" == \>* ]]; then
    if [[ "$line" == 0* ]]; then
        id=$(echo "$line" | sed 's/.*>\([^|]*\)|.*/\1/')
    else
        dupl=$(echo "$line" | sed 's/.*>\([^|]*\)|.*/\1/')
        echo -e "$dupl\t$id" >> $outduplic
    fi
fi
done < ${tmpcdhit}.clstr

echo "Table with duplicates finished."

echo "4 tables created. Containing duplicates, clusters, genomes and strains. Final clean up and upload to PostgreSQL database"

# -------------------- clean up

rm -r ${outdir}.tmp
