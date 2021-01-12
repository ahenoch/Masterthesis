#!/bin/bash
#by Alexander Henoch

set -e

usage () { 
	#echo "Usage: $0 [-u <user name>] [-d <database name>] [-o <output directory>] [-g <segment filter>] [-s <species filter>] [-r <primer>] [-x <output fasta name>] [-f <include five prime UTR>] [-p <include ORF>] [-t <include three prime UTR>] [-n <include negative strain>]" 1>&2; exit 1; 
	echo -e "Usage: $0 [Options]\n\nOptions\n\n   -u\texisting PostgreSQL username\n   -d\texisting PostgreSQL ISU database\n   -o\tpath to output directory\n   -g\tfilter for segments\t\t\t(e.g. -g 1 -g 2)\n   -s\tfilter for specific species\t\t(e.g. -s C)\n   -r\tcustom primer sequence between segments\n   -x\toutput fasta name\n   -f\tinclude segments five prime UTR\t\t(switch)\n   -p\tinclude segments ORF\t\t\t(switch)\n   -t\tinclude segments three prime UTR\t(switch)\n   -n\tinclude segments negative strain\t(switch)\n\n   Script for ISU database search and subsequent FASTA conversion.\n\n   Made by Alexander Henoch." 1>&2; exit 1; 
}

while getopts ":u:d:o:g:s:r:x:fptn" option ; do
        case ${option} in
        u) user=${OPTARG} ;;
        d) database=${OPTARG} ;;
        o) outdir=${OPTARG} ;;
        g) genomes=( ${genomes[@]} ${OPTARG} ) ;;
        s) species=( ${species[@]} ${OPTARG} ) ;;
        r) primer=${OPTARG} ;;
        x) name=${OPTARG} ;;
        f) five=1 ;;
        p) positive=1 ;;
        t) three=1 ;;
        n) negative=1 ;;
        *) usage ;;
        esac
done
shift $((OPTIND -1))

if ! [ "${outdir: -1}" = "/" ]; then
	outdir+=/
fi

#Input format check

if [ -z "$user" ] || [ -z "$database" ] || [ -z "$name" ] || [ -z "$primer" ] || [ -z "$genomes" ] || [ -z "$species" ] || ! [ -d "$outdir" ]; then
	echo -e "Some input value missing or output directory not existing.\n"
	usage
fi

#if ! [[ $random =~ ^[0-9]+$ ]]; then 
#	echo "Invalid number of random nucleotides."
#	usage
#fi

for g in ${genomes[@]}; do
	if ! [[ $g =~ ^[1-8]+$ ]]; then 
		echo -e "Invalid segment names [1-8].\n"
		usage
	fi
done

for s in ${species[@]}; do
	if ! [[ $s =~ ^[A-C]+$ ]]; then 
		echo -e "Invalid species names [A-C].\n"
		usage
	fi
done

#Finished format check

#sel="select '>gb:'||v.id||'|'"
sel="select distinct '>'"
for g in ${genomes[@]}; do
	sel+="||'Segment '||coalesce(g$g.segment, '8')||':'||coalesce(g$g.accession, 'none')||'|'"
done

sel+="||'Strain:'||s${genomes[0]}.strain||'|'||'Species:'||s${genomes[0]}.species||'|'||'Subtype:'||coalesce(s${genomes[0]}.subtype, 'none')||';'"
for g in ${genomes[@]}; do
	if [ "$g" = "${genomes[-1]}" ]; then
		if ! [ -z "$five" ]; then
			sel+="||coalesce(g$g.utr_five, '')"
		fi
		if ! [ -z "$positive" ]; then
			sel+="||coalesce(g$g.orf, '')"
		fi
		if ! [ -z "$three" ]; then
			sel+="||coalesce(g$g.utr_three, '')"
		fi
		if ! [ -z "$negative" ]; then
			sel+="||coalesce(g$g.genome, '')"
		fi
	else
		if ! [ -z "$five" ]; then
			sel+="||coalesce(g$g.utr_five, '')"
		fi
		if ! [ -z "$positive" ]; then
			sel+="||coalesce(g$g.orf, '')"
		fi
		if ! [ -z "$three" ]; then
			sel+="||coalesce(g$g.utr_three, '')"
		fi
		if ! [ -z "$negative" ]; then
			sel+="||coalesce(g$g.genome, '')"
		fi
		sel+="||'${primer}'"
	fi
done

join=""
for g in ${genomes[@]}; do
	if [ "$g" = "${genomes[0]}" ]; then
		join+="left outer join genomes g$g on v.segment_$g = g$g.accession left outer join strains s$g on g$g.strain = s$g.strain"
	else
		join+=" left outer join genomes g$g on v.segment_$g = g$g.accession"
	fi
done

where="where"
for s in ${species[@]}; do
	if [ "$s" = "${species[-1]}" ]; then
		where+=" s${genomes[0]}.species = '$s'"	
	else
		where+=" s${genomes[0]}.species = '$s' or"
	fi
done

psql -q -d $database -U $user -c "drop view if exists fasta;"

psql -q -d $database -U $user -c "
	create view fasta as
	$sel
	from variations v
	$join
	$where
;"

psql -q -d $database -U $user -t -c "\copy (select * from fasta) to '${outdir}tmp' csv"

psql -q -d $database -U $user -c "drop view if exists fasta;"

sed -i 's/;/\n/g' ${outdir}tmp

awk '/^>/{sub(">", ">gb:"++i"|")}1' ${outdir}tmp > ${outdir}${name}

rm ${outdir}tmp
