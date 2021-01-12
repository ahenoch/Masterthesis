#!/bin/bash
#by Alexander Henoch

set -e

usage () { 
	echo "Usage: $0 [-u <user>] [-d <database>] [-p <procsses>] [-f <fasta>] [-o <output directory>] [-x <optional dump name>] [-a <with alignment>]" 1>&2; exit 1; 
}

while getopts ":u:d:p:f:o:x:a" option ; do
        case ${option} in
        u) user=${OPTARG} ;;
        d) database=${OPTARG} ;;
        p) proc=${OPTARG} ;;
        f) infsta=${OPTARG} ;;
        o) outdir=${OPTARG} ;;
        x) dump=${OPTARG} ;;
        a) align=1 ;;
        *) usage ;;
        esac
done
shift $((OPTIND -1))

if ! [ "${outdir: -1}" = "/" ]; then
	outdir+=/
fi

if [ -z "$user" ] || [ -z "$database" ] || [ -z "$proc" ] || [ -z "$dump" ] || ! [ -f "$infsta" ] || ! [ -d "${outdir}" ]; then
	usage
fi

scriptdir=$(realpath $(dirname $0))

if ! [ -z "$align" ]; then
	$scriptdir/database.sh -f $infsta -c $scriptdir/Countrys.csv -o $outdir -p $proc -a
else
	$scriptdir/database.sh -f $infsta -c $scriptdir/Countrys.csv -o $outdir -p $proc
fi


psql -U $user -d $database -c "drop table if exists clusters;"

psql -U $user -d $database -c "drop table if exists duplicates;"

psql -U $user -d $database -c "drop table if exists genomes;"

psql -U $user -d $database -c "drop table if exists strains;"


psql -U $user -d $database -c "
    create table strains(
        strain text primary key not null,
        species varchar(1) not null,
        city text,
        subcountry text,
        country text,
        year integer,
        subtype text,
        host text
    );
"

psql -U $user -d $database -c "\copy strains from '${outdir}strains.tsv' delimiter E'\t' null as 'null' csv header;"


psql -U $user -d $database -c "
    create table genomes(
        accession text primary key not null,
        strain text not null references strains(strain),
        segment integer not null,
        genome text not null
    );
"

psql -U $user -d $database -c "\copy genomes from '${outdir}genomes.tsv' delimiter E'\t' null as 'null' csv header;"


psql -U $user -d $database -c "
    create table duplicates(
        accession_duplicate text primary key not null references genomes(accession),
        accession_original text not null
    );
"

psql -U $user -d $database -c "\copy duplicates from '${outdir}duplicates.tsv' delimiter E'\t' null as 'null' csv header;"


psql -U $user -d $database -c "
    create table clusters(
        accession_cluster text primary key not null references genomes(accession),
        cluster integer not null
    );
"

psql -U $user -d $database -c "\copy clusters from '${outdir}clusters.tsv' delimiter E'\t' null as 'null' csv header;"


if ! [ -z "$dump" ]; then

	pg_dump $database > ${outdir}${dump}.sql
	
fi

echo "Finished."
