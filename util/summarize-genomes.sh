usage (){
    echo "Generate summary files for a set of genomes"
    echo "REQUIRED ARGUMENTS"
    echo "  -i Input data directory"
    echo "  -o Output directory"
    exit 0
}

# print help with no arguments
[[ $# -eq 0 ]] && usage

while getopts "hi:o:" opt; do
    case $opt in
        h)
            usage ;;
        i) 
            idir=${OPTARG%/} ;;
        o) 
            odir=${OPTARG%/} ;;
    esac 
done

[[ -d $odir ]] || mkdir $odir

# input: set of full genomes
scaflen=$odir/scaffold-lengths.tab
ls $idir/*fna | parallel "smof stat -q {} > $odir/{/}.tab "
for j in $odir/*fna.tab
do
    s=${j%%.fna.tab}
    s=`basename $s`
    sed "s/^/$s\t/" $j
    rm $j
done > $scaflen
