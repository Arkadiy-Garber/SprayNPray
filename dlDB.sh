if [ "$#" == 0 ] || [ $1 == "-h" ]; then
  printf "Usage:\t dlDB.sh location_for_db \n"
  exit
fi

wget --directory-prefix=$1 ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
