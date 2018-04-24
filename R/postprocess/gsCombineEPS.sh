indir=$1
outputFileName=$2
/oak/stanford/groups/andrewg/users/szmamie/software/ghostscript-9.23-linux-x86_64/gs-923-linux-x86_64 -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER \
    -r220 -g2660x2300 \
    -sOutputFile=$outputFileName $indir/*.eps
