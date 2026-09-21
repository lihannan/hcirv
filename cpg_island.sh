#!/bin/bash
# CpG island detection using EMBOSS cpgplot
# Dependency: EMBOSS cpgplot, GNU parallel

set -e

ls *fasta | awk -F "." '{print $1}' | parallel -j 1 'cpgplot -sequence {}.fasta -window 50 -minlen 80 -outfile {}.cpgplot -minoe 0.6 -minpc 60 -graph png -gtitle "{}" -outfeat {}.gff -goutfile {}.png'

echo "Done."
