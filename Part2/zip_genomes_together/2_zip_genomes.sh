#!/usr/bin/bash

set -e

cd ./mid/1_grouped_genomes
for SAG_DIR in */
do
	ID=${SAG_DIR%/}
	DST="/mnt/scgc_nfs/lab/ggavelis/analysis/siavash_hgt/aug9_second_round/mid/2_gzipped_inputs/${ID}.tar.gz"
	#DST="../mid/2_zipped_inputs/${ID}.tar.gz"
	#if [ ! -f /mnt/scgc_nfs/data/deliverables/230912_VH00511_108_AAC777CHV/klewis/${ID}.tar.gz ] | [ "$OVERWRITE_TAR" = TRUE ]; then echo "zipping ${ID}"; tar cf - ${ID} | pigz --fast -p 8 > /mnt/scgc_nfs/data/deliverables/230912_VH00511_108_AAC777CHV/klewis/${ID}.tar.gz; fi; done
	if [ ! -f $DST ]
	then
		echo "zipping ${ID}"
		tar cf - ${ID} | pigz --fast -p 8 > $DST
	fi
done
