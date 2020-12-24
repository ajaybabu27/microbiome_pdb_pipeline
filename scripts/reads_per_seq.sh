#!/usr/bin/env bash
# reads_per_seq
#
# Run 8 jobs in parallel on many BAM files:
#
#   parallel -j8 reads_per_seq ::: *.bam
#
# Bit   Description
# ===   ===========
# 0x1   template having multiple fragments in sequencing
# 0x2   each fragment properly aligned according to the aligner
# 0x4   fragment unmapped
# 0x8   next fragment in the template unmapped
# 0x10  SEQ being reverse complemented
# 0x20  SEQ of the next fragment in the template being reversed 0x40 the first fragment in the template
# 0x80  the last fragment in the template
# 0x100 secondary alignment
# 0x200 not passing quality controls 0x400 PCR or optical duplicate

set -eou pipefail

if [[ $# -ne 1 || "$1" != *.bam ]]; then
  echo "usage: $(basename $0) file.bam"
  echo
  echo "  Produce a tab-delimited file.seq with number per sequence of"
  echo "  primary alignments with proper read pairs."
  echo
  exit 1
fi

bam="$1"

samtools view -f 2 -F 104 "$bam" \
  | cut -f3 | sort | uniq -c | sort -k1nr \
  | awk 'OFS="\t" {print $1,$2}' \
> "${bam%.bam}.seq"
