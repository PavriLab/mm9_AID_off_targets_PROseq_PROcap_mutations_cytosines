#!/usr/bin/env bash
if [ -e "$2" ]; then
  if [[ $1 == 'procap' ]]; then
    for file in *.bam; do bedtools genomecov -5 -ibam "$file" -bg -strand - > "$file"_minus.bedgraph; done
    for file in *.bam; do bedtools genomecov -5 -ibam "$file" -bg -strand + > "$file"_plus.bedgraph; done
    rename 's/'.bam_minus'/_minus/g' *
    rename 's/'.bam_plus'/_plus/g' *

    for file in *minus*; do ./../TSS_scripts/TSS_analysis/sign_correction.py "$file" corrected_bg/corrected_"$file"; done
    for file in *plus*; do ./../TSS_scripts/TSS_analysis/sign_correction.py "$file" corrected_bg/corrected_"$file"; done
  elif [[ $1 == 'proseq' ]]; then
    for file in *.bam; do bedtools genomecov -3 -ibam "$file" -bg -strand - > "$file"_minus.bedgraph; done
    for file in *.bam; do bedtools genomecov -3 -ibam "$file" -bg -strand + > "$file"_plus.bedgraph; done
    rename 's/'.bam_minus'/_minus/g' *
    rename 's/'.bam_plus'/_plus/g' *

    for file in *minus*; do ./../TSS_scripts/TSS_analysis/sign_correction.py "$file" corrected_bg/corrected_"$file"; done
    for file in *plus*; do ./../TSS_scripts/TSS_analysis/sign_correction.py "$file" corrected_bg/corrected_"$file"; done
  elif [[ $1 == 'quantseq' ]]; then
    for file in *.bam; do bedtools genomecov -5 -ibam "$file" -bg -strand + > "$file"_minus.bedgraph; done
    for file in *.bam; do bedtools genomecov -5 -ibam "$file" -bg -strand - > "$file"_plus.bedgraph; done
    rename 's/'.bam_minus'/_minus/g' *
    rename 's/'.bam_plus'/_plus/g' *

    for file in *minus*; do ./../TSS_scripts/TSS_analysis/sign_correction.py "$file" corrected_bg/corrected_"$file"; done
    for file in *plus*; do ./../TSS_scripts/TSS_analysis/sign_correction.py "$file" corrected_bg/corrected_"$file"; done

  else
    echo 'Not a valid option.'
    echo 'Usage 1: ./split_and_correct_bg.sh proseq genome_file'
    echo 'Usage 2: ./split_and_correct_bg.sh procap genome_file'
    echo 'Usage 3: ./split_and_correct_bg.sh quantseq genome_file'
  fi
else
  echo 'No genome file provided'
  echo 'Usage 1: ./split_and_correct_bg.sh proseq genome_file'
  echo 'Usage 2: ./split_and_correct_bg.sh procap genome_file'
  echo 'Usage 3: ./split_and_correct_bg.sh quantseq genome_file'
fi
