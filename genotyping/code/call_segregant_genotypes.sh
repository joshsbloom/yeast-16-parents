#Restart here 10/19/17
ulimit -n 100000

# global variables 
genome='/data/rr/Parents/Reference/sacCer3.fasta'
chrs=('chrI' 'chrII' 'chrIII' 'chrIV' 'chrV' 'chrVI' 'chrVII' 'chrVIII' 'chrIX' 'chrX' 'chrXI' 'chrXII' 'chrXIII' 'chrXIV' 'chrXV' 'chrXVI')

crosses=('2999' '3000' '3001' '3049' '3003' '3004' '3043' '3028' '375' 'A'    '376'  'B'   '377'  '393'  '381'  '3008')
#note extra spaces are toxic here for 'key' values or between '=' and 'value'
declare -A RR_p1
RR_p1=( ["2999"]="YPS1009x" ["3000"]="I14a"  ["3001"]="Y10x"  ["3049"]="PW5a"  ["3003"]="273614xa" ["3004"]="YJM981x"   ["3043"]="CBS2888a"   ["3028"]="CLIB219x"   ["375"]="M22"  ["A"]="BYa"  ["376"]="RMx"        ["B"]="YPS163a"    ["377"]="YJM145x"    ["393"]="CLIB413a"   ["381"]="YJM978x"    ["3008"]="YJM454a")    
declare -A RR_p2
RR_p2=( ["2999"]="I14a"     ["3000"]="Y10x"  ["3001"]="PW5a" ["3049"]="273614xa" ["3003"]="YJM981x" ["3004"]="CBS2888a" ["3043"]="CLIB219x"   ["3028"]="M22"        ["375"]="BYa"  ["A"]="RMx"  ["376"]="YPS163a"    ["B"]="YJM145x"    ["377"]="CLIB413a"   ["393"]="YJM978x"    ["381"]="YJM454a"   ["3008"]="YPS1009x")

# split bam, goal is to remove duplicates and rerun gatk haplotypecaller
bamtools split -in /media/jbloom/d1/rr/rr16_segs.bam -tag RG -stub /data/rrv2/genotyping/segregants/ss/rr16

# remove duplicates and index bams
parallel --progress --jobs 70  'java -jar -Xms2g -Xmx2g /home/jbloom/Local/bin/picard.jar MarkDuplicates \
    I=/data/rrv2/genotyping/segregants/ss/rr16.TAG_RG_{}.bam \
    O=/data/rrv2/genotyping/segregants/segd/{}.bam \
    M=/data/rrv2/genotyping/segregants/dupstats/{}.stats \
    REMOVE_DUPLICATES=TRUE \
    CREATE_INDEX=TRUE \
    ASSUME_SORTED=TRUE' :::: /data/rrv2/genotyping/segregants/seg_lists/seg_bam_list.txt

# killed restart
#::: `tail -n +8161 /data/rrv2/genotyping/segregants/seg_lists/seg_bam_list.txt`

grep -P '\*' /data/rrv2/genotyping/parents/parents.vcf > star.variants
#cross='2999'

# extract only variants that differ for two parents for each cross into separate vcfs
for cross in "${crosses[@]}"
do
    echo $cross
    vcftools --vcf /data/rrv2/genotyping/parents/parents.vcf --indv ${RR_p1[${cross}]} --indv ${RR_p2[${cross}]} \
    --not-chr chrM \
    --non-ref-ac-any 1 \
    --max-non-ref-ac 1 \
    --exclude-positions /data/rrv2/genotyping/parents/star.variants \
    --recode-INFO-all --recode --out /data/rrv2/genotyping/parents/cross_vcf/${cross}.vcf
done

# split bams into groups by cross--------------------------------------------------
for cross in "${crosses[@]}"
do
    echo $cross
    grep $cross /data/rrv2/genotyping/segregants/seg_lists/seg_bam_list.txt | sed 's/$/.bam/g' - | sed 's#^#/data/rrv2/genotyping/segregants/segd/#g' -   > /data/rrv2/genotyping/segregants/seg_lists/${cross}_seg.list
done
#-----------------------------------------------------------------------------


#aggreagte bams and move to ssd, convert to cram
for cross in "${crosses[@]}"
do
    mf='java -jar -Xms64g -Xmx64g /home/jbloom/Local/bin/picard.jar MergeSamFiles '
    input="/data/rrv2/genotyping/segregants/seg_lists/${cross}_seg.list"
    while IFS= read -r var
    do
        s
    done < "$input"
    mf+="O=/home/jbloom/segregants/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=10000000"
    eval $mf
    (java -jar /home/jbloom/Local/cramtools/cramtools-3.0.jar cram -I /home/jbloom/segregants/${cross}.bam -R $genome -O /home/jbloom/segregants/${cross}.cram --capture-all-tags --lossy-quality-score-spec \*8;  samtools index /home/jbloom/segregants/${cross}.cram;   rm -f /home/jbloom/segregants/${cross}.bam; rm -f /home/jbloom/segregants/${cross}.bai ) &
done



# also converts to cram but compression ratio isn't as good
#samtools view -T $genome -C -@ 16 -o /home/jbloom/segregants/${cross}c.cram /home/jbloom/segregants/${cross}c.bam ;   

#Parallelize over chromosomes and crosses 
parallel  --jobs 32 'java -d64 -Xms16g -Xmx16g -jar /home/jbloom/Local/GenomeAnalysisTK-3.8/GenomeAnalysisTK.jar \
     -R /data/rr/Parents/Reference/sacCer3.fasta \
     -T HaplotypeCaller \
     -I /home/jbloom/segregants/{1}.cram \
     -o /data/rrv2/genotyping/segregants/vcfs/{1}_{2}.vcf \
     --genotyping_mode GENOTYPE_given_ALLELES \
     --alleles /data/rrv2/genotyping/parents/cross_vcf/{1}.vcf.recode.vcf \
     -L  {2} \
     --sample_ploidy 1'  ::: ${crosses[@]} ::: ${chrs[@]}

#replaces ${crosses[@]] with '3000'


# now concatenate vcf files 
cd /data/rrv2/genotyping/segregants/vcfs/
for cross in "${crosses[@]}"
do
  mf='vcf-concat '
  for it in  ${chrs[@]}
       do
          mf+="${cross}_${it}.vcf "
      done 
  mf+=" > ${cross}.vcf"
  eval $mf

    vcftools --vcf /data/rrv2/genotyping/segregants/vcfs/${cross}.vcf \
       --keep /data/rrv2/genotyping/segregants/seg_lists/${cross}_seg_name_only.list \
       --recode-INFO-all --recode --out /data/rrv2/genotyping/segregants/vcfs/${cross}.vcf

    bgzip /data/rrv2/genotyping/segregants/vcfs/${cross}.vcf.recode.vcf 
    tabix /data/rrv2/genotyping/segregants/vcfs/${cross}.vcf.recode.vcf.gz 
done
# extract only info for real segregants



# recreate parent vcfs, incorporating structural variant information
vcf-concat /data/rrv2/genotyping/parents/parents.vcf.gz /data/rrv2/genotyping/parents/structural_variants/final/final_reordered.vcf.gz > /data/rrv2/genotyping/parents/parents_w_svar.vcf
vcf-sort parents_w_svar.vcf >  parents_w_svar_sorted.vcf
bgzip  parents_w_svar_sorted.vcf
tabix  parents_w_svar_sorted.vcf.gz

for cross in "${crosses[@]}"
do
    #--exclude-positions /data/rrv2/genotyping/parents/star.variants \
    echo $cross
    vcftools --gzvcf /data/rrv2/genotyping/parents/parents_w_svar_sorted.vcf.gz --indv ${RR_p1[${cross}]} --indv ${RR_p2[${cross}]} \
    --not-chr chrM \
    --non-ref-ac-any 1 \
    --max-non-ref-ac 1 \
    --recode-INFO-all --recode --out /data/rrv2/genotyping/parents/cross_vcf/${cross}_w_svar.vcf
    bgzip /data/rrv2/genotyping/parents/cross_vcf/${cross}_w_svar.vcf.recode.vcf 
    tabix /data/rrv2/genotyping/parents/cross_vcf/${cross}_w_svar.vcf.recode.vcf.gz

done


# Run haplotype caller. Is further parallelization necessary?
# this annoyingly outputs blanks for all of the segregants not part of cross of interest
parallel  'java -d64 -Xms24g -Xmx24g -jar /home/jbloom/Local/GenomeAnalysisTK-3.8/GenomeAnalysisTK.jar \
     -R /data/rr/Parents/Reference/sacCer3.fasta \
     -T HaplotypeCaller \
     -I /data/rrv2/genotyping/segregants/seg_lists/{}_seg.list \
     -o /data/rrv2/genotyping/segregants/vcfs/{}.vcf \
     --genotyping_mode GENOTYPE_given_ALLELES \
     --alleles /data/rrv2/genotyping/parents/cross_vcf/{}.vcf.recode.vcf \
     -L  /data/rrv2/genotyping/parents/cross_vcf/{}.vcf.recode.vcf \
      -ip 100 \
      -nct 2 \
     --sample_ploidy 1'  ::: ${crosses[@]}

#not sure why 3049 failed 
# try running manually
parallel  'java -d64 -Xms24g -Xmx24g -jar /home/jbloom/Local/GenomeAnalysisTK-3.8/GenomeAnalysisTK.jar \
     -R /data/rr/Parents/Reference/sacCer3.fasta \
     -T HaplotypeCaller \
     -I /data/rrv2/genotyping/segregants/seg_lists/{}_seg.list \
     -o /data/rrv2/genotyping/segregants/vcfs/{}.vcf \
     --genotyping_mode GENOTYPE_given_ALLELES \
     --alleles /data/rrv2/genotyping/parents/cross_vcf/{}.vcf.recode.vcf \
     --sample_ploidy 1'  ::: '3049' 
#${crosses[@]}

#sed 's/.*\/\(.*\).bam/\1/g' /data/rrv2/genotyping/segregants/seg_lists/2999_seg.list
for cross in "${crosses[@]}"
do
    echo $cross
    grep $cross /data/rrv2/genotyping/segregants/seg_lists/seg_bam_list.txt > /data/rrv2/genotyping/segregants/seg_lists/${cross}_seg_name_only.list
done


#vcftools --vcf /data/rrv2/genotyping/segregants/vcfs/2999.vcf \
#    --keep /data/rrv2/genotyping/segregants/seg_lists/2999_seg_name_only.list \
#    --recode-INFO-all --recode --out /data/rrv2/genotyping/segregants/vcfs/2999.vcf
#
#"${crosses[@]}" 
#`ls /data/rrv2/genotyping/parents/Ldir`


