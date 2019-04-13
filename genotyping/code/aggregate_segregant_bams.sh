# bash 4 associative array
crosses=('2999' '3000' '3001' '3049' '3003' '3004' '3043' '3028' '375' 'A'    '376'  'B'   '377'  '393'  '381'  '3008')
#note extra spaces are toxic here for 'key' values or between '=' and 'value'
declare -A RR_p1
RR_p1=( ["2999"]="YPS1009x" ["3000"]="I14a"  ["3001"]="Y10x"  ["3049"]="PW5a"  ["3003"]="273614xa" ["3004"]="YJM981x"   ["3043"]="CBS2888a"   ["3028"]="CLIB219x"   ["375"]="M22"  ["A"]="BYa"  ["376"]="RMx"        ["B"]="YPS163a"    ["377"]="YJM145x"    ["393"]="CLIB413a"   ["381"]="YJM978x"    ["3008"]="YJM454a")    
declare -A RR_p2
RR_p2=( ["2999"]="I14a"     ["3000"]="Y10x"  ["3001"]="PW5a" ["3049"]="273614xa" ["3003"]="YJM981x" ["3004"]="CBS2888a" ["3043"]="CLIB219x"   ["3028"]="M22"        ["375"]="BYa"  ["A"]="RMx"  ["376"]="YPS163a"    ["B"]="YJM145x"    ["377"]="CLIB413a"   ["393"]="YJM978x"    ["381"]="YJM454a"   ["3008"]="YPS1009x")

# Now for each block of segregants, merge bams
##################################################################################################################################
cross='2999'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
picard ReorderSam I=/data/rr/02/N702/bam/2999_G1.bam O=${crossdir}/2999_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/04/N701/bam/2999_G2.bam O=${crossdir}/2999_G2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/04/N702/bam/2999_G3.bam O=${crossdir}/2999_G3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/04/N703/bam/2999_G4.bam O=${crossdir}/2999_G4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/04/N704/bam/2999_G5.bam O=${crossdir}/2999_G5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/05/N701/bam/2999_R1.bam O=${crossdir}/2999_R1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/09/N701/bam/2999_R2.bam O=${crossdir}/2999_R2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/09/N702/bam/2999_R3.bam O=${crossdir}/2999_R3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/09/N703/bam/2999_R4.bam O=${crossdir}/2999_R4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/09/N704/bam/2999_R5.bam O=${crossdir}/2999_R5.bam CREATE_INDEX=TRUE REFERENCE=$genome &

picard MergeSamFiles I=${crossdir}/2999_G1.bam \
                     I=${crossdir}/2999_G2.bam \
                     I=${crossdir}/2999_G3.bam \
                     I=${crossdir}/2999_G4.bam \
                     I=${crossdir}/2999_G5.bam \
                     I=${crossdir}/2999_R1.bam \
                     I=${crossdir}/2999_R2.bam \
                     I=${crossdir}/2999_R3.bam \
                     I=${crossdir}/2999_R4.bam \
                     I=${crossdir}/2999_R5.bam \
                     O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000
rm -f ${crossdir}/2999_G1.bam  ${crossdir}/2999_G2.bam ${crossdir}/2999_G3.bam  ${crossdir}/2999_G4.bam  ${crossdir}/2999_G5.bam  ${crossdir}/2999_R1.bam  ${crossdir}/2999_R2.bam  ${crossdir}/2999_R3.bam  ${crossdir}/2999_R4.bam   ${crossdir}/2999_R5.bam 


##################################################################################################################################
cross='3000'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
picard ReorderSam I=/data/rr/07/N711/bam/3000_G1.bam O=${crossdir}/3000_G1.bam  CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/11/N701/bam/3000_G2.bam O=${crossdir}/3000_G2.bam  CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/11/N702/bam/3000_G3.bam O=${crossdir}/3000_G3.bam  CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/11/N703/bam/3000_G4.bam O=${crossdir}/3000_G4.bam  CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/12/N701/bam/3000_G5.bam O=${crossdir}/3000_G5.bam  CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/10/N701/bam/3000_R1.bam O=${crossdir}/3000_R1.bam  CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/12/N702/bam/3000_R2.bam O=${crossdir}/3000_R2.bam  CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/12/N703/bam/3000_R3.bam O=${crossdir}/3000_R3.bam  CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/12/N704/bam/3000_R4.bam O=${crossdir}/3000_R4.bam  CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/12/N705/bam/3000_R5.bam O=${crossdir}/3000_R5.bam  CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/14/N711/bam/3000_G1.bam O=${crossdir}/3000_G1b.bam CREATE_INDEX=TRUE REFERENCE=$genome &

picard MergeSamFiles I=${crossdir}/3000_G1.bam  \
                     I=${crossdir}/3000_G2.bam  \
                     I=${crossdir}/3000_G3.bam  \
                     I=${crossdir}/3000_G4.bam  \
                     I=${crossdir}/3000_G5.bam  \
                     I=${crossdir}/3000_R1.bam  \
                     I=${crossdir}/3000_R2.bam  \
                     I=${crossdir}/3000_R3.bam  \
                     I=${crossdir}/3000_R4.bam  \
                     I=${crossdir}/3000_R5.bam  \
                     I=${crossdir}/3000_G1b.bam \
                     O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000
rm -f ${crossdir}/3000_G1.bam ${crossdir}/3000_G2.bam ${crossdir}/3000_G3.bam ${crossdir}/3000_G4.bam ${crossdir}/3000_G5.bam ${crossdir}/3000_R1.bam ${crossdir}/3000_R2.bam ${crossdir}/3000_R3.bam ${crossdir}/3000_R4.bam ${crossdir}/3000_R5.bam ${crossdir}/3000_G1b.bam 
##################################################################################################################################
cross='3001'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
picard ReorderSam I=/data/rr/02/N704/bam/3001_G1.bam O=${crossdir}/3001_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome & 
picard ReorderSam I=/data/rr/03/N701/bam/3001_G2.bam O=${crossdir}/3001_G2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/03/N702/bam/3001_G3.bam O=${crossdir}/3001_G3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/03/N703/bam/3001_G4.bam O=${crossdir}/3001_G4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/03/N704/bam/3001_G5.bam O=${crossdir}/3001_G5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/03/N707/bam/3001_R1.bam O=${crossdir}/3001_R1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/05/N702/bam/3001_R2.bam O=${crossdir}/3001_R2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/05/N703/bam/3001_R3.bam O=${crossdir}/3001_R3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/05/N704/bam/3001_R4.bam O=${crossdir}/3001_R4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/05/N705/bam/3001_R5.bam O=${crossdir}/3001_R5.bam CREATE_INDEX=TRUE REFERENCE=$genome &

picard MergeSamFiles \
I=${crossdir}/3001_G1.bam \
I=${crossdir}/3001_G2.bam \
I=${crossdir}/3001_G3.bam \
I=${crossdir}/3001_G4.bam \
I=${crossdir}/3001_G5.bam \
I=${crossdir}/3001_R1.bam \
I=${crossdir}/3001_R2.bam \
I=${crossdir}/3001_R3.bam \
I=${crossdir}/3001_R4.bam \
I=${crossdir}/3001_R5.bam \
O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000
rm -f ${crossdir}/3001_G1.bam  ${crossdir}/3001_G2.bam  ${crossdir}/3001_G3.bam  ${crossdir}/3001_G4.bam  ${crossdir}/3001_G5.bam  ${crossdir}/3001_R1.bam  ${crossdir}/3001_R2.bam  ${crossdir}/3001_R3.bam  ${crossdir}/3001_R4.bam  ${crossdir}/3001_R5.bam  

gatk -T DepthOfCoverage \
-R $genome  \
-I /data/rr/Segregants/3001/3001.bam  \
-o /data/rr/Segregants/3001/3001.coverage  \
-l chrM \
--minMappingQuality 1 
##################################################################################################################################
cross='3003'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
picard ReorderSam I=/data/rr/02/N705/bam/3003_G1.bam O=${crossdir}/3003_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/06/N701/bam/3003_G2.bam O=${crossdir}/3003_G2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/06/N702/bam/3003_G3.bam O=${crossdir}/3003_G3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/06/N703/bam/3003_G4.bam O=${crossdir}/3003_G4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/06/N704/bam/3003_G5.bam O=${crossdir}/3003_G5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/06/N705/bam/3003_R1.bam O=${crossdir}/3003_R1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/09/N705/bam/3003_R2.bam O=${crossdir}/3003_R2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/10/N702/bam/3003_R3.bam O=${crossdir}/3003_R3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/10/N703/bam/3003_R4.bam O=${crossdir}/3003_R4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/10/N704/bam/3003_R5.bam O=${crossdir}/3003_R5.bam CREATE_INDEX=TRUE REFERENCE=$genome &

picard MergeSamFiles \
I=${crossdir}/3003_G1.bam \
I=${crossdir}/3003_G2.bam \
I=${crossdir}/3003_G3.bam \
I=${crossdir}/3003_G4.bam \
I=${crossdir}/3003_G5.bam \
I=${crossdir}/3003_R1.bam \
I=${crossdir}/3003_R2.bam \
I=${crossdir}/3003_R3.bam \
I=${crossdir}/3003_R4.bam \
I=${crossdir}/3003_R5.bam \
O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000

rm -f ${crossdir}/3003_G1.bam ${crossdir}/3003_G2.bam ${crossdir}/3003_G3.bam ${crossdir}/3003_G4.bam ${crossdir}/3003_G5.bam ${crossdir}/3003_R1.bam ${crossdir}/3003_R2.bam ${crossdir}/3003_R3.bam ${crossdir}/3003_R4.bam ${crossdir}/3003_R5.bam 
##################################################################################################################################
cross='3004'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
#picard ReorderSam I=/data/rr/02/N706/bam/3004_G1.bam O=${crossdir}/3004_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
#    #picard ReorderSam I=/data/rr/14/N707/bam/3004_G1.bam O=${crossdir}/3004_G1b.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/15/bam/3004_G1.bam O=${crossdir}/3004_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome
picard ReorderSam I=/data/rr/08/N708/bam/3004_G2.bam O=${crossdir}/3004_G2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/08/N709/bam/3004_G3.bam O=${crossdir}/3004_G3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/08/N710/bam/3004_G4.bam O=${crossdir}/3004_G4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/07/N712/bam/3004_G5.bam O=${crossdir}/3004_G5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/06/N706/bam/3004_R1.bam O=${crossdir}/3004_R1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/10/N709/bam/3004_R2.bam O=${crossdir}/3004_R2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/10/N710/bam/3004_R3.bam O=${crossdir}/3004_R3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/10/N711/bam/3004_R4.bam O=${crossdir}/3004_R4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/10/N712/bam/3004_R5.bam O=${crossdir}/3004_R5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/14/N712/bam/3004_G5.bam O=${crossdir}/3004_G5b.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/14/N706/bam/3004_G1.bam O=${crossdir}/3004_G1b.bam CREATE_INDEX=TRUE REFERENCE=$genome &

#I=${crossdir}/3004_G1.bam  \
#I=${crossdir}/3004_G1b.bam \

picard MergeSamFiles \
I=${crossdir}/3004_G1.bam  \
I=${crossdir}/3004_G2.bam  \
I=${crossdir}/3004_G3.bam  \
I=${crossdir}/3004_G4.bam  \
I=${crossdir}/3004_G5.bam  \
I=${crossdir}/3004_R1.bam  \
I=${crossdir}/3004_R2.bam  \
I=${crossdir}/3004_R3.bam  \
I=${crossdir}/3004_R4.bam  \
I=${crossdir}/3004_R5.bam  \
I=${crossdir}/3004_G5b.bam \
O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000
rm -f ${crossdir}/3004_G1.bam  ${crossdir}/3004_G2.bam  ${crossdir}/3004_G3.bam  ${crossdir}/3004_G4.bam  ${crossdir}/3004_G5.bam  ${crossdir}/3004_R1.bam  ${crossdir}/3004_R2.bam  ${crossdir}/3004_R3.bam  ${crossdir}/3004_R4.bam  ${crossdir}/3004_R5.bam  ${crossdir}/3004_G5b.bam ${crossdir}/3004_G1b.bam 
##################################################################################################################################
cross='3008'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
#picard ReorderSam I=/data/rr/02/N707/bam/3008_G1.bam O=${crossdir}/3008_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/15/bam/3008_G1.bam O=${crossdir}/3008_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome
picard ReorderSam I=/data/rr/04/N705/bam/3008_G2.bam O=${crossdir}/3008_G2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/04/N706/bam/3008_G3.bam O=${crossdir}/3008_G3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/04/N707/bam/3008_G4.bam O=${crossdir}/3008_G4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/04/N708/bam/3008_G5.bam O=${crossdir}/3008_G5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/05/N706/bam/3008_R1.bam O=${crossdir}/3008_R1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/08/N711/bam/3008_R2.bam O=${crossdir}/3008_R2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/08/N712/bam/3008_R3.bam O=${crossdir}/3008_R3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/09/N706/bam/3008_R4.bam O=${crossdir}/3008_R4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/09/N707/bam/3008_R5.bam O=${crossdir}/3008_R5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/14/N707/bam/3008_G1.bam O=${crossdir}/3008_G1b.bam CREATE_INDEX=TRUE REFERENCE=$genome &

#I=${crossdir}/3008_G1b.bam \

picard MergeSamFiles \
I=${crossdir}/3008_G1.bam  \
I=${crossdir}/3008_G2.bam  \
I=${crossdir}/3008_G3.bam  \
I=${crossdir}/3008_G4.bam  \
I=${crossdir}/3008_G5.bam  \
I=${crossdir}/3008_R1.bam  \
I=${crossdir}/3008_R2.bam  \
I=${crossdir}/3008_R3.bam  \
I=${crossdir}/3008_R4.bam  \
I=${crossdir}/3008_R5.bam  \
O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000

##################################################################################################################################
cross='3028'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
picard ReorderSam I=/data/rr/02/N708/bam/3028_G1.bam O=${crossdir}/3028_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/04/N709/bam/3028_G2.bam O=${crossdir}/3028_G2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/04/N710/bam/3028_G3.bam O=${crossdir}/3028_G3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/04/N711/bam/3028_G4.bam O=${crossdir}/3028_G4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/04/N712/bam/3028_G5.bam O=${crossdir}/3028_G5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/03/N709/bam/3028_R1.bam O=${crossdir}/3028_R1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/06/N707/bam/3028_R2.bam O=${crossdir}/3028_R2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/06/N708/bam/3028_R3.bam O=${crossdir}/3028_R3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/06/N709/bam/3028_R4.bam O=${crossdir}/3028_R4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/06/N710/bam/3028_R5.bam O=${crossdir}/3028_R5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/14/N708/bam/3028_G1.bam O=${crossdir}/3028_G1b.bam CREATE_INDEX=TRUE REFERENCE=$genome &

picard MergeSamFiles \
I=${crossdir}/3028_G1.bam  \
I=${crossdir}/3028_G2.bam  \
I=${crossdir}/3028_G3.bam  \
I=${crossdir}/3028_G4.bam  \
I=${crossdir}/3028_G5.bam  \
I=${crossdir}/3028_R1.bam  \
I=${crossdir}/3028_R2.bam  \
I=${crossdir}/3028_R3.bam  \
I=${crossdir}/3028_R4.bam  \
I=${crossdir}/3028_R5.bam  \
I=${crossdir}/3028_G1b.bam \
O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000
##################################################################################################################################

cross='3043'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
picard ReorderSam I=/data/rr/11/N711/bam/3043_G1.bam O=${crossdir}/3043_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/12/N706/bam/3043_G2.bam O=${crossdir}/3043_G2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/12/N707/bam/3043_G3.bam O=${crossdir}/3043_G3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/12/N708/bam/3043_G4.bam O=${crossdir}/3043_G4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/12/N709/bam/3043_G5.bam O=${crossdir}/3043_G5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/11/N712/bam/3043_R1.bam O=${crossdir}/3043_R1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/12/N710/bam/3043_R2.bam O=${crossdir}/3043_R2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/12/N711/bam/3043_R3.bam O=${crossdir}/3043_R3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/12/N712/bam/3043_R4.bam O=${crossdir}/3043_R4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/13/N701/bam/3043_R5.bam O=${crossdir}/3043_R5.bam CREATE_INDEX=TRUE REFERENCE=$genome &

picard MergeSamFiles \
I=${crossdir}/3043_G1.bam  \
I=${crossdir}/3043_G2.bam  \
I=${crossdir}/3043_G3.bam  \
I=${crossdir}/3043_G4.bam  \
I=${crossdir}/3043_G5.bam  \
I=${crossdir}/3043_R1.bam  \
I=${crossdir}/3043_R2.bam  \
I=${crossdir}/3043_R3.bam  \
I=${crossdir}/3043_R4.bam  \
I=${crossdir}/3043_R5.bam  \
O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000
##################################################################################################################################
cross='3049'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
picard ReorderSam I=/data/rr/02/N710/bam/3049_G1.bam O=${crossdir}/3049_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/03/N706/bam/3049_G2.bam O=${crossdir}/3049_G2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/05/N707/bam/3049_G3.bam O=${crossdir}/3049_G3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/05/N708/bam/3049_G4.bam O=${crossdir}/3049_G4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/05/N709/bam/3049_G5.bam O=${crossdir}/3049_G5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/03/N710/bam/3049_R1.bam O=${crossdir}/3049_R1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/05/N710/bam/3049_R2.bam O=${crossdir}/3049_R2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/05/N711/bam/3049_R3.bam O=${crossdir}/3049_R3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/05/N712/bam/3049_R4.bam O=${crossdir}/3049_R4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/06/N711/bam/3049_R5.bam O=${crossdir}/3049_R5.bam CREATE_INDEX=TRUE REFERENCE=$genome &

picard MergeSamFiles \
I=${crossdir}/3049_G1.bam  \
I=${crossdir}/3049_G2.bam  \
I=${crossdir}/3049_G3.bam  \
I=${crossdir}/3049_G4.bam  \
I=${crossdir}/3049_G5.bam  \
I=${crossdir}/3049_R1.bam  \
I=${crossdir}/3049_R2.bam  \
I=${crossdir}/3049_R3.bam  \
I=${crossdir}/3049_R4.bam  \
I=${crossdir}/3049_R5.bam  \
O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000
##################################################################################################################################
cross='375'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
picard ReorderSam I=/data/rr/02/N711/bam/375_G1.bam O=${crossdir}/375_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/07/N701/bam/375_G2.bam O=${crossdir}/375_G2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/07/N702/bam/375_G3.bam O=${crossdir}/375_G3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/07/N703/bam/375_G4.bam O=${crossdir}/375_G4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/07/N704/bam/375_G5.bam O=${crossdir}/375_G5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/03/N711/bam/375_R1.bam O=${crossdir}/375_R1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/07/N705/bam/375_R2.bam O=${crossdir}/375_R2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/07/N706/bam/375_R3.bam O=${crossdir}/375_R3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/07/N707/bam/375_R4.bam O=${crossdir}/375_R4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/07/N708/bam/375_R5.bam O=${crossdir}/375_R5.bam CREATE_INDEX=TRUE REFERENCE=$genome &

picard MergeSamFiles \
I=${crossdir}/375_G1.bam \
I=${crossdir}/375_G2.bam \
I=${crossdir}/375_G3.bam \
I=${crossdir}/375_G4.bam \
I=${crossdir}/375_G5.bam \
I=${crossdir}/375_R1.bam \
I=${crossdir}/375_R2.bam \
I=${crossdir}/375_R3.bam \
I=${crossdir}/375_R4.bam \
I=${crossdir}/375_R5.bam \
O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000
##################################################################################################################################
cross='376'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
picard ReorderSam I=/data/rr/02/N712/bam/376_G1.bam O=${crossdir}/376_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/09/N709/bam/376_G2.bam O=${crossdir}/376_G2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/09/N710/bam/376_G3.bam O=${crossdir}/376_G3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/09/N711/bam/376_G4.bam O=${crossdir}/376_G4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/09/N712/bam/376_G5.bam O=${crossdir}/376_G5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/03/N712/bam/376_R1.bam O=${crossdir}/376_R1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/08/N701/bam/376_R2.bam O=${crossdir}/376_R2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/08/N702/bam/376_R3.bam O=${crossdir}/376_R3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/08/N703/bam/376_R4.bam O=${crossdir}/376_R4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/06/N712/bam/376_R5.bam O=${crossdir}/376_R5.bam CREATE_INDEX=TRUE REFERENCE=$genome &

picard MergeSamFiles \
I=${crossdir}/376_G1.bam \
I=${crossdir}/376_G2.bam \
I=${crossdir}/376_G3.bam \
I=${crossdir}/376_G4.bam \
I=${crossdir}/376_G5.bam \
I=${crossdir}/376_R1.bam \
I=${crossdir}/376_R2.bam \
I=${crossdir}/376_R3.bam \
I=${crossdir}/376_R4.bam \
I=${crossdir}/376_R5.bam \
O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000
##################################################################################################################################
cross='377'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
picard ReorderSam I=/data/rr/09/N708/bam/377_G1.bam O=${crossdir}/377_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/11/N705/bam/377_G2.bam O=${crossdir}/377_G2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/11/N706/bam/377_G3.bam O=${crossdir}/377_G3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/11/N707/bam/377_G4.bam O=${crossdir}/377_G4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/11/N708/bam/377_G5.bam O=${crossdir}/377_G5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/11/N704/bam/377_R1.bam O=${crossdir}/377_R1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/13/N710/bam/377_R2.bam O=${crossdir}/377_R2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/13/N711/bam/377_R3.bam O=${crossdir}/377_R3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/13/N712/bam/377_R4.bam O=${crossdir}/377_R4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/14/N710/bam/377_R5.bam O=${crossdir}/377_R5.bam CREATE_INDEX=TRUE REFERENCE=$genome &

picard MergeSamFiles \
I=${crossdir}/377_G1.bam \
I=${crossdir}/377_G2.bam \
I=${crossdir}/377_G3.bam \
I=${crossdir}/377_G4.bam \
I=${crossdir}/377_G5.bam \
I=${crossdir}/377_R1.bam \
I=${crossdir}/377_R2.bam \
I=${crossdir}/377_R3.bam \
I=${crossdir}/377_R4.bam \
I=${crossdir}/377_R5.bam \
O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000
##################################################################################################################################
cross='381'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
picard ReorderSam I=/data/rr/11/N709/bam/381_G1.bam O=${crossdir}/381_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/13/N702/bam/381_G2.bam O=${crossdir}/381_G2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/13/N703/bam/381_G3.bam O=${crossdir}/381_G3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/13/N704/bam/381_G4.bam O=${crossdir}/381_G4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/13/N705/bam/381_G5.bam O=${crossdir}/381_G5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/11/N710/bam/381_R1.bam O=${crossdir}/381_R1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/13/N706/bam/381_R2.bam O=${crossdir}/381_R2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/13/N707/bam/381_R3.bam O=${crossdir}/381_R3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/13/N708/bam/381_R4.bam O=${crossdir}/381_R4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/13/N709/bam/381_R5.bam O=${crossdir}/381_R5.bam CREATE_INDEX=TRUE REFERENCE=$genome &

picard MergeSamFiles \
I=${crossdir}/381_G1.bam \
I=${crossdir}/381_G2.bam \
I=${crossdir}/381_G3.bam \
I=${crossdir}/381_G4.bam \
I=${crossdir}/381_G5.bam \
I=${crossdir}/381_R1.bam \
I=${crossdir}/381_R2.bam \
I=${crossdir}/381_R3.bam \
I=${crossdir}/381_R4.bam \
I=${crossdir}/381_R5.bam \
O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000
##################################################################################################################################
cross='393'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
picard ReorderSam I=/data/rr/02/N701/bam/393_G1.bam O=${crossdir}/393_G1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/08/N704/bam/393_G2.bam O=${crossdir}/393_G2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/08/N705/bam/393_G3.bam O=${crossdir}/393_G3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/08/N706/bam/393_G4.bam O=${crossdir}/393_G4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/07/N709/bam/393_G5.bam O=${crossdir}/393_G5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/08/N707/bam/393_R1.bam O=${crossdir}/393_R1.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/10/N705/bam/393_R2.bam O=${crossdir}/393_R2.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/10/N706/bam/393_R3.bam O=${crossdir}/393_R3.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/10/N707/bam/393_R4.bam O=${crossdir}/393_R4.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/10/N708/bam/393_R5.bam O=${crossdir}/393_R5.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/14/N709/bam/393_G5.bam O=${crossdir}/393_G5b.bam CREATE_INDEX=TRUE REFERENCE=$genome &

picard MergeSamFiles \
I=${crossdir}/393_G1.bam  \
I=${crossdir}/393_G2.bam  \
I=${crossdir}/393_G3.bam  \
I=${crossdir}/393_G4.bam  \
I=${crossdir}/393_G5.bam  \
I=${crossdir}/393_R1.bam  \
I=${crossdir}/393_R2.bam  \
I=${crossdir}/393_R3.bam  \
I=${crossdir}/393_R4.bam  \
I=${crossdir}/393_R5.bam  \
I=${crossdir}/393_G5b.bam \
O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000
##################################################################################################################################
cross='B'
crossdir=/data/rr/Segregants/${cross}/
mkdir $crossdir
picard ReorderSam I=/data/rr/00/N701/bam/B_01.bam O=${crossdir}/B_01.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/00/N702/bam/B_02.bam O=${crossdir}/B_02.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/00/N703/bam/B_03.bam O=${crossdir}/B_03.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/00/N704/bam/B_04.bam O=${crossdir}/B_04.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/00/N705/bam/B_05.bam O=${crossdir}/B_05.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/00/N706/bam/B_06.bam O=${crossdir}/B_06.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/00/N707/bam/B_07.bam O=${crossdir}/B_07.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/00/N708/bam/B_08.bam O=${crossdir}/B_08.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/00/N709/bam/B_09.bam O=${crossdir}/B_09.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard ReorderSam I=/data/rr/00/N710/bam/B_10.bam O=${crossdir}/B_10.bam CREATE_INDEX=TRUE REFERENCE=$genome &
picard MergeSamFiles \
I=${crossdir}/B_01.bam \
I=${crossdir}/B_02.bam \
I=${crossdir}/B_03.bam \
I=${crossdir}/B_04.bam \
I=${crossdir}/B_05.bam \
I=${crossdir}/B_06.bam \
I=${crossdir}/B_07.bam \
I=${crossdir}/B_08.bam \
I=${crossdir}/B_09.bam \
I=${crossdir}/B_10.bam \
O=${crossdir}/${cross}.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000


# workflow for BYxRM (Cross A)
# remerged A01_01-A11_96 using samtools merge -r -h etc...
# note, need to resest ulimit -n
# reorder BAM file
# clean BAM file (MAPQ =0 error)
# Ignoring SAM validation error: ERROR: Record 10963391, Read name HWI-ST387_0114:5:47:7351:85174#0, MAPQ should be 0 for unmapped read.
# resort (this shouldn't be necessary ... wtf .. CleanSam bug?)
# this could potentially have removed the offending reads but concerned about picard strictness settings
# samtools view -bF 4 A.bam > A2.bam
# clean and compress 1056 BAM
picard CleanSam I=A.bam O=A3.bam CREATE_INDEX=TRUE MAX_RECORDS_IN_RAM=10000000
picard ReorderSam I=/data/rr/Segregants/A/A3.bam O=/data/rr/Segregants/A/A4.bam CREATE_INDEX=TRUE REFERENCE=$genome
mv /data/rr/Segregants/A/A4.bam /data/rr/Segregants/A/A.bam
samtools index /data/rr/Segregants/A/A.bam

# Now merge giant bams into super giant bam
picard MergeSamFiles \
I=/data/rr/Segregants/375/375.bam \
I=/data/rr/Segregants/A/A.bam \
I=/data/rr/Segregants/376/376.bam \
I=/data/rr/Segregants/B/B.bam \
I=/data/rr/Segregants/377/377.bam \
I=/data/rr/Segregants/393/393.bam \
I=/data/rr/Segregants/381/381.bam \
I=/data/rr/Segregants/3008/3008.bam \
I=/data/rr/Segregants/2999/2999.bam \
I=/data/rr/Segregants/3000/3000.bam \
I=/data/rr/Segregants/3001/3001.bam \
I=/data/rr/Segregants/3049/3049.bam \
I=/data/rr/Segregants/3003/3003.bam \
I=/data/rr/Segregants/3004/3004.bam \
I=/data/rr/Segregants/3043/3043.bam \
I=/data/rr/Segregants/3028/3028.bam \
O=/data/rrv2/genotyping/rr16_segs.bam CREATE_INDEX=TRUE ASSUME_SORTED=TRUE USE_THREADING=TRUE MAX_RECORDS_IN_RAM=1000000
