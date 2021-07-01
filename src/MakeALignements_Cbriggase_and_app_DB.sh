##For C. briggsae

#In Ibex

cd /ibex/scratch/velazqam/piRNAi_master/c_briggsae/parts
for num in {00..99}; do cat ${num}/Histogram.txt >> ../20mer_cds_unrevts.Histogram.Cbriggsae.txt; done
cd ..
scp 20mer_cds_unrevts.Histogram.Cbriggsae.txt velazqam@dm:/home/velazqam/Desktop/Files

mkdir Process_histogram


#In my computer
cd /home/velazqam/Documents/Projects/piRNAi_March_2020/c_briggsae/Process_histogram

scp velazqam@dm:/home/velazqam/Desktop/Files/20mer_cds_unrevts.Histogram.Cbriggsae.txt .

awk -F"\t" '{print $2}' ../20mer_per_gene.tab | perl -pe 'tr/ATGCatgc/TACGtacg/' | rev | paste ../20mer_per_gene.tab - > 20merRevCom_per_gene.tab

awk '{if(NF == 3){array[$3]=$1"\t"$2}else{print array[$1]"\t"$0}}' 20merRevCom_per_gene.tab 20mer_cds_unrevts.Histogram.Cbriggsae.txt | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > Unique20Mers_CDS_WS270.tab


awk '{if(/>/){split($0,toto,">"); split(toto[2],info, " "); split(info[2],gene,"="); printf "\n"gene[2]"\t"info[1]"\t"}else{printf $0}}END{print ""}' ../c_briggsae.PRJNA10731.WS270.CDS_transcripts.fa | tail -n +2 > Gene_transcript_cds.tab


awk -F"\t" '{OFS="\t"; for(i=1; i<=(length($3)-19); i++){seq=substr($3,i,20); print $1,$2,seq,i}}' Gene_transcript_cds.tab > Gene_transcript_cds_pos.tab

#~8GB ram


##What format I want for the output fdile

awk -F"\t" '{if(NF==4){hash[$1"-"$3]=$0}else{print hash[$1"-"$2]"\t"$0}}' Gene_transcript_cds_pos.tab Unique20Mers_CDS_WS270.tab | awk -F"\t" '{OFS="\t"; sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq=$7;k=length(seq); for (i=1;i<=k;i++) {if (substr(seq,i,1)=="T") sumT+=1; else if (substr(seq,i,1)=="A") sumA+=1; else if (substr(seq,i,1)=="G") sumG+=1; else if (substr(seq,i,1)=="C") sumC+=1; else if (substr(seq,i,1)=="N") sumN+=1}; GCcontent=(sumC+sumG)/k*100; print $1,$2,$4,$7";"GCcontent,$8,$9,$10,$11,$12,$13}' > Master_piRNA_unsorted.tsv

sort -k1,1 -k2,2 -k3,3n Master_piRNA_unsorted.tsv | awk -F"\t" '{if(/WBGene/){print $0}}' - > Master_piRNAi_Cbr_WS270.tsv



##For app
mkdir for_app
cd for_app

awk -F"\t" '{if($2 =="WormBase"){if($3 =="gene"){print $0}}}' ../../c_briggsae.PRJNA10731.WS270.annotations.gff3 | awk -F"\t" '{split($9,alias,"Alias="); split($9,info,";"); split(info[1],gene,":"); print gene[2]"\t"alias[2]}' > WBID_Alias.tsv

awk -F"\t" '{OFS="\t"; mm=1; for(i=6; i<= NF; i++){if( $i > 0){break}else{mm++}}; print $1,$2,$3,$4,mm}' ../Master_piRNAi_Cbr_WS270.tsv > Master_piRNAi_Cbr_WS270_forDB.tsv 



###

##Old FOrmat
#5107845	GAGTTCCAGGTGTTGTGCTC;55	4	5
#PosGenomic	Seq;gc	MM_gen	MM_exo

##
#FIle name
#WBGene00000001_Y110A7A.10.1_aap-1.txt
#WBID_TranscriptID_locus.txt

##New Format
#3	GAGTTCCAGGTGTTGTGCTC;55	2
#POsTranscript	Seq;gc	MM_exo

#WBGene00000001_Y110A7A.10.1.txt
#WBID_TranscriptID.txt

mkdir DB

awk -F"\t" '{OFS="\t"; print $3,$4,$5 >> "DB/"$1"_"$2".txt"}' Master_piRNAi_Cbr_WS270_forDB.tsv

ls DB/*.txt | wc -l

cp DB/*.txt /home/velazqam/Documents/Projects/Git_repositories/piRNAi_dev/DataBase/

awk -F"\t" '{hash[$2]=$1} END{for(key in hash){print hash[key]"\t"key}}' Master_piRNAi_Cbr_WS270_forDB.tsv | awk -F"\t" '{if(hash[$1]==""){hash[$1]=$2}else{print $0"\t"hash[$1]}}' WBID_Alias.tsv - | sort -k1,1 > Gene_DB.tsv

cat Gene_DB.tsv >> /home/velazqam/Documents/Projects/Git_repositories/piRNAi_dev/Gene_DB.tsv


##IN APP to finalize

awk -F"\t" '{split($3,locus,","); for(i=1; i <= length(locus); i++){print $1"\t"$2"\t"locus[i]}}' Gene_DB.tsv > Main_DB.tsv

cd ~/Documents/Projects/piRNAi_March_2020/c_briggsae/Process_histogram

awk -F"\t" '{print $1"\t"$2"\t"length($3)}' Gene_transcript_cds.tab > CDS_sizes.tab

cat CDS_sizes.tab >> /home/velazqam/Documents/Projects/Git_repositories/piRNAi_dev/WorkingSpace/CDS_sizes.tab 

##Now entry of genes for tracks
###But for that I require to make an mRNA DB as I did for C. elegans

awk -F"\t" '{OFS="\t"; mm=0; for(i=5; i<= 5; i++){mm=mm+$i}; if(mm < 2){print $1,$2,$3,$4,($6+1)}}' Master_piRNAi_Cel_WS270.tsv > 1MM.seqs.txt

awk -F"\t" '{OFS="\t"; mm=0; for(i=5; i<= 6; i++){mm=mm+$i}; if(mm < 2){print $1,$2,$3,$4,($7+1)}}' Master_piRNAi_Cel_WS270.tsv > 2MM.seqs.txt

awk -F"\t" '{OFS="\t"; mm=0; for(i=5; i<= 7; i++){mm=mm+$i}; if(mm < 2){print $1,$2,$3,$4,($8+1)}}' Master_piRNAi_Cel_WS270.tsv > 3MM.seqs.txt

awk -F"\t" '{OFS="\t"; mm=0; for(i=5; i<= 8; i++){mm=mm+$i}; if(mm < 2){print $1,$2,$3,$4,($9+1)}}' Master_piRNAi_Cel_WS270.tsv > 4MM.seqs.txt

awk -F"\t" '{OFS="\t"; mm=0; for(i=5; i<= 9; i++){mm=mm+$i}; if(mm < 2){print $1,$2,$3,$4,($10+1)}}' Master_piRNAi_Cel_WS270.tsv > 5MM.seqs.txt

awk -F"\t" '{OFS="\t"; mm=0; for(i=5; i<= 10; i++){mm=mm+$i}; if(mm < 2){print $1,$2,$3,$4,($11+1)}}' Master_piRNAi_Cel_WS270.tsv > 6MM.seqs.txt


awk -F"\t" '{print $1}' 0MM.seqs.txt | sort | uniq | wc -l

awk -F"\t" '{print $1}' 1MM.seqs.txt | sort | uniq | wc -l

awk -F"\t" '{print $1}' 2MM.seqs.txt | sort | uniq | wc -l

awk -F"\t" '{print $1}' 3MM.seqs.txt | sort | uniq | wc -l

awk -F"\t" '{print $1}' 4MM.seqs.txt | sort | uniq | wc -l




awk -F"\t" '{split($4,info,";"); print info[1]}' 4MM.seqs.txt | perl -pe 'tr/ATGCatgc/TACGtacg/' | rev | paste - 4MM.seqs.txt | awk -F"\t" '{print $1"\t"$6}' > 4MM_seqs.txt

awk -F"\t" '{split($4,info,";"); print info[1]}' 4MM.seqs.txt | perl -pe 'tr/ATGCatgc/TACGtacg/' | rev | paste - 4MM.seqs.txt | awk -F"\t" '{print $1"\t"$6}' > 4MM_seqs.txt

awk -F"\t" '{split($4,info,";"); print info[1]}' 3MM.seqs.txt | perl -pe 'tr/ATGCatgc/TACGtacg/' | rev | paste - 3MM.seqs.txt | awk -F"\t" '{print $1"\t"$6}' > 3MM_seqs.txt

awk -F"\t" '{split($4,info,";"); print info[1]}' 2MM.seqs.txt | perl -pe 'tr/ATGCatgc/TACGtacg/' | rev | paste - 2MM.seqs.txt | awk -F"\t" '{print $1"\t"$6}' > 2MM_seqs.txt

awk -F"\t" '{split($4,info,";"); print info[1]}' 1MM.seqs.txt | perl -pe 'tr/ATGCatgc/TACGtacg/' | rev | paste - 1MM.seqs.txt | awk -F"\t" '{print $1"\t"$6}' > 1MM_seqs.txt

awk -F"\t" '{split($4,info,";"); print info[1]}' 0MM.seqs.txt | perl -pe 'tr/ATGCatgc/TACGtacg/' | rev | paste - 0MM.seqs.txt | awk -F"\t" '{print $1"\t"$6}' > 0MM_seqs.txt

