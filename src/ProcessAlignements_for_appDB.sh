##For app
cd Process_Histogram/for_app

awk -F"\t" '{if($2 =="WormBase"){if($3 =="gene"){print $0}}}' c_elegans.PRJNA13758.WS270.annotations.gff3 | awk -F"\t" '{split($9,alias,"Alias="); split($9,info,";"); split(info[1],gene,":"); print gene[2]"\t"alias[2]}' > WBID_Alias.tsv

awk -F"\t" '{OFS="\t"; mm=1; for(i=6; i<= NF; i++){if( $i > 0){break}else{mm++}}; print $1,$2,$3,$4,mm}' Master_piRNAi_Cel_WS270.tsv > Master_piRNAi_Cel_WS270_forDB.tsv 


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

awk -F"\t" '{OFS="\t"; print $3,$4,$5 >> "DB/"$1"_"$2".txt"}' Master_piRNAi_Cel_WS270_forDB.tsv

ls DB/*.txt | wc -l

cp DB/*.txt piRNAi_dev/DataBase/

awk -F"\t" '{hash[$2]=$1} END{for(key in hash){print hash[key]"\t"key}}' Master_piRNAi_Cel_WS270_forDB.tsv | awk -F"\t" '{if(hash[$1]==""){hash[$1]=$2}else{print $0"\t"hash[$1]}}' WBID_Alias.tsv - | sort -k1,1 > Gene_DB.tsv

cp Gene_DB.tsv piRNAi_dev

cd /piRNAi_March_2020/Process_Histogram

awk -F"\t" '{print $1"\t"$2"\t"length($3)}' Gene_transcript_cds.tab > CDS_sizes.tab

cp CDS_sizes.tab piRNAi_dev/WorkingSpace/

##Now entry of genes

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



