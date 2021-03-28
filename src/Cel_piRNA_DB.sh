###C. briggsae piRNA DB; a copy of Cel DB commands for the piRNA app
setwd /home/velazqam/Documents/Projects/piRNA_Jun_2020/Cbri


##Obtain sequence of piRNA sites and named it into a bed

mkdir nbed

#for file in `ls *.bed`; do
#echo ${file%.bed}_seq.bed;
#bedtools getfasta -fi c_briggsae.PRJNA10731.WS274.genomic.fa  -bed ${file} -s -tab | awk -F"\t" '{sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq=$2;k=length(seq); for (i=1;i<=k;i++) {if (substr(seq,i,1)=="T") sumT+=1; else if (substr(seq,i,1)=="A") sumA+=1; else if (substr(seq,i,1)=="G") sumG+=1; else if (substr(seq,i,1)=="C") sumC+=1; else if (substr(seq,i,1)=="N") sumN+=1}; GCcontent=(sumC+sumG)/k*100; print $1"\t"$2"\t"GCcontent}' | perl -pe 's/\(\-\)/-m/g' | perl -pe 's/\(\+\)/-p/g' | awk -F"\t" '{OFS="\t"; split($1,ch,":"); split(ch[2],info,"-"); if(info[3]==m){st="-"}else{st="+"}; print ch[1],info[1],info[2],$2";"$3,".",st}' > nbed/${file%.bed}_seq.bed; done

for file in `ls *piRNA-seqs.bed`; do echo ${file%.bed}_seq.bed; bedtools getfasta -fi c_briggsae.PRJNA10731.WS274.genomic.fa  -bed ${file} -s -tab | awk -F"\t" '{sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq=$2;k=length(seq); for (i=1;i<=k;i++) {if (substr(seq,i,1)=="T") sumT+=1; else if (substr(seq,i,1)=="A") sumA+=1; else if (substr(seq,i,1)=="G") sumG+=1; else if (substr(seq,i,1)=="C") sumC+=1; else if (substr(seq,i,1)=="N") sumN+=1}; GCcontent=(sumC+sumG)/k*100; print $1"\t"$2"\t"GCcontent}' | perl -pe 's/\(\-\)/-m/g' | perl -pe 's/\(\+\)/-p/g' | awk -F"\t" '{OFS="\t"; split($1,ch,":"); split(ch[2],info,"-"); if(info[3]=="m"){st="-"}else{st="+"}; print ch[1],info[1],info[2],$2";"$3,".",st}' > nbed/${file%.bed}_seq.bed; done



##Obtain gene annotations
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS274/species/c_briggsae/PRJNA10731/c_briggsae.PRJNA10731.WS274.annotations.gff3.gz

##Obtain CDSs and asign name of its transcript
zcat c_briggsae.PRJNA10731.WS274.annotations.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="CDS"){print $0}}}' - | awk -F"\t|;|=" '{print $1"\t"$4"\t"$5"\t"$12"\t.\t"$7}' | perl -pe 's/Transcript://g' | awk -F"\t" '{OFS="\t"; split($4,trans,","); for(i=1;i<=length(trans);i++){print $1,$2,$3,trans[i],$5,$6}}' > CDSs.bed

cd nbed
mkdir genomic

bedtools intersect -a ../CDSs.bed -b Cbri_gen_5MM_piRNA-seqs_seq.bed -wao -S | awk -F"\t" '{if($13 == 20){print $0"\t"5MM}}' > genomic/Cbr_5MM.txt
bedtools intersect -a ../CDSs.bed -b Cbri_gen_4MM_piRNA-seqs_seq.bed -wao -S | awk -F"\t" '{if($13 == 20){print $0"\t"4MM}}' > genomic/Cbr_4MM.txt
bedtools intersect -a ../CDSs.bed -b Cbri_gen_3MM_piRNA-seqs_seq.bed -wao -S | awk -F"\t" '{if($13 == 20){print $0"\t"3MM}}' > genomic/Cbr_3MM.txt
bedtools intersect -a ../CDSs.bed -b Cbri_gen_2MM_piRNA-seqs_seq.bed -wao -S | awk -F"\t" '{if($13 == 20){print $0"\t"2MM}}' > genomic/Cbr_2MM.txt
bedtools intersect -a ../CDSs.bed -b Cbri_gen_1MM_piRNA-seqs_seq.bed -wao -S | awk -F"\t" '{if($13 == 20){print $0"\t"1MM}}' > genomic/Cbr_1MM.txt
bedtools intersect -a ../CDSs.bed -b Cbri_gen_0MM_piRNA-seqs_seq.bed -wao -S | awk -F"\t" '{if($13 == 20){print $0"\t"0MM}}' > genomic/Cbr_0MM.txt

mkdir exomic
bedtools intersect -a ../CDSs.bed -b Cbri_exo_5MM_piRNA-seqs_seq.bed -wao -S | awk -F"\t" '{if($13 == 20){print $0"\t"5MM}}' > exomic/Cbr_5MM.txt
bedtools intersect -a ../CDSs.bed -b Cbri_exo_4MM_piRNA-seqs_seq.bed -wao -S | awk -F"\t" '{if($13 == 20){print $0"\t"4MM}}' > exomic/Cbr_4MM.txt
bedtools intersect -a ../CDSs.bed -b Cbri_exo_3MM_piRNA-seqs_seq.bed -wao -S | awk -F"\t" '{if($13 == 20){print $0"\t"3MM}}' > exomic/Cbr_3MM.txt
bedtools intersect -a ../CDSs.bed -b Cbri_exo_2MM_piRNA-seqs_seq.bed -wao -S | awk -F"\t" '{if($13 == 20){print $0"\t"2MM}}' > exomic/Cbr_2MM.txt
bedtools intersect -a ../CDSs.bed -b Cbri_exo_1MM_piRNA-seqs_seq.bed -wao -S | awk -F"\t" '{if($13 == 20){print $0"\t"1MM}}' > exomic/Cbr_1MM.txt
bedtools intersect -a ../CDSs.bed -b Cbri_exo_0MM_piRNA-seqs_seq.bed -wao -S | awk -F"\t" '{if($13 == 20){print $0"\t"0MM}}' > exomic/Cbr_0MM.txt

###Produce genebodies for db
cd ..
zcat c_briggsae.PRJNA10731.WS274.annotations.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="CDS"){print $0}}}' - | awk -F"\t|;|=" '{print $1"\t"$4"\t"$5"\t"$12"\t.\t"$7}' | perl -pe 's/Transcript://g' | awk -F"\t" '{OFS="\t"; split($4,trans,","); for(i=1;i<=length(trans);i++){print $1,$2,$3,trans[i],$5,$6}}' | awk -F"\t" '{if(array[$4] != 0){if(start[$4] > $2){start[$4]=$2};if(end[$4] < $3){end[$4]=$3};array[$4]=$1"\t"start[$4]"\t"end[$4]"\t"$4"\t"$5"\t"$6}else{start[$4]=$2;end[$4]=$3;array[$4]=$0}}END{for(key in array){print array[key]}}' - | sort > GeneBodies.bed
zcat c_briggsae.PRJNA10731.WS274.annotations.gff3.gz | grep -v "#" | grep "WormBase" | awk -F"\t" '{if ($3=="mRNA"){print $9}}'| awk -F";" '{split($1,ts,":"); split($2,gs,":"); split($5,loc,"="); print ts[2]"\t"gs[2]"_"ts[2]"_"loc[2]}' | awk -F"\t" '{OFS="\t"; print $1,$2,$2,$1}' > mRNAID_2.txt
awk -F"\t" '{OFS="\t"; if(array[$4]==0){array[$4]=$2}else{print $1,$2,$3,array[$4],$5,$6}}' mRNAID_2.txt GeneBodies.bed > GeneBodies_WBID_2.bed
awk -F"\t" '{OFS="\t"; if(array[$4]==0){array[$4]=$2}else{print $1,$2,$3,array[$4],$5,$6,$4}}' mRNAID_2.txt GeneBodies.bed > GeneBodies_WBID.bed


#Do analysis for genomic
cd genomic
awk -F"\t" '{if(array[$4"-"$8"-"$10] == ""){array[$4"-"$8"-"$10]=$0; track[$4"-"$8"-"$10]=$14}else{track[$4"-"$8"-"$10]=$14}} END{for(key in array){print array[key]"\t"track[key]}}' Cbr_0MM.txt Cbr_1MM.txt Cbr_2MM.txt Cbr_3MM.txt Cbr_4MM.txt Cbr_5MM.txt > combine.txt

sort -k1,1 -k2,2n combine.txt | awk -F"\t" '{if(NF==7){data[$7]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}else{print data[$4]"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}}' ../../GeneBodies_WBID.bed - > temp.txt

mkdir DB
cd DB
awk -F"\t" '{OFS="\t"; array[$4]=array[$4]""$8"\t"$10"\t"$15"\n"}END{for(keys in array){print array[keys] > keys".txt"}}' ../temp.txt
cd ..
#Final gene table
awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$4"\t"$6]++}END{for(keys in array){print keys}}' temp.txt  | sort -k1,1 -k2,2n | awk -F"_" '{OFS="\t"; print $1,$2,$3}' > Genes.txt

cd ..
##Now for exomic
cd exomic
awk -F"\t" '{if(array[$4"-"$8"-"$10] == ""){array[$4"-"$8"-"$10]=$0; track[$4"-"$8"-"$10]=$14}else{track[$4"-"$8"-"$10]=$14}} END{for(key in array){print array[key]"\t"track[key]}}' Cbr_0MM.txt Cbr_1MM.txt Cbr_2MM.txt Cbr_3MM.txt Cbr_4MM.txt Cbr_5MM.txt > combine.txt

sort -k1,1 -k2,2n combine.txt | awk -F"\t" '{if(NF==7){data[$7]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}else{print data[$4]"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15}}' ../../GeneBodies_WBID.bed - > temp.txt

mkdir DB
cd DB
awk -F"\t" '{OFS="\t"; array[$4]=array[$4]""$8"\t"$10"\t"$15"\n"}END{for(keys in array){print array[keys] > keys".txt"}}' ../temp.txt
cd ..
#Final gene table
awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$4"\t"$6]++}END{for(keys in array){print keys}}' temp.txt  | sort -k1,1 -k2,2n | awk -F"_" '{OFS="\t"; print $1,$2,$3}' > Genes.txt

cd ..

#Now combine them
mkdir combine
cd exomic
mkdir strange
cd DB
for file in `ls *.txt`; do awk -F"\t" '{print "comb\t"$0}' ${file} > ../strange/${file}; done
cd ../strange

##For each file in here get its partner
for file in `ls *.txt`; do awk -F"\t" '{OFS="\t"; if(NF == 3){ stata[$2]=$1; gen[$2]=$3; exo[$2]=-1}; if(NF == 4){  if(stata[$3] == ""){stata[$3]=$2; gen[$3]=-1; exo[$3]=$4}else{exo[$3]=$4}  }} END{for(keys in stata){print stata[keys]"\t"keys"\t"gen[keys]"\t"exo[keys]}}' ../../genomic/DB/${file} ${file} | sort -k1,1n > ../../combine/${file}; done

cd ../../

cp exomic/Genes.txt .

cd ..

##Now make Introns and alias
zcat c_briggsae.PRJNA10731.WS274.annotations.gff3.gz | awk -F"\t" '{OFS="\t"; if($2=="WormBase"){if($3=="intron"){print $1,$4,$5,$9,$8,$7}}}' | awk -F"\t" '{OFS="\t"; split($4,info,";"); split(info[1],dat,":"); print $1,$2,$3,dat[2],$5,$6}' | sort -k1,1 -k2,2n | awk -F"\t" '{if(st[$4] == ""){chr[$4]=$1; st[$4]=$2; ed[$4]= $3; str[$4]=$6}else{st[$4]=st[$4]","$2; ed[$4]=ed[$4]","$3;}} END{for (key in st){print chr[key]"\t"st[key]"\t"ed[key]"\t"key"\t.\t"str[key]}}' | sort -k1,1 -k4,4 > Introns.bed
awk -F"\t" '{OFS="\t"; print $4,$2,$3}' Introns.bed > Introns_min.tsv

zcat c_briggsae.PRJNA10731.WS274.annotations.gff3.gz | awk -F"\t" '{OFS="\t"; if($2=="WormBase"){if($3=="gene"){print $9}}}' | awk -F";" '{print $1"\t"$NF}' | awk -F"\t" '{split($1,info,":"); split($2,alas,"="); split(alas[2],array,","); for(extr in array){print info[2]"\t"array[extr]}}' | awk -F"\t" '{count[$2]++; data[$2]=$1} END{for(key in count){if(count[key]==1){print key"\t"data[key]}}}' | sort -k2,2 > Alias.txt

##Strange naming for chromosomes
cd nbed/genomic/
awk -F"\t" '{array[$1"\t"$2"\t"$3"\t"$4"\t"$6]++}END{for(keys in array){print keys}}' temp.txt  | sort -k1,1 -k2,2n | awk -F"\t" '{split($4,info,"_"); OFS="\t"; print $1,$2,$3,info[1],info[2],info[3],$5}' > Genes.txt

cp Genes.txt ../

