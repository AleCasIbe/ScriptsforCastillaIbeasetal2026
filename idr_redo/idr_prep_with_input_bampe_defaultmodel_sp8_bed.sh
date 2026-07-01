#!/bin/bash


### Se crea una carpeta temporal, que se borrara al final del programa, para almacenar el input y output: SCRATCH
id=$3
SCRATCH=./scratch$id

mkdir -p $SCRATCH || exit $?

### Guardo en run una carpeta por punto, con los .bed de las replicas. Se recoge el nombre de los archivos y los archivos mismos.
### Ademas se recoge el working directory para grabar los datos.

wd=$PWD
folder=$1
name_folder=${folder##*/}
file1=$folder/rochp_sp8_rep1_sorted_filt_chr.bam
file2=$folder/rochp_sp8_rep2_sorted_filt_chr.bam
input1=$folder/roinp_sp8_rep1_sorted_filt_chr.bam
input2=$folder/roinp_sp8_rep2_sorted_filt_chr.bam

### Cuando se llama al programa se le da el tamanio del genoma como argumento. 

gsize=$2

### Se copian los datos requeridos en el archivo temporal SCRATCH

cp $file1 $SCRATCH/1.bam
cp $file2 $SCRATCH/2.bam
cp $input1 $SCRATCH/i1.bam
cp $input2 $SCRATCH/i2.bam

### Se crea la carpeta results en SCRATCH. 

cd $SCRATCH
cat i1.bam i2.bam > input.bam

mkdir results
cd results

### Peak calling en las replicas y en el pool de replicas: IDR conservative peaks --> high confidence peaks, that represent
### events across true biological replciates and that account for true biological and technical noise.

#activate the virtual env for python 2 which has MACS2
#conda activate IDR_Juan


echo 'Starting proper replicate peak calling'
mkdir CONS_PEAKS
macs2 callpeak -f BAM -t ../1.bam -c ../i1.bam -p 0.02 --keep-dup auto -g ${gsize} -n ${name_folder}_rep1 --outdir ./CONS_PEAKS &
macs2 callpeak -f BAM -t ../2.bam -c ../i2.bam -p 0.02 --keep-dup auto -g ${gsize} -n ${name_folder}_rep2 --outdir ./CONS_PEAKS &
macs2 callpeak -f BAM -t ../1.bam ../2.bam -c ../input.bam -p 0.02 --keep-dup auto -g ${gsize} -n ${name_folder}_pooled_rep --outdir ./CONS_PEAKS &
wait %1 %2 %3
echo 'Done with proper replicate peak calling'

### Se filtran los 5E5 mejores picos (p-value en 8 columna) --> nr indica numerico y orden reverso

echo 'Extracting top 500000 peaks'
sort -k8,8nr ./CONS_PEAKS/${name_folder}_rep1_peaks.narrowPeak | head -n 500000 > ./CONS_PEAKS/${name_folder}_rep1_toppeaks.bed &
sort -k8,8nr ./CONS_PEAKS/${name_folder}_rep2_peaks.narrowPeak | head -n 500000 > ./CONS_PEAKS/${name_folder}_rep2_toppeaks.bed &
sort -k8,8nr ./CONS_PEAKS/${name_folder}_pooled_rep_peaks.narrowPeak | head -n 500000 > ./CONS_PEAKS/${name_folder}_pool_toppeaks.bed &
wait %1 %2 %3
echo 'Finished extracting top 500000 peaks'

### Un segundo set de picos son los denominados picos optimos --> high-confidence peaks, that represent
### reproducible events and that accounting for read sampling noise. The optimal set is more sensitive, 
### especially when one of the replicates has lower data quality than the other.
### En este caso se han de comparar pseudoreplicas. Para cada replica biologica, se obtienen 2 pseudoreplicas
### para lo cual se barajan los datos de dicha replica y se divide en dos archivos con el mismo numero de lineas.
### Los datos que se barajan son los originales, los nucfree.bed, que contiene coordenadas de reads, pero estan de una en una,
### no agrupadas con su dato de frecuencia (como hace MACS2). Por eso es posible barajar y dividir por la mitad: si una zona
### tiene un pico de verdad, al dividir las reads en dos, si es suficientemente fuerte, los dos ficheros tendran suficientes 
### reads en esa zona como para dar lugar a un pico.

## Creacion de pseudoreplicas

echo 'Creating pseudoreplicates'
mkdir OPT_PEAKS

nlines1=$( cat ../1.bam | wc -l ) 
halfnlines1=$(( (nlines1 + 1) / 2 )) 
cat ../1.bam | shuf | split -d -l ${halfnlines1} - ./OPT_PEAKS/${name_folder}_rep1_ps  ## -d indica que se pondran sufijos numericos

nlines2=$( cat ../2.bam | wc -l ) 
halfnlines2=$(( (nlines2 + 1) / 2 )) 
cat ../2.bam | shuf | split -d -l ${halfnlines2} - ./OPT_PEAKS/${name_folder}_rep2_ps

cat ../1.bam ../2.bam | shuf > ./pooled.bam
nlinesp=$( cat ./pooled.bam | wc -l ) 
halfnlinesp=$(( (nlinesp + 1) / 2 )) 
cat ./pooled.bam | shuf | split -d -l ${halfnlinesp} - ./OPT_PEAKS/${name_folder}_pool_ps 

echo 'Done creating pseudoreplicates'

## Peak calling de las pseudoreplicas

echo 'Starting pseudoreplicate peak calling'
macs2 callpeak -f BAM -t ./OPT_PEAKS/${name_folder}_rep1_ps00 -c ../i1.bam --keep-dup auto -p 0.02 -g ${gsize} -n ${name_folder}_rep1_ps00 --outdir ./OPT_PEAKS &
macs2 callpeak -f BAM -t ./OPT_PEAKS/${name_folder}_rep1_ps01 -c ../i1.bam --keep-dup auto -p 0.02 -g ${gsize} -n ${name_folder}_rep1_ps01 --outdir ./OPT_PEAKS &

macs2 callpeak -f BAM -t ./OPT_PEAKS/${name_folder}_rep2_ps00 -c ../i2.bam --keep-dup auto -p 0.02 -g ${gsize} -n ${name_folder}_rep2_ps00 --outdir ./OPT_PEAKS &
macs2 callpeak -f BAM -t ./OPT_PEAKS/${name_folder}_rep2_ps01 -c ../i2.bam --keep-dup auto -p 0.02 -g ${gsize} -n ${name_folder}_rep2_ps01 --outdir ./OPT_PEAKS &

macs2 callpeak -f BAM -t ./OPT_PEAKS/${name_folder}_pool_ps00 -c ../input.bam --keep-dup auto -p 0.02 -g ${gsize} -n ${name_folder}_pool_ps00 --outdir ./OPT_PEAKS &
macs2 callpeak -f BAM -t ./OPT_PEAKS/${name_folder}_pool_ps01 -c ../input.bam --keep-dup auto -p 0.02 -g ${gsize} -n ${name_folder}_pool_ps01 --outdir ./OPT_PEAKS &

wait %1 %2 %3 %4 %5 %6
echo 'Done with pseudoreplicate peak calling'

echo 'Extracting top 500000 peaks'
sort -k8,8nr ./OPT_PEAKS/${name_folder}_rep1_ps00_peaks.narrowPeak | head -n 500000 > ./OPT_PEAKS/${name_folder}_rep1_ps00_toppeaks.bed &
sort -k8,8nr ./OPT_PEAKS/${name_folder}_rep1_ps01_peaks.narrowPeak | head -n 500000 > ./OPT_PEAKS/${name_folder}_rep1_ps01_toppeaks.bed &

sort -k8,8nr ./OPT_PEAKS/${name_folder}_rep2_ps00_peaks.narrowPeak | head -n 500000 > ./OPT_PEAKS/${name_folder}_rep2_ps00_toppeaks.bed &
sort -k8,8nr ./OPT_PEAKS/${name_folder}_rep2_ps01_peaks.narrowPeak | head -n 500000 > ./OPT_PEAKS/${name_folder}_rep2_ps01_toppeaks.bed &

sort -k8,8nr ./OPT_PEAKS/${name_folder}_pool_ps00_peaks.narrowPeak | head -n 500000 > ./OPT_PEAKS/${name_folder}_pool_ps00_toppeaks.bed &
sort -k8,8nr ./OPT_PEAKS/${name_folder}_pool_ps01_peaks.narrowPeak | head -n 500000 > ./OPT_PEAKS/${name_folder}_pool_ps01_toppeaks.bed &

wait %1 %2 %3 %4 %5 %6
echo 'Finished extracting top 500000 peaks'
