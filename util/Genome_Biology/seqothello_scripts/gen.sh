#!/bin/bash

TOOLCHAIN_PATH_DEFAULT=`pwd`/build/bin
read -e -p "Please enter the folder contains SeqOthello toolchain [default: $TOOLCHAIN_PATH_DEFAULT]: " INPUT
TOOLCHAIN_PATH="${INPUT:-$TOOLCHAIN_PATH_DEFAULT}"
while [ ! -f ${TOOLCHAIN_PATH}/Build ]; do 
    echo ${TOOLCHAIN_PATH}'/Build not found.'
    read -e -p "Please enter the folder contains SeqOthello toolchain [default: $TOOLCHAIN_PATH_DEFAULT]: " INPUT
    TOOLCHAIN_PATH="${INPUT:-$TOOLCHAIN_PATH_DEFAULT}"
done

KMER_PATH_DEFAULT=`pwd`/example/kmer
read -e -p "Please enter the folder that contains all Jellyfish generated Kmer files [default: $KMER_PATH_DEFAULT]: " INPUT
KMER_PATH="${INPUT:-$KMER_PATH_DEFAULT}"
while [ ! -d ${KMER_PATH} ]; do 
    echo ${KMER_PATH} 'not found.'
    read -e -p "Please enter the folder that contains all Jellyfish generated Kmer files [default: $KMER_PATH_DEFAULT]: " INPUT
    KMER_PATH="${INPUT:-$KMER_PATH_DEFAULT}"
done

KMER_FLIST_DEFAULT=${KMER_PATH}/flist
read -e -p "Please enter a file that contains all filenames of the Kmer files [default: ${KMER_FLIST_DEFAULT}]: " INPUT
KMER_FLIST="${INPUT:-$KMER_FLIST_DEFAULT}"

while [ ! -f ${KMER_FLIST} ]; do 
    echo ${KMER_FLIST} 'does not exist.'
    read -e -p "Please enter a file that contains all filenames of the Kmer files [default: ${KMER_FLIST_DEFAULT}]: " INPUT
    KMER_FLIST="${INPUT:-$KMER_FLIST_DEFAULT}"
done

read -e -p "Where can we keep some temporary files? [default: `pwd`/example]: " INPUT
TEMP_FOLDER="${INPUT:-`pwd`/example}"

if [ ! -d ${TEMP_FOLDER}/bin/ ]; then 
    echo 'folder' ${TEMP_FOLDER}'/bin/ does not exist, creating.'
    mkdir -p ${TEMP_FOLDER}/bin
fi


read -e -p "Please enter the value of k [default : 20]: " INPUT
k=${INPUT:-20}

NUM_CONVERTED=0
CONVERT_TO_BINARY=example/ConvertToBinary.sh
BINARY_LIST=example/BinaryList
BINARY_LIST_PREFIX=BinaryList.Part.

if [ -f ${BINARY_LIST} ]; then
    echo ${BINARY_LIST} 'exists. rebuilding this file.'
    rm -rf ${BINARY_LIST}*
fi
if [ -f ${CONVERT_TO_BINARY} ]; then
    echo ${CONVERT_TO_BINARY} 'exists. rebuilding this file.'
    rm -rf ${CONVERT_TO_BINARY}
fi

for kmerfile in `cat ${KMER_FLIST}`; do
    binfile="${kmerfile%.*}".Bin
    echo ${TOOLCHAIN_PATH}/PreProcess --in=${KMER_PATH}/${kmerfile} --out=${TEMP_FOLDER}/bin/$binfile --k=${k} >> ${CONVERT_TO_BINARY}
    echo $binfile >>  ${BINARY_LIST}
    ((NUM_CONVERTED++))
done

echo 'Step 1: Convert the kmer files to binaryKmer files. '
echo '        Prepared the script of converting '${NUM_CONVERTED}' kmer files' as ${CONVERT_TO_BINARY}
echo '        These binaryKmer files are listed in '  ${BINARY_LIST}

echo 'Step 2: Group the binaryKmer files.'
echo ${TEMP_FOLDER}
FILE_PER_GROUP=30
while [ $((FILE_PER_GROUP * FILE_PER_GROUP)) -le $NUM_CONVERTED  ] ; do
    ((FILE_PER_GROUP++))
done


split -l ${FILE_PER_GROUP} -d ${BINARY_LIST} ${TEMP_FOLDER}/${BINARY_LIST_PREFIX}
echo '        Each group contains at most '$FILE_PER_GROUP 'files.'
echo '        These group description files are' ${TEMP_FOLDER}/${BINARY_LIST_PREFIX}'*' 

if [ ! -d ${TEMP_FOLDER}/grp ]; then 
    echo 'folder' ${TEMP_FOLDER}'/grp/ does not exist, creating.'
    mkdir  -p ${TEMP_FOLDER}/grp
fi

MAKE_GROUP=example/MakeGroup.sh
if [ -f ${MAKE_GROUP} ]; then
    echo ${MAKE_GROUP} 'exists. rebuilding this file.'
    rm -rf ${MAKE_GROUP}
fi

GRPLIST=`pwd`/example/GrpList
if [ -d ${GRPLIST} ]; then
    echo ${GRPLIST} 'exists.' rebuilding this file.
    rm -rf ${GRPLIST}
fi
GRP_CONVERTED=0
echo ${TEMP_FOLDER}
for flist in `ls -m1 ${TEMP_FOLDER}/${BINARY_LIST_PREFIX}*`; do 
    echo ${TOOLCHAIN_PATH}/Group --flist=$flist --folder=${TEMP_FOLDER}/bin/ --output=${TEMP_FOLDER}/grp/Grp"${flist##*.}" >> ${MAKE_GROUP}
    echo Grp"${flist##*.}" >> ${GRPLIST}
    ((GRP_CONVERTED++))
done

echo '        Prepared the script of making '${GRP_CONVERTED}' groups' as ${MAKE_GROUP}
echo '        These group files are listed in '  ${GRPLIST}
echo '        These group are located in in '  ${TEMP_FOLDER}/grp

echo 'Step 3: Build the SeqOthello structure'

BUILD_SCRIPT=example/BuildSeqOthello.sh

EXPORT_FOLDER_DEFAULT=`pwd`/example/out
read -e -p "Please enter the folder to put the SeqOthello files. [default: $EXPORT_FOLDER_DEFAULT]: " INPUT
EXPORT_FOLDER="${INPUT:-$EXPORT_FOLDER_DEFAULT}"

if [ ! -d ${EXPORT_FOLDER} ]; then 
    echo 'folder' ${EXPORT_FOLDER} 'does not exist, creating.'
    mkdir -p ${EXPORT_FOLDER}
fi

while [ ! -d ${KMER_PATH} ]; do 
    echo ${KMER_PATH} 'not found.'
    read -e -p "Please enter the folder that contains all Jellyfish generated Kmer files [default: $KMER_PATH_DEFAULT]: " INPUT
done
echo ${TOOLCHAIN_PATH}/Build --flist=${GRPLIST} --folder=${TEMP_FOLDER}/grp/ --out-folder=$EXPORT_FOLDER/ '>' $EXPORT_FOLDER/Build.log > ${BUILD_SCRIPT}

chmod +x ${CONVERT_TO_BINARY} ${BUILD_SCRIPT} ${MAKE_GROUP}
echo
echo 'Summary: Generated the following three scripts to build the SeqOthello structure'
echo '         1.  ' ${CONVERT_TO_BINARY} ': ' `wc ${CONVERT_TO_BINARY} -l | cut -d' ' -f1` 'lines.'
echo '         2.  ' ${MAKE_GROUP} ': ' `wc ${MAKE_GROUP} -l | cut -d' ' -f1` 'lines.'
echo '         3.  ' ${BUILD_SCRIPT} ': ' `wc ${BUILD_SCRIPT} -l | cut -d' ' -f1` 'lines.'
echo
echo ' Please run these three scripts one by one, for example, using the following commands. '
echo '         # ./'${CONVERT_TO_BINARY}
echo '         # ./'${MAKE_GROUP}
echo '         # ./'${BUILD_SCRIPT}
echo
echo ' Note: the commands within each script can be executed in parallel.'
echo '       For example, using GNU Parallel, you can run the scripts as:'
echo '         # cat '${CONVERT_TO_BINARY}' | parallel -j 16'
echo '         # cat '${MAKE_GROUP}' | parallel -j 16'
echo '         # ./'${BUILD_SCRIPT}


