#!/bin/bash

m8_path=/home/safpt/OperationalChain/LSASAF_MSG-IODC_Products/FRP-PIXEL
m11_path=/home/safpt/OperationalChain/LSASAF_Products/FRP-PIXEL


m11_Qproduct=HDF5_LSASAF_MSG_FRP-PIXEL-QualityProduct_MSG-Disk_
m8_Qproduct=NetCDF_LSASAF_MSG-IODC_FRP-PIXEL-QualityProduct_IODC-Disk_
m11_Lproduct=HDF5_LSASAF_MSG_FRP-PIXEL-ListProduct_MSG-Disk_
m8_Lproduct=/NetCDF_LSASAF_MSG-IODC_FRP-PIXEL-ListProduct_IODC-Disk_


LocalPath=/home/mobaxterm/Documents/Earth_Observation/sandbox/Data

while read File; do
    m8_QP=${m8_Qproduct}${File}
    m8_LP=${m8_Lproduct}${File}
    m11_QP=${m11_Qproduct}${File}
    m11_LP=${m11_Lproduct}${File}


    if [ ! -f "${LocalPath}/${m8_QP}" ] ; then
	scp rhnguyen@193.137.20.100:${m8_path}/${m8_QP}.bz2 ${LocalPath}/
        bunzip2 ${LocalPath}/${m8_QP}.bz2
    fi
    if [ ! -f "${LocalPath}/${m8_LP}" ] ; then
         scp rhnguyen@193.137.20.100:${m8_path}/${m8_LP}.bz2 ${LocalPath}/
	 bunzip2 ${LocalPath}/${m8_LP}.bz2
    fi
    if [ ! -f "${LocalPath}/${m11_QP}" ] ; then
        scp rhnguyen@193.137.20.100:${m11_path}/${m11_QP}.bz2 ${LocalPath}/
	bunzip2 ${LocalPath}/${m11_QP}.bz2
    fi
    if [ ! -f "${LocalPath}/${m11_LP}" ] ; then
        scp rhnguyen@193.137.20.100:${m11_path}/${m11_LP}.bz2 ${LocalPath}/
	bunzip2 ${LocalPath}/${m11_LP}.bz2
    fi

    

done < ./DownloadList.txt
rm ./DownloadList.txt



