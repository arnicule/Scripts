#!/bin/bash
i=0;

filepath="/Users/bass/Desktop/APOE_Positions/APOE44_14mon_061719_positions_MWM/"

while IFS=, read -r Test Animal Treatment Code Stage Trial test1;
do
    if [ $i -gt 0 ]
    then
	      # Use if data file is missing certain tests
	      if [ $i -eq 109 ] || [ $i -eq 110 ] || [ $i -eq 111 ]
	      then
	          i=112
	      fi
	      oldFile="APOE44_14mon_061719_positions_MWM_-_Test_"$i.csv
        newFile=$Animal"_"$Stage"_T"$Trial"_positions".csv
        echo $oldFile
        echo $newFile
        mv "${filepath}$oldFile" "${filepath}$newFile"
        i=$((i+1))
    else
        i=$((i+1))
    fi
done < ~/Desktop/APOE_Data/MWM/APOE44_14mon_061019MWM.csv;