#!/bin/bash
i=0;

for file in ~/Desktop/APOE_Positions/APOE44_14mon_061719_positions_MWM/*.csv;
do
    newfilename=$(echo $file | sed -e "s/\ /_/g")
    echo $file $newfilename
    mv "$file" "$newfilename"
done
#mv ${subjpath}/$file ${subjpath}/$newFile
i=1