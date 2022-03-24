#!/bin/bash
read_ID=""
echo -e "ID\tCHR\tMIN\tMAX\tDISTANCE"
while read line;
do
arr=($line)
if [ "$read_ID" != "" ]; then
   if [ "$read_ID" != "${arr[3]}" ]; then
      echo -e $read_ID"\t"$chr"\t"$pos1"\t"$pos2"\t"$((pos2-pos1))
      read_ID=${arr[3]}
      pos1=${arr[1]}
      pos2=${arr[1]}
      chr=${arr[0]}
   else 
      if [ $pos1 -gt ${arr[1]} ]; then
         pos1=${arr[1]}
      fi
      if [ $pos2 -lt ${arr[1]} ]; then
         pos2=${arr[1]}
      fi
   fi
else
   read_ID=${arr[3]}
   pos1=${arr[1]}
   pos2=${arr[1]}
   chr=${arr[0]}
fi
done < $1
