#!/bin/bash

str1="Please enter the variable"
str2="is located in these parent directories"
str6="Total excluding comments:"

if [[ "$1" != "" ]]; then
	direct="$1"
else
	direct="../ACME/components/clm"
fi

#Call for input
echo $str1
read varIn

#Confirm input variable
printf "\n"
str22="$varIn $str2"
echo $str22

#All top level Directories containing variable
#CHANGE cut to -f1-9 for tmux session, -f1-6 for normal
find $direct -name '*.F90' | xargs grep -l "\s$varIn\>" | cut -d/ -f1-6 |sort -u

#************************************Use Count*********************************************************

#Count all files containing variable
#find . -name '*.F90' | xargs grep -l "\s$varIn\>" | sort -u 
varLoc=$(find $direct -name '*.F90' | xargs grep -l "\b$varIn\b" | sort -u | wc -l)
echo $varLoc "Files"

#Total Count without Comments
printf "\n"
echo $str6
varTot=$(find $direct -name '*.F90' | xargs sed 's/!.*//' | grep -o "\b$varIn\b" | wc -l)
echo $varTot

#Comment Count
varTotc=$(find $direct -name '*.F90' | xargs grep -o "\b$varIn\b" | wc -l)
varCom=$(($varTotc-varTot))

#Return Definition count
varDef=$(find $direct -name '*.F90' | xargs sed 's/!.*//' | grep "\b$varIn\b" | grep -o "\suse\s" |wc -l)

#Return Assignment count
varAs=$(find $direct -name '*.F90' | xargs sed 's/!.*//' | grep "\b$varIn\b" | grep -o "=>" |wc -l)

#Return Association count
varAsc2=$(find $direct -name '*.F90' | xargs sed 's/!.*//' | grep "\b$varIn\b" | grep -o "\s=\s" | wc -l)
varAsc3=$(find $direct -name '*.F90' | xargs sed 's/!.*//' | grep "\b$varIn\b" | grep -o "\s==\s" |wc -l)
varAsc=$(($varAsc1+$varAsc2+$varAsc3))

#Return Other count
varOth=$(($varTot-$varAs-$varDef-$varAsc))


printf "\n"
echo "Variable usage Types"
echo "Definition  | Assignment | Association | Other | Comments"
echo "      " $varDef "         " $varAs "         " $varAsc "       " $varOth "       " $varCom


