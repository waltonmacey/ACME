#!/bin/bash


str1="Please enter the variable"
str2="is located in these locations"
str3="Count as an Assignment"
str4="Count as a Definition"
str5="Count as an Association"
str6="Total occurences:"
str6c="Total occurences with comments:"
str7="Count as Other Type"

#Call for input
echo $str1
read varIn

#Confirm input variable
printf "\n"
str22="$varIn $str2"
echo $str22

#All top level Directories containing variable
find ../ACME/components/clm -name '*.F90' | xargs grep -l "\s$varIn\>" | cut -d/ -f1-6 |sort -u

#Alternative displaying all files
#varDir=$(find . -name '*.F90' | xargs grep -l "\s$varIn\>" | sort -u)
#cut -d/ -f1,2,3,4,5 <<< "${varDir}" | sort -u

#************************************Use Count*********************************************************
#Count all files containing variable
#find . -name '*.F90' | xargs grep -l "\s$varIn\>" | sort -u 
varLoc=$(find ../ACME/components/clm -name '*.F90' | xargs grep -l "\s$varIn\>" | sort -u | wc -l)
echo $varLoc "Files"

#Return total count excluding comments
printf "\n"
echo $str6
varTot=$(find ../ACME/components/clm -name '*.F90' | xargs sed 's/!.*//' | grep -o "\s$varIn\>" | wc -l)
echo $varTot

#Return total count including comments
printf "\n"
echo $str6c
varTotc=$(find ../ACME/components/clm -name '*.F90' | xargs grep -o "\s$varIn\>" | wc -l)
echo $varTotc

#Return Definition count
printf "\n"
echo $str4
varDef=$(find ../ACME/components/clm -name '*.F90' | xargs sed 's/!.*//' | grep "\s$varIn\>" | grep -o "use" |wc -l)
echo $varDef

#Return Assignment count
printf "\n"
echo $str3
varAs=$(find ../ACME/components/clm -name '*.F90' | xargs sed 's/!.*//' | grep "\s$varIn\>" | grep -o "=>" |wc -l)
echo $varAs

#Return Association count
printf "\n"
echo $str5
varAsc1=$(find ../ACME/components/clm -name '*.F90' | xargs sed 's/!.*//' | grep "\s$varIn\>" | grep -o "\s=\s" |wc -l)
varAsc2=$(find ../ACME/components/clm -name '*.F90' | xargs sed 's/!.*//' | grep "\s$varIn\>" | grep -o "\s==\s" |wc -l)
varAsc=$(($varAsc1+$varAsc2))
echo $varAsc

#Return Other count
printf "\n"
echo $str7
varOth=$(($varTot-$varAs-$varDef-$varAsc))
echo $varOth




#****************************************Replace Variable**********************************************

#Give options
strA='Please select from the following: '
echo $strA
options=("Replace all variable instances" "Replace variable instances in select directory" "Replace variable instance by usage type" "Quit")
select opt in "${options[@]}"
do
  case $opt in
      "Replace all variable instances")
	      
              #Call for input
	      echo 'Enter new replacement variable'
	      read varReplace
              strB='Replacing variable:'
	      strC='with variable:'
	      strconfirm="$strB $varIn $strC $varReplace"
	      echo $strconfirm
	      
	      find ../ACME/components/clm -name '*.F90' | xargs sed -i "s/\b$varIn\b/$varReplace/g"

          ;;
      "Replace variable instances in select directory")
		#All top level Directories containing variable
		varDir=$(find ../ACME/components/clm -name '*.F90' | xargs grep -l "\s$varIn\>" | cut -d/ -f1-6 | sort -u)
		echo 'Enter new replacement variable'
		read varReplace
		
		for fold in $varDir; do
			printf "\n"
			echo $fold
			read -r -p "Would you like to replace all variables in this location? [y/N] " response
			if [[ "$response" =~ ^([yY][eE][sS]|[yY])+$ ]]
			then
				find $fold -name '*.F90' | xargs sed -i "s/\b$varIn\b/$varReplace/g"
				echo $varIn "replaced with" $varReplace
			else
				echo "The variable was not replaced for" $fold
			fi
			
		done

          ;;
      "Replace variable instance by usage type")
      		#Give options
		strA='Which type would you like to replace '
		echo $strA
		options=("Definition" "Association" "Assignment" "Other" "Quit")
		select opt in "${options[@]}"
		do
			case $opt in
	        		"Definition")
					#Call for input
					echo 'Enter new replacement variable'
			                read varReplace
					strB='Replacing variable:'
			                strC='with variable:'
			                strconfirm="$strB $varIn $strC $varReplace"
					echo $strconfirm
					#find ../ACME/components/clm -name '*.F90' | xargs sed -i "s/\b$varIn\b/$varReplace/g"
				;;
				"Association")
					echo "Association"
				;;
				"Assignment")
					echo "Assignment"
				;;
				"Other")
					echo "other"
				;;

              			"Quit")
					break
					;;
				*) echo invalid option;;
			esac
		done
          ;;
      "Quit")
          break
          ;;
      *) echo invalid option;;
   esac
done

#Replace based on directory

#Replace based on Type



#***********************************Match to Function name**********************************************
#find ../ACME/components/clm -name '*.F90' | xargs grep -n "\s$varIn\>" | cut -d : -f 1,2

#fileLine=$(find ../ACME/components/clm -name '*.F90' | xargs grep -n "\s$varIn\>" | cut -d : -f 1,2)
#echo $fileLine




#Look at all matches of input variable
#find . -name '*.F90' | xargs grep "\s$varIn\>"

#All matches with comments
#find ../ACME/components/clm -name '*.F90' | xargs sed 's/[^!]*!//' | grep "\s$varIn\>"



# "Return functions containing variable")

#              find ../ACME/components/clm -name '*.F90' | xargs grep -n "\s$varIn\>" | cut -d : -f 1,2
	      



