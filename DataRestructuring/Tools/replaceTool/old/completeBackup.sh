#!/bin/bash

str1="Please enter the variable"
str2="is located in these locations"
str3="Count as an Assignment"
str4="Count as a Definition"
str5="Count as an Association"
str6="Total occurences excluding comments:"
str6c="Count as Comments:"
str7="Count as Other/Comments Type"

#Call for input
echo $str1
read varIn

#Confirm input variable
printf "\n"
str22="$varIn $str2"
echo $str22

#All top level Directories containing variable
find ../ACME/components/clm -name '*.F90' | xargs grep -l "\s$varIn\>" | cut -d/ -f1-6 |sort -u

#************************************Use Count*********************************************************
#Count all files containing variable
#find . -name '*.F90' | xargs grep -l "\s$varIn\>" | sort -u 
varLoc=$(find ../ACME/components/clm -name '*.F90' | xargs grep -l "\b$varIn\b" | sort -u | wc -l)
echo $varLoc "Files"

#Return total count excluding comments
printf "\n"
echo $str6
varTot=$(find ../ACME/components/clm -name '*.F90' | xargs sed 's/!.*//' | grep -o "\b$varIn\b" | wc -l)
echo $varTot

#Return total count of comments
varTotc=$(find ../ACME/components/clm -name '*.F90' | xargs grep -o "\b$varIn\b" | wc -l)
varCom=$(($varTotc-varTot))

#Return Definition count
varDef=$(find ../ACME/components/clm -name '*.F90' | xargs sed 's/!.*//' | grep "\b$varIn\b" | grep -o "\suse\s" |wc -l)

#Return Assignment count
varAs=$(find ../ACME/components/clm -name '*.F90' | xargs sed 's/!.*//' | grep "\b$varIn\b" | grep -o "=>" |wc -l)

#Return Association count
varAsc2=$(find ../ACME/components/clm -name '*.F90' | xargs sed 's/!.*//' | grep "\b$varIn\b" | grep -o "\s=\s" | wc -l)
varAsc3=$(find ../ACME/components/clm -name '*.F90' | xargs sed 's/!.*//' | grep "\b$varIn\b" | grep -o "\s==\s" |wc -l)
varAsc=$(($varAsc1+$varAsc2+$varAsc3))

#Return Other count
varOth=$(($varTot-$varAs-$varDef-$varAsc))


printf "\n"
echo "Variable usage Types"
echo "Definition  | Assignment | Association | Other"
echo "      " $varDef "         " $varAs "         " $varAsc "         " $varOth
printf "\n"
echo $str6c
echo $varCom

#Unhide to see all occurences
#find ../ACME/components/clm -name '*.F90' | xargs sed 's/!.*//' | grep "\b$varIn\b"
 

#****************************************Replace Variable**********************************************

#Give options
printf "\n"
strA='Please select from the following: '
echo $strA
options=("Replace all variable instances" "Replace by directory" "Replace by usage type" "Quit")
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

      "Replace by directory")
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

      "Replace by usage type")
      		#Give options
		strA='Which type would you like to replace '
		echo $strA
		options=("Definition" "Association" "Assignment" "Replace all remaining" "Quit")
		select opt in "${options[@]}"
		do
			case $opt in
	        		"Definition")
					#Call for input
					echo 'Enter new replacement variable'
			                read varReplace
					find ../ACME/components/clm -name '*.F90' | xargs sed -i "/\suse\s/s/\b$varIn\b/$varReplace/g"					
				;;
				"Association")
					#Call for input
					echo 'Enter new replacement variable'
					read varReplace
					find ../ACME/components/clm -name '*.F90' | xargs sed -i "/=>/s/\b$varIn\b/$varReplace/g"
				;;
				"Assignment")
					#Call for input
                                        echo 'Enter new replacement variable'
					read varReplace
					find ../ACME/components/clm -name '*.F90' | xargs sed -i "/\s==\s/s/\b$varIn\b/$varReplace/g"
					find ../ACME/components/clm -name '*.F90' | xargs sed -i "/\s=\s/s/\b$varIn\b/$varReplace/g"
				echo "Assignment"
				;;
				"Replace all remaining")
					#Call for input
					echo 'Enter new replacement variable'
					read varReplace
					find ../ACME/components/clm -name '*.F90' | xargs sed -i "s/\b$varIn\b/$varReplace/g"
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
	      



