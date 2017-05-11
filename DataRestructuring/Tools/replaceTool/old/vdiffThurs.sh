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
echo "      " $varDef "         " $varAs "         " $varAsc "         " $varOth "       " $varCom


#Unhide to see all occurences
#find ../ACME/components/clm -name '*.F90' | xargs sed 's/!.*//' | grep "\b$varIn\b"
 

#****************************************Replace Variable**********************************************

#Give options
printf "\n"
strA='Please select from the following: '
echo $strA
options=("Replace all" "Replace by usage type" "Preview replacement by type" "Quit")
select opt in "${options[@]}"
do
  case $opt in
      "Replace all")
	      
              #Call for input
	      echo 'Enter new replacement variable'
	      read varReplace
              strB='Replacing variable:'
	      strC='with variable:'
	      strconfirm="$strB $varIn $strC $varReplace"
	      echo $strconfirm
	      
	      find $direct -name '*.F90' | xargs sed -i "s/\b$varIn\b/$varReplace/g"

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
					find $direct -name '*.F90' | xargs sed -i "/\suse\s/s/\b$varIn\b/$varReplace/g"					
				;;
				"Association")
					#Call for input
					echo 'Enter new replacement variable'
					read varReplace
					find $direct -name '*.F90' | xargs sed -i "/=>/s/\b$varIn\b/$varReplace/g"
				;;
				"Assignment")
					#Call for input
                                        echo 'Enter new replacement variable'
					read varReplace
					find $direct -name '*.F90' | xargs sed -i "/\s==\s/s/\b$varIn\b/$varReplace/g"
					find $direct -name '*.F90' | xargs sed -i "/\s=\s/s/\b$varIn\b/$varReplace/g"
				echo "Assignment"
				;;
				"Replace all remaining")
					#Call for input
					echo 'Enter new replacement variable'
					read varReplace
					find $direct -name '*.F90' | xargs sed -i "s/\b$varIn\b/$varReplace/g"
				;;

              			"Quit")
					break
					;;
				*) echo invalid option;;
			esac
		done
          ;;

      "Preview replacement by type")
		#Give options
		strA='Which type would you like to see a preview'
		echo $strA
		options=("Definition" "Association" "Assignment" "Replace all" "Quit")
		select opt in "${options[@]}"
			do
			case $opt in
				"Definition")
					#Call for input
					echo 'Enter new replacement variable'
					read varReplace
					find $direct -name '*.F90' | xargs sed -n "/\suse\s/s/\b$varIn\b/$varReplace/gp" | less
				;;
				"Association")					
					#Call for input
					echo 'Enter new replacement variable'
					read varReplace
					find $direct -name '*.F90' | xargs sed -n "/=>/s/\b$varIn\b/$varReplace/gp" | less	
				;;
				"Assignment")
					#Call for input
                                      	echo 'Enter new replacement variable'
					read varReplace
					find $direct -name '*.F90' | xargs sed -n "/\s==\s/s/\b$varIn\b/$varReplace/gp" | less
					find $direct -name '*.F90' | xargs sed -n "/\s=\s/s/\b$varIn\b/$varReplace/gp" | less
				;;
				"Replace all")
					#Call for input
					echo 'Enter new replacement variable'
					read varReplace
					find $direct -name '*.F90' | xargs sed -n "s/\b$varIn\b/$varReplace/gp" | less					
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
	      



