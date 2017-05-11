#!/bin/bash

str1="Please enter the variable"
str2="is located in these parent directories"
str6="Total excluding comments:"

if [[ "$1" != "" ]]; then
	direct="$1"
else
	direct="../ACME/components/clm"
fi

varIn=${2:-grc}
rpl=$(echo $varIn|sed 's/$/_pp/')
varReplace=${3:-$rpl}
opt="Replace by usage type"

  case $opt in
      "Replace by usage type")
      		#Give options
		strA='Which type would you like to replace '
		echo $strA
		options=("Definition" "Association" "Assignment" "Replace All" "Quit")
		select opt1 in "${options[@]}"
		do
			case $opt1 in
	        		"Definition")
					find $direct -name '*.F90' | xargs sed -i "/\suse\s/s/\b$varIn\b/$varReplace/g"					
					strB='Replacing variable:'
					strC='with variable:'
					strconfirm="$strB $varIn $strC $varReplace"
					echo $strconfirm "when used as" $opt1
				;;
				"Association")
					find $direct -name '*.F90' | xargs sed -i "/=>/s/\b$varIn\b/$varReplace/g"
					strB='Replacing variable:'
					strC='with variable:'
					strconfirm="$strB $varIn $strC $varReplace"
					echo $strconfirm "when used as" $opt1
				;;
				"Assignment")
					find $direct -name '*.F90' | xargs sed -i "/\s==\s/s/\b$varIn\b/$varReplace/g"
					find $direct -name '*.F90' | xargs sed -i "/\s=\s/s/\b$varIn\b/$varReplace/g"
					strB='Replacing variable:'
					strC='with variable:'
					strconfirm="$strB $varIn $strC $varReplace"
					echo $strconfirm "when used as" $opt1
				;;
				"Replace All")
					find $direct -name '*.F90' | xargs sed -i "s/\b$varIn\b/$varReplace/g"
					strB='Replacing variable:'
					strC='with variable:'
					strconfirm="$strB $varIn $strC $varReplace"
					echo $strconfirm "when used as" $opt1
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

	      



