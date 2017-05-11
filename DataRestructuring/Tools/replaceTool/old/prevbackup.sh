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
opt="Preview replacement by type"

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
		select opt1 in "${options[@]}"
		do
			case $opt1 in
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
					#find $direct -name '*.F90' | xargs sed -i "s/\b$varIn\b/$varReplace/g"
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
		echo $strA $varIn $varReplace
		options=("Definition" "Association" "Assignment" "Replace all" "Quit")
		select opt1 in "${options[@]}"
			do
			case $opt1 in
				"Definition")
					find $direct -name '*.F90' | xargs perl -e 'for my $f (@ARGV){ my $f1=$f;$f1=~s|.*'$direct'||;$. = 0; open A, "$f"; while (<A>) { if ($_=~ m/\suse\s/){if ($_=~m/\b'$varIn'\b/){ print "$f1:$.\n\<$_\>"; $_=~s/\b'$varIn'\b/'$varReplace'/g;print $_;}}}}' | less
				;;
				"Association")					
					#find $direct -name '*.F90' | xargs sed -n "/=>/s/\b$varIn\b/$varReplace/gp" | less	
					find $direct -name '*.F90' | xargs perl -e 'for my $f (@ARGV){ my $f1=$f;$f1=~s|.*'$direct'||;$. = 0; open A, "$f"; while (<A>) { if ($_=~ m/=>/){if ($_=~m/\b'$varIn'\b/){ print "$f1:$.\n\<$_\>"; $_=~s/\b'$varIn'\b/'$varReplace'/g;print $_;}}}}' | less
				;;
				"Assignment")
					#find $direct -name '*.F90' | xargs sed -n "/\s==\s/s/\b$varIn\b/$varReplace/gp" | less
					#find $direct -name '*.F90' | xargs sed -n "/\s=\s/s/\b$varIn\b/$varReplace/gp" | less
					find $direct -name '*.F90' | xargs perl -e 'for my $f (@ARGV){ my $f1=$f;$f1=~s|.*'$direct'||;$. = 0; open A, "$f"; while (<A>) { if ($_=~ m/\s(==|=)\s/){if ($_=~m/\b'$varIn'\b/){ print "$f1:$.\n\<$_\>"; $_=~s/\b'$varIn'\b/'$varReplace'/g;print $_;}}}}' | less
				;;
				"Replace all")
					#find $direct -name '*.F90' | xargs sed -n "s/\b$varIn\b/$varReplace/gp" | less					
					find $direct -name '*.F90' | xargs perl -e 'for my $f (@ARGV){ my $f1=$f;$f1=~s|.*'$direct'||;$. = 0; open A, "$f"; while (<A>) {if ($_=~m/\b'$varIn'\b/){ print "$f1:$.\n\<$_\>"; $_=~s/\b'$varIn'\b/'$varReplace'/g;print $_;}}}' | less
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

	      



