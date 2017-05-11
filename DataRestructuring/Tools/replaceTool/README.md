
General form
	./script.sh <directory> <initial variable> <replacement variable>

Defaults
	<directory>
        	"../ACME/components/clm"
		The '..' assumes the script is being run from the directory 'replaceTool'
	<initial variable>
		prev.sh and replace.sh:
		"grc"
		summary.sh:
		"col_pp"
	<replacement variable>
		prev.sh:
			initial variable + "_pp"
			In this case "grc_pp"
		replace.sh:
			initial variable or "grc"
		summary.sh
			not applicable



Examples:
	./summary.sh ../ACME/components/clm/src/biogeophys ecophyscon
		"Returns summary information in this format"
			<subdirectories containing variable>
			Number of total files containing variable

			Count of occurences

			Count of occurences for various usage types (Definition, Assignment, Association, Other and Comments)

		Example:
			../ACME/components/clm/src/biogeophys
			5 Files

			Total excluding comments:
			42

			Variable usage Types
			Definition  | Assignment | Association | Other | Comments
			       1           19           22         0         0	

	./prev.sh ../ACME/components/clm/src/biogeophys ecophyscon veg_pp
		"Returns preview of replacement in this format"
		/<Filename>.F90:<linenumber>
		< original line containing variable
		> replacement line containing variable

		Example:

		/PhotosynthesisMod.F90:271
		<         i_vcmax       => ecophyscon%i_vc       
		>         i_vcmax       => veg_pp%i_vc        
		/PhotosynthesisMod.F90:272
		<         s_vcmax       => ecophyscon%s_vc        
		>         s_vcmax       => veg_pp%s_vc           



	./replace.sh ../ACME/components/src/clm ecophyscon veg_pp

		Example:

			Which type would you like to replace
			1) Definition
			2) Association
			3) Assignment
			4) Replace All
			5) Quit
			#? 1
			Replacing variable: ecophyscon with variable: veg_pp when used as Definition
			#? 5	
	


	More examples:
	./prev.sh 
		dir = ../ACME/components/clm
		initial variable = grc
		replacement variable = grc_pp

	./replace.sh ../ACME/components/clm grc_pp grc	
