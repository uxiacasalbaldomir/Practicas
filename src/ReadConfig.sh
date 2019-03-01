#!/bin/bash
# Read parameters from a configuration file
# Tamara Prieto 2017 at uvigo


if [ "$#" -ne 1 ]
then
        echo "You must specify absolute path to a configuration file as argument. More info about how to create this file: ReadConfig.sh --help"
	echo ""
#	echo "--------------------------------------------------------"
	exit
elif [ "$1" = "--help" ]
then
        echo ""
#        echo "------------------------------------------------------------------------------------------------------------------------"
        echo ""
        echo "INSTRUCTIONS IN HOW TO CREATE A CONFIG FILE"
	echo ""
        echo "The config file must contain following arguments and absolute directories."
	echo ""
	echo "s | sample_list: File containing one sample name per row. This file must be located in the resources directory."
        echo "c | control_list: File containing the name of the control samples."
        echo "ori | original_directory: Absolute path of the folder containing raw data and sample lists."
        echo "work | working_directory: Absolute path to the working directory."
        echo "ref | reference_name: Name of the reference genome."
        echo "res | resources_directory: Absolute path of the reference and other resources directory."
	echo "script | scripts_directory: Absolute path to the scripts directory"
	echo "ad | library_adapters: Tab delimited file containing the name of library in the first column and the sequence of the adapters in the other columns."
	echo "lib | library: Name of the sequencing library."
	echo "wga_lib | whole_genome_amplification_library: Name of the WGA kit among AMPLI-1,PICOPLEX or MALBAC."
        echo ""
#	echo "Example of lines in config file:"
#	echo ""
#	echo "sample_list=mysamples.txt"
#	echo "c=mycontrols.txt"
#	echo "resources_directory=/mnt/lustre/scratch/home/uvi/be/phylocancer/RESOURCES/"
#	echo ""
#	echo "These variables will correspond to:"
#        echo '$IDLIST,$SAMPLELIST,$CONTROL,$ORIDIR,$WORKDIR,$RESDIR,$SCRIPTDIR,$REF,$TARGET,$LIBRARY,$WGA_LIBRARY,$PLATFORM,$ADAPTERS,$QUEUE,$EMAIL,$GENDER'
#	echo "-----------------------------------------------------------------------------------------------------------------------"
        echo ""
        exit
else
	if [ -e $1 ]
	then
	while read i; do
		case $i in
		id=*|id_list=*)
		IDLIST="${i#*=}"
		;;
    		s=*|sample_list=*)
    		SAMPLELIST="${i#*=}"
    		;;
    		c=*|control_list=*)
    		CONTROL="${i#*=}"
    		;;
    		ori=*|original_directory=*)
    		ORIDIR="${i#*=}"
    		;;
    		work=*|working_directory=*)
    		WORKDIR="${i#*=}"
    		;;
    		res=*|resources_directory=*)
    		RESDIR="${i#*=}"
    		;;
    		script=*|scripts_directory=*)
		SCRIPTDIR="${i#*=}"
		;;
		ref=*|reference_name=*)
    		REF="${i#*=}"
    		;;
    		targ=*|target_capture=*)
    		TARGET="${i#*=}"
    		;;
                ad=*|library_adapters=*)
                ADAPTERS="${i#*=}"
                ;;
		lib=*|library=*)
		LIBRARY="${i#*=}"
		;;
		pl=*|platform=*)
		PLATFORM="${i#*=}"
		;;
		wga_lib=*|whole_genome_amplification_library=*)
                WGA_LIBRARY="${i#*=}"
                ;;
		e=*|email=*)
                EMAIL="${i#*=}"
                ;;
		p=*|queue=*)
		QUEUE="${i#*=}"
		;;
		g=*|gender=*)
		GENDER="${i#*=}"
		;;
    		*)
	    	echo "$i is not a valid parameter"
            	exit
    		;;
		esac
	done < $1
	else
	echo "File $1 is not in current directory or does not exist"
	exit
	fi
fi
