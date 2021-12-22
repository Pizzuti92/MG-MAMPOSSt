#!/bin/bash


Help()
{
  echo "***********************************"
  echo "script_runmam"
  echo
  echo "Run the MG-MAMPOSSt"
  echo
  echo "Input/Output files are selected in gomamposst_x.inp"
  echo "Additional Options in Options.txt"
  echo
  echo "Syntax: ./scripts/script_runmam.sh [-d|-h]"
  echo "Or: sh scripts/script_Lib.sh [-d|-h]"
  echo "Options:"
  echo "-d | --directory [<directory_path>]: select the working directory"
  echo "-h | --help: displays help"
  echo "***********************************"
	
}

cartellacl=$PWD/

while [ -n "$1" ]; do
  case $1 in
    -d | --directory)
     cartellacl="$2"     
      ;;
    -h | --help) 
     Help
     exit;;
  esac
  shift
done





./gomamposstopt.e < $cartellacl\gomamposst_x.inp   
python3 plot.py               

