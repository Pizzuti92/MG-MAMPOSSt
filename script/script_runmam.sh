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
  echo "where the gomamposst_x.inp file is located."
  echo "-t | --time: print time of execution of the code"
  echo "-h | --help: displays help"
  echo "-ts | --test: perform an automated test for reproducibility"
  echo "***********************************"
	
}

cartellacl=$PWD/
printime=0
while [ -n "$1" ]; do
  case $1 in
    -d | --directory)
     cartellacl="$2"     
      ;;
    -h | --help) 
     Help
     exit;;
    -t | --time)
     printime=1 
     ;;
    -ts | --test)
     testt=1
     echo $testt
     ;;
  esac
  shift
done


start=`date +%s`
if [ $testt -ge 1 ]; then
 echo 'Performing basic smoke test on DHOST gravity.'
 echo 'This can take several minuts, please wait' 
 gomamposstopt < $cartellacl\gomamposst_x.inp  > test_record.txt
 python3 test.py
 echo 'records of the run written in test_record.txt'
 exit
else 

 gomamposstopt < $cartellacl\gomamposst_x.inp  
fi 
gomamposstopt < $cartellacl\gomamposst_x.inp  
 
end=`date +%s`

runtime=$((end-start))

if [ $printime -ge 1 ]; then
    echo "******************************"
	echo "time of execution: "  $runtime " seconds"
	echo "******************************"
fi	
python3 plot.py               

