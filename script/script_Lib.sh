#!/bin/bash
FCOMP=f95

Help()
{
  echo "***********************************"
  echo "script_Lib"
  echo
  echo "Install the necessary libraries and dependecies for MG-MAMPOSSt"
  echo
  echo "Syntax: ./scripts/script_Lib.sh [-f|-d|-h]"
  echo "Or: sh scripts/script_Lib.sh [-f|-d|-h]"
  echo "Options:"
  echo "-f | --fortrancompiler [f95 | gfort | ifort]: Default is f95" 
  echo "-d | --directory [<directory_path>]: select the working directory"
  echo "-h | --help: displays help"
  echo "***********************************"
	
}

 cartellacl=$PWD/

while [ -n "$1" ]; do
  case $1 in
    -f | --fortrancompiler) 
     FCOMP="$2"
     echo $FCOMP      
      ;;
    -d | --directory)
     cartellacl="$2" 
     ;; 
    -h | --help) 
     Help
     exit;;
  esac
  shift
done
  
 
 cd $cartellacl\Utili/

 $FCOMP -c *.f
 ar -r libUtil.a  *.o
 cd  ../GamI/
 
 
 
# $FCOMP -c $cartellacl\GamI/*.f
# $FCOMP -c $cartellacl\GamI/*.f90
# ar -r $cartellacl\GamI/libGAM.a $cartellacl\GamI/*.o

 $FCOMP -c *.f
 $FCOMP -c *.f90
 ar -r libGAM.a *.o

cd  ../Newuoa/
$FCOMP -c *.f
ar -r libNewuoa.a *.o

cd ../Powell/
$FCOMP -c *.f
ar -r libPowell.a *.o

cd ../JJin/
$FCOMP -c *.for
$FCOMP -c *.f
ar -r libJin.a *.o

cd ..
 
echo ''
echo '*******************************'
echo 'Libriaries created sucessfully'
echo '*******************************'

 
       


