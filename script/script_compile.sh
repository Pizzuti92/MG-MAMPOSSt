#!/bin/bash

Help()
{
  echo "***********************************"
  echo "script_compile"
  echo
  echo "Compile the MG-MAMPOSSt source code"
  echo
  echo "Syntax: ./scripts/script_compile.sh [-f|-h]"
  echo "Or: sh scripts/script_compile.sh [-f|-h]"
  echo "Options:"
  echo "-f | --fortrancompiler [f95 | gfort | ifort]: Default is f95" 
  echo "-h | --help: displays help"
  echo "***********************************"
	
}





FCOMP=f95

while [ -n "$1" ]; do
  case $1 in
    -f | --fortrancompiler) 
     FCOMP="$2"
     echo $FCOMP      
      ;;
    -h | --help) 
     Help
     exit;;
  esac
  shift
done



$FCOMP -o gomamposstopt.e gomamposstoptS.f -L GamI/ -lGAM -L Utili/ -lUtil -L Newuoa/ -lNewuoa -L Powell/ -lPowell -L JJin/ -lJin  -lm
