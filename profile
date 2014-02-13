#!/bin/bash

make all

if [ $? -eq 0 ] 
  then 

  valgrind --tool=callgrind ./box

  if [ $? -eq 0 ] 
    then 

    mv callgrind.out.* profile.dat

    kcachegrind profile.dat 

  fi

fi

make clean
  
rm ft*.dat




