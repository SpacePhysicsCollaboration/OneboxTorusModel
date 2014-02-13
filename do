#!/bin/bash

make all

if [ $? -eq 0 ] 
  then 

  time ./box > runlog

  if [ $? -eq 0 ] 
    then 
  
    gnuplot gplot
  
    display gplot.jpeg

  fi

fi

make clean
  
rm ft*.dat




