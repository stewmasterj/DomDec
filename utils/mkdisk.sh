#!/bin/bash

awk 'BEGIN{
 dx = 0.5;
 lx = 10.0;
 ly = 10.0;
 lz = 2.0;
 print lx" "ly" "lz;
 for ( x=0.5*(dx-lx); x<=0.5*lx; x+=dx ) {
  for ( y=0.5*(dx-ly); y<=0.5*ly; y+=dx ) {
    for ( z=0.5*(dx-lz); z<=0.5*lz; z+=dx ) {
       r=sqrt(x**2+y**2); 
       if (r<5.1) {
            sx = (rand()-0.5)*dx*0.1; # 10% random perturbation
            sy = (rand()-0.5)*dx*0.1; # 10% random perturbation
            sz = (rand()-0.5)*dx*0.1; # 10% random perturbation
            print "Cu "x+sx" "y+sy" "z+sz
       }
    }
  }
 }
}' >tmp

wc -l tmp | awk '{print $1-1}' >tmp1
#echo "" >>tmp1
cat tmp1 tmp >../disk.xyz
rm tmp tmp1
