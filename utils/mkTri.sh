#!/bin/bash

awk 'BEGIN{
 dx = 0.2;
 lx = 10.0;
 ly = 10.0;
 lz = 1.0;
 print lx" "ly" "lz;
 for ( x=0.5*(dx-lx); x<=0.5*lx; x+=dx ) {
  for ( y=0.5*(dx-ly); y<=0.5*ly; y+=dx ) {
    for ( z=0.5*(dx-lz); z<=0.5*lz; z+=dx ) {
       if (x+0.5*lx<(1-(y/ly+0.5))*lx) {
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
cat tmp1 tmp >../Tri.xyz
rm tmp tmp1
