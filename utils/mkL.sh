#!/bin/bash

awk 'BEGIN{
 dx = 0.5;
 lx = 10.0;
 ly = 16.0;
 lz = 2.0;
 w  = 2.0; # width
 print lx" "ly" "lz;
 # upper part of L
 for ( x=0.5*(dx-w); x<=0.5*w; x+=dx ) {
  for ( y=0.5*(dx-ly); y<=0.5*ly-w; y+=dx ) {
    for ( z=0.5*(dx-lz); z<=0.5*lz; z+=dx ) {
            sx = (rand()-0.5)*dx*0.1; # 10% random perturbation
            sy = (rand()-0.5)*dx*0.1; # 10% random perturbation
            sz = (rand()-0.5)*dx*0.1; # 10% random perturbation
            print "Cu "x+sx-0.5*lx+0.5*w" "y+sy+w" "z+sz
    }
  }
 }
 # lower part of L
 for ( x=0.5*(dx-lx); x<=0.5*lx; x+=dx ) {
  for ( y=0.5*(dx-w); y<=0.5*w; y+=dx ) {
    for ( z=0.5*(dx-lz); z<=0.5*lz; z+=dx ) {
            sx = (rand()-0.5)*dx*0.1; # 10% random perturbation
            sy = (rand()-0.5)*dx*0.1; # 10% random perturbation
            sz = (rand()-0.5)*dx*0.1; # 10% random perturbation
            print "Cu "x+sx" "y+sy-0.5*ly+0.5*w" "z+sz
    }
  }
 }
}' >tmp

wc -l tmp | awk '{print $1-1}' >tmp1
#echo "" >>tmp1
cat tmp1 tmp >../L.xyz
rm tmp tmp1
