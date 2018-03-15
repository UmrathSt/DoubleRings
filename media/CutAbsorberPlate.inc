// POV-Ray include file "jerusalem.inc"
//-------------------------------------
#include "shapes3.inc"
#ifndef( CutAbsorberPlate_Inc)
    #declare CutAbsorberPlate_Inc = version;
    #version 3.1;
    //-------------------------------------
#macro CutAbsorberPlate(L1,  // the plates 1st edge length
                        L2,  // the plates 2nd edge length which is 
                        lz,  // supposed to be smaller than L1
                        origin, 
                        rotation)  
object{ prism{0, lz, 6,
             <-L2/2, -L1/2>,
             <+L1/2, -L1/2>,
             <+L1/2, +L1/2>,
             <-L1/2, +L1/2>,
             <-L1/2, -L2/2>,
             <-L2/2, -L1/2>
        texture{ Copper_Metal finish {phong 1}
               }
        rotate<-90,0,0>
        rotate<0, 0, rotation> 
        translate origin
            }
        }
          
 #end// ------------------- end of macro
 
#version CutAbsorberPlate_Inc;
#end 
                                                                                          // ---------------- end of include file 
