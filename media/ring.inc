// POV-Ray include file "jerusalem.inc"
//-------------------------------------
#include "shapes3.inc"
#ifndef( ring_Inc)
    #declare ring_Inc = version;
    #version 3.1;
    //-------------------------------------
#macro ring( Ri,  // column height
             Ra,
             lz,  // base radius
             origin,
             rotation)  
object{ Segment_of_CylinderRing(Ra, Ri, lz, 360)
        texture{ Copper_Metal finish {phong 1}
               }
        rotate<rotation,0,0>
        translate origin 
      }
          
 #end// ------------------- end of macro
 
#version ring_Inc;
#end 
                                                                                          // ---------------- end of include file 
