// POV-Ray include file "jerusalem.inc"
//-------------------------------------
#ifndef( jerusalem_Inc)
    #declare jerusalem_Inc = version;
    #version 3.1;
    //-------------------------------------
#macro jerusalem( lx,  // column height
                  ly,
                  lz,  // base radius
                  wx,// border radius
                  wy,  // top scale
                  origin,
                  rotation  
                )
union{ box{ <lx/2, ly/2,-lz/2>, <-lx/2, ly/2+wy, lz/2>
            texture {Copper_Metal}
            translate origin
            rotate rotation}
       box{ <lx/2, -ly/2, -lz/2>, <-lx/2, -ly/2-wy, lz/2>
            texture {Copper_Metal}
            translate origin
            rotate rotation}
       box{ <-wx/2, -ly/2, lz/2>, <wx/2, ly/2, -lz/2>
            texture {Copper_Metal}
            translate origin
            rotate rotation}
     }// end Column
 #end// ------------------- end of macro
 
#version jerusalem_Inc;
#end 
                                                                                          // ---------------- end of include file 
