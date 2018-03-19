#ifndef(plotball_Inc)
    #declare plotball_Inc = version;
    #version 3.1;
#macro plotball(X,Y,Z,R,period)
 sphere{<X,Y,Z>,R
   texture{
     pigment{
       gradient <0,1,0>
       color_map {[0 color rgb<1,1,1>]  // White
                  [0.25 color rgb<1,0,0>] // Red
                  [0.5 color rgb <1,1,1>]
                  [0.75 color rgb <0,0,1>]
                  [1 color rgb<1,1,1>]
                     
              }
        scale period
             }
    }// end of texture
 }// end of sphere
       
 #end // macro end
#version plotball_Inc;
#end

