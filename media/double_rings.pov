#include "colors.inc"
#include "textures.inc"
#include "metals.inc"
#include "jerusalem.inc"
#include "ring.inc"
#include "analytical_g.inc" 
#declare FontSize=12;
#declare Xo = -0;
#declare Yo = 0;
#declare Zo = 8;
#declare Larrow = 30;
#declare CopperLz = 0.018;

camera {
    location <0,-120, 50>
    look_at  <0, 0, 0>
    rotate <0,0,25>

}


// box { <-1.5, 1, -3.5>, <0.5, -1, 1.5>
//      texture { Copper_Metal }}

// sphere {
//     <1, 0.5, 0>, 1
//    texture {Copper_Metal }
//     finish { phong 1
//              ambient 0.1 
//              brilliance 16.0
//              diffuse 0.9}
// }

//jerusalem(10, 12, 1, 1, 1, <0, 0, 0>, <0, 0, 0>)
//jerusalem(10, 12, 1, 1, 1, <0, 0, 0>, <0, 0, 90>)
#declare BPsize = 51;
box { <-BPsize, -BPsize, -0.018> <BPsize, BPsize, 0> 
       texture {Copper_Metal}
   }


#for (dy, -40, 40, 20)
#for (dx, -40, 40, 20)
ring(9.8, 8.3, CopperLz, <dx, dy, 2.00>, 90)
ring(5.1, 4.6, CopperLz, <dx, dy, 2.00>, 90)

box { <-10+dx, -10+dy, 0> <10+dx, 10+dy, 2> 
       texture {pigment {color Yellow transmit 0.2}}
   }

  
background{White}

#end
#end
plane {
     <0,0,1>, -CopperLz           //This represents the plane 0x+0y+z=0
     texture { pigment{ color White} }  //The texture comes from the file "metals.inc"
     finish { phong 1 reflection 1 }
   }
object{ Vector( <Xo,Yo,Zo>,<Xo,Yo,Larrow>, 1) pigment{ color Green }  }
#declare X = text{
     ttf "timrom.ttf" "X" 0.05, 0 
     pigment { color Black}
     finish { reflection 0 specular 1 }
           };
object{X rotate <90, 0,180>
         scale 10 
         translate<Xo-Larrow,Yo,Zo>
     } 
#declare Y = text{
     ttf "timrom.ttf" "Y" .05, 0 
     pigment { color Black}
     finish { reflection 0 specular 1 }
           };
object{Y rotate <90, 0,180>
         scale 10 
         translate<Xo+5,Yo+Larrow,Zo>
     }
#declare Z = text{
     ttf "timrom.ttf" "Z" .05, 0 
     pigment { color Black}
     finish { reflection 0 specular 1 }
           };
object{Z rotate <90, 0,180>
         scale 10 
         translate<Xo,Yo,Zo+Larrow-5>
     }

object{ Vector( <Xo,Yo,Zo>,<Xo,Yo+Larrow,Zo>, 1) pigment{ color Blue}  }
object{ Vector( <Xo,Yo,Zo>,<Xo-Larrow,Yo,Zo>, 1) pigment{ color Red }  }

light_source { <0, -50, 100> color White}


