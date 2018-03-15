#include "colors.inc"
#include "textures.inc"
#include "metals.inc"
#include "jerusalem.inc"
#include "CutAbsorberPlate.inc"
#include "analytical_g.inc" 
#declare FontSize=12;
#declare Xo = -0;
#declare Yo = 0;
#declare Zo = 8;
#declare Larrow = 15;
#declare CopperLz = 0.018;
#declare Luc = 14;
#declare L1 = 4;
#declare L2 = 3;
#declare gap = 0.6;
#declare length = 2*L1+gap;
#declare L3 = sqrt(2)*(L1-L2) + gap/2;
#declare dphi = 45;
#declare reswidth = 0.6;
#declare reslength = 1;
#declare resheight = 0.2;
#declare rho= 2;
camera {
    location <0,-60, 50>
    look_at  <0, 0, 0>
    rotate <0,0,25>

}

#declare BPsize = 25;
box { <-BPsize, -BPsize, -0.018> <BPsize, BPsize, 0> 
       texture {Copper_Metal}
   }

// generate multiples of the unit-cell
#for (dy, -14, 14, 14)
#for (dx, -14, 14, 14)
box { <-7+dx, -7+dy, 0> <7+dx, 7+dy, 3.2> 
       texture {pigment {color Yellow transmit 0.2}}
   }

box { <dx-L3/2, dy-reswidth/2, 3.2> 
      <dx-L3/2-reslength,dy+reswidth/2,3.2+resheight>
       texture {pigment {color Black transmit 0.5}}
   }   
box { <dx+L3/2, dy-reswidth/2, 3.2> 
      <dx+L3/2+reslength,dy+reswidth/2,3.2+resheight>
       texture {pigment {color Black transmit 0.5}}
   }   
box { <dx-reswidth/2, dy-L3/2, 3.2> 
      <dx+reswidth/2,dy-L3/2-reslength,3.2+resheight>
       texture {pigment {color Black transmit 0.5}}
   }   
box { <dx-reswidth/2, dy+L3/2, 3.2> 
      <dx+reswidth/2,dy+L3/2+reslength,3.2+resheight>
       texture {pigment {color Black transmit 0.5}}
   }   
// outer set of resistors
box { <-reslength/2, -reswidth/2, 0> 
      <reslength/2,+reswidth/2,resheight>
       texture {pigment {color Black transmit 0.5}}
       rotate<0,0,-45>
       translate<dx+rho, dy+rho, 3.2>
   }   
box { <-reslength/2, -reswidth/2, 0> 
      <reslength/2,+reswidth/2,resheight>
       texture {pigment {color Black transmit 0.5}}
       rotate<0,0,135>
       translate<dx-rho, dy-rho, 3.2>
   }   
box { <-reslength/2, -reswidth/2, 0> 
      <reslength/2,+reswidth/2,resheight>
       texture {pigment {color Black transmit 0.5}}
       rotate<0,0,45>
       translate<dx+rho, dy-rho, 3.2>
   }   
box { <-reslength/2, -reswidth/2, 0> 
      <reslength/2,+reswidth/2,resheight>
       texture {pigment {color Black transmit 0.5}}
       rotate<0,0,225>
       translate<dx-rho, dy+rho, 3.2>
   }   
CutAbsorberPlate(4, 1, CopperLz, <dx, dy+L1/sqrt(2)+gap, 3.218>, 
                45)
CutAbsorberPlate(4, 1, CopperLz, <dx, dy-L1/sqrt(2)-gap, 3.218>, 
                225)
CutAbsorberPlate(4, 1, CopperLz, <dx+L1/sqrt(2)+gap, dy, 3.218>, 
                -45)
CutAbsorberPlate(4, 1, CopperLz, <dx-L1/sqrt(2)-gap, dy, 3.218>, 
                135)
box { <dx-L3/2, dy-L3/2, 3.218> 
      <dx+L3/2, dy+L3/2, 3.2> 
       texture { Copper_Metal finish {phong 1}
               }
   }
  
background{White}
#end
#end


plane {
     <0,0,1>, -CopperLz           //This represents the plane 0x+0y+z=0
     texture { pigment{ color White} }  //The texture comes from the file "metals.inc"
     finish { phong 1 reflection 1 }
   }
#declare X = text{
     ttf "timrom.ttf" "X" 0.05, 0 
     pigment { color Black}
     finish { reflection 0 specular 1 }
           };
object{X rotate <90, 0,180>
         scale 5 
         translate<Xo-Larrow,Yo,Zo>
     } 
#declare Y = text{
     ttf "timrom.ttf" "Y" .05, 0 
     pigment { color Black}
     finish { reflection 0 specular 1 }
           };
object{Y rotate <90, 0,180>
         scale 5 
         translate<Xo+5,Yo+Larrow,Zo>
     }
#declare Z = text{
     ttf "timrom.ttf" "Z" .05, 0 
     pigment { color Black}
     finish { reflection 0 specular 1 }
           };
object{Z rotate <90, 0,180>
         scale 5 
         translate<Xo,Yo,Zo+Larrow>
     }

object{ Vector( <Xo,Yo,Zo>,<Xo,Yo,Zo+Larrow>, 0.5) pigment{ color Green }  }
object{ Vector( <Xo,Yo,Zo>,<Xo,Yo+Larrow,Zo>, 0.5) pigment{ color Blue}  }
object{ Vector( <Xo,Yo,Zo>,<Xo-Larrow,Yo,Zo>, 0.5) pigment{ color Red }  }

light_source { <0, -50, 100> color White}


