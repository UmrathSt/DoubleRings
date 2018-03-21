#include "colors.inc"
#include "shapes.inc"
#include "shapes3.inc"
#include "textures.inc"
#include "glass.inc"
#include "ball.inc"
#include "analytical_g.inc"
#include "ball.inc"
#declare arrowl = 12.5;
#declare L = 10;
#declare alph = 25;
#declare distfact = .5;
#declare Z = L*1.2;  // start value Z
#declare EndZ = 40;  //   end value Z
#declare Step = 0.1;// step value
// loop start Z:
#while ( Z < EndZ + Step)

  #declare X = -10;  // start value X
  #declare EndX = 20;//   end value X
  //loop start X:
  #while ( X < EndX + Step)

 plotball(X, sin(2*Z), Z, 0.15, 4) 

  #declare X = X + Step;//next X
  #end // ----------- loop end X

#declare Z = Z + Step; //next Z
#end // ------------ loop end Z
object{
    cylinder { <0,L/2,0>,<0,L/2,L/sqrt(2)>,0.1
            color Black 
            rotate <0,alph,0>
            }
     }
object{
    cylinder { <0,L/2,0>,<0,L/2,L/sqrt(2)>,.1
    //        pigment { checker
    //                color White
    //                color Red}
        pigment {color Red    }
        rotate <0,0,0>
             }
      }
object{ Segment_of_Torus( L/sqrt(2)*0.9, 0.15, alph)
    texture { pigment{color Black}
              finish { phong 1 }  
            } // end of texture

            rotate<0,-90,0>
            translate <0,L/2,0>
} // end of Segment_of_Torus(...) ----


background{White}
camera {
    location <30, 30, 15>
    look_at <0, 0, 0>
//    focal_point 0 blur_samples 50 aperture .2
}
                                         
light_source {
    <5, 20, 25>
    color White 
}



object{
Round_Box (
    <-L/2, -L/2, -L/2>,
    <L/2, L/2, L/2>, 0.15, 0)
    pigment {Gray transmit 0.1}
finish{
        conserve_energy
        diffuse 1
        ambient 0.2 
        specular 0.85
        reflection{ 0 0.63 fresnel metallic off}
        roughness 0.001
        phong 2
        phong_size 100
    }
    rotate <0,alph,0>
    interior{ ior 1.32}
}
#declare X = text{
 ttf "timrom.ttf" "X" .05, 0
 pigment { color Black}
 finish { reflection 0 specular 1 }
       };
#declare Y = text{
 ttf "timrom.ttf" "Y" .05, 0
 pigment { color Black}
 finish { reflection 0 specular 1 }
       };

#declare Z = text{
 ttf "timrom.ttf" "Z" .05, 0
 pigment { color Black}
 finish { reflection 0 specular 1 }
       };


object{X rotate <0, 40,0>
         scale 2 
         translate<arrowl,distfact*L,0>
     }
object{Y rotate <0, 40,0>
         scale 2 
         translate<0,distfact*L+arrowl,0>
     }
object{Z rotate <180, 40,0>
         scale 2 
         translate<0,distfact*L+2,-arrowl>
     }
object{ Vector( <0,distfact*L,0>,
                <arrowl,distfact*L,0>, 0.25) pigment{ 
        color Green }  }

object{ Vector( <0,distfact*L,0>,
                        <0,distfact*L+arrowl,0>, 0.25) pigment{ 
        color Blue }  }

object{ Vector( <0,distfact*L,0>,
                <0,distfact*L,-arrowl>, 0.25) pigment{ 
        color Red }  }



plane {
    <0, 1, 0>,L*1000 
    pigment {
    White 
    }
    finish {
        specular 0.9 
        roughness 0.1 
        ambient 0.15 
        diffuse 3.5 
    }

}

