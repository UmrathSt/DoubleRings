#include "colors.inc"
#include "shapes.inc"
#include "textures.inc"
#include "glass.inc"
#include "ball.inc"

#include "ball.inc"

#declare L = 10;
#declare Z = L*1.4;  // start value Z
#declare EndZ = 40;  //   end value Z
#declare Step = 0.1;// step value
// loop start Z:
#while ( Z < EndZ + Step)

  #declare X = -5;  // start value X
  #declare EndX = 5;//   end value X
  //loop start X:
  #while ( X < EndX + Step)

 plotball(X, sin(2*Z), Z, 0.15, 4) 

  #declare X = X + Step;//next X
  #end // ----------- loop end X

#declare Z = Z + Step; //next Z
#end // ------------ loop end Z
union {
	cylinder { <0,2,20>,<0,2, L>,0.15 }
	cone { <0,2,L-0.75>,0,<0,2,L>,0.5 }
	pigment { color Black }
}


background{White}
camera {
    location <20, 20, 25>
    look_at <0, 0, 0>
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
    rotate <0,45,0>
    interior{ ior 1.32}
}



plane {
    <0, 1, 0>,-L/2 
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

