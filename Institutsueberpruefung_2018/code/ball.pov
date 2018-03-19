#include "colors.inc"
#include "ball.inc"


#declare Z = -5;  // start value Z
#declare EndZ = 5;  //   end value Z
#declare Step = 0.1;// step value
// loop start Z:
#while ( Z < EndZ + Step)

  #declare X = -5;  // start value X
  #declare EndX = 5;//   end value X
  //loop start X:
  #while ( X < EndX + Step)

 plotball(X, sin(2*X), Z, 0.15, 4) 

  #declare X = X + Step;//next X
  #end // ----------- loop end X

#declare Z = Z + Step; //next Z
#end // ------------ loop end Z

background{White}
camera {
     location <-1,4, 12>
     look_at <0, 1, 0>
}
 
 
 
 
light_source {
     <-50, 40, 40>
     color White
        }

