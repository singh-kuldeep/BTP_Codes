// both the convection and diffusion 
#include "math.h"
class delta_t
{
public:
	float dt ; 
	delta_t(vector<float> variablesvector, float delta_S, float CFL){
		float net_velocity = sqrt(pow(variablesvector[1],2) + pow(
			variablesvector[2],2) + pow(variablesvector[3],2))/
			variablesvector[0] ;
		float sound_speed = sqrt((0.56*variablesvector[4]/variablesvector[0])+
			0.56*0.5*pow(net_velocity,2)) ; 
		float max_char_speed =  net_velocity + sound_speed ;  
		dt =  CFL * delta_S/max_char_speed;
	};

	// ~delta_t();
	
};