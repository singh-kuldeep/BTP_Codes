/*This call calulates the numerical diffusion flux, which is intern makes the
solver second order accureate in the space*/

#include "math.h"
#include "iostream"
#include "interface.h"

using namespace std ;

class diffusion_flux
{
	public:
	float diffusionfluxvector[5] ;

	diffusion_flux( vector<float>& U_l_m, vector<float>& U_l, vector<float>& 
		U_r, vector<float>& U_r_p, vector<float>& A_l, vector<float>& A_r,
		vector<float>& A_r_p, float V_l_m, float V_l, float V_r, float V_r_p,
		float dt)
	{
		float thetai[5];
		float thetaiplus[5];
		float gvector[5] ;
		float gvectorplus[5] ;
		float signale ;
		float betainterface[5] ;
		float phiinterface[5] ;
		float Z_interface[5];
		float pshiinterface[5] ;
	
		interface left(U_l_m, U_l, A_l, V_l_m, V_l, dt) ;
		interface right(U_l, U_r, A_r, V_l, V_r, dt) ;
		interface rightplus(U_r, U_r_p, A_r_p, V_r, V_r_p, dt);

		// here I am somehow trying to defie the gi and gi+1
		for (int i = 0; i < 5; ++i)
		{
			if (right.g_interface[i] >= 0)
			{
				signale = 1.0;
			}
			else
			{
				signale = -1.0 ;
			}

			// this is pura jugaad so change karna hai baad me 
			//"gvector[i]" ka expression
			float tempB = (left.g_interface[i]*signale)  ;
			float tempA = fabs(right.g_interface[i]) ;
			tempA = min(tempA,tempB) ;
			float temp_zero=0.0;
			tempB = max(temp_zero, tempA) ;
			gvector[i] = signale *  tempB ;

			// gvector[i] = signale * max(0.0, min(fabs(right.g_interface[i]),
			// 		left.g_interface[i]*signale ) ) ;
		}

		for (int i = 0; i < 5; ++i)
		{
			if (rightplus.g_interface[i] >= 0)
			{
				signale = 1.0;
			}
			else
			{
				signale = -1.0 ;
			}

			// this is pura jugaad so change karna hai baad me 
			//"gvector[i]" ka expression
			float tempA = fabs(rightplus.g_interface[i]) ;
			float tempB = (right.g_interface[i]*signale)  ;
			tempA = min(tempA,tempB) ;
			float temp_zero=0.0;
			tempB = max(temp_zero, tempA) ;
			gvectorplus[i] = signale * tempB ;
			// gvector[i] = signale * max(0.0, min( fabs
			// (right.g_interface[i]),
			// left.g_interface[i]*signale ) ) ;
		}

		// theta i term defying 
		for (int i = 0; i < 5; ++i)
			{
				//here I have doubt about not equal to symbol because it cant 
				//be exect 0.00000 so most of the time 
				//we endup choosing theta i = 0.0l
				if((fabs(right.alpha_interface[i]) + 
					fabs(left.alpha_interface[i]))!=0.0)
				{
					thetai[i] = fabs(right.alpha_interface[i] - left.
						alpha_interface[i])/(fabs(right.alpha_interface[i]) +
						fabs(left.alpha_interface[i])) ;
				}
				else
				{
					thetai[i] = 0.0;
				}
			}

		// theta i plus 
		for (int i = 0; i < 5; ++i)
		{
			/* here I have doubt about not equal to symbol because it can't 
			be exect 0.00000 so most of the
			time we endup choosing theta i = 0.0 */

			if((fabs(rightplus.alpha_interface[i])+
				fabs(right.alpha_interface[i]))!=0.0)
			{
				thetaiplus[i]=fabs(rightplus.alpha_interface[i]-right.
					alpha_interface[i])/(fabs(rightplus.alpha_interface[i]) +
					fabs(right.alpha_interface[i])) ;
			}
			else
			{
				thetaiplus[i] = 0.0 ;
			}
		}

		// beta defying
		// before that we nee to define the omega values which are constant 
		float omega[5] = {0.25, 1.0, 1.0, 1.0, 0.25} ;
		for (int i = 0; i < 5; ++i)
		{
			betainterface[i] = 1.0 + omega[i]*(max(thetai[i],thetaiplus[i]));
		}
		/////////////////////////////////
		// phi defying
		for (int i = 0; i < 5; ++i)
		{
			if (right.alpha_interface[i] != 0)
			{
				phiinterface[i] = (gvectorplus[i] - gvector[i]) / 
				right.alpha_interface[i] ;
			}
			else
			{
				phiinterface[i] = 0.0 ;
			}
		}

		// now re-defying the Zinterface
		for (int i = 0; i < 5; ++i)
		{
			Z_interface[i] = right.Z_interface[i] + betainterface[i]*
			phiinterface[i] ;
		}

		// // pshiinterface re-defying (for invisid wall flow)
		// for (int i = 0; i < 5; ++i)
		// {
		// 	pshiinterface[i] = pow(Z_interface[i],2) + 0.25 ;
		// }

		// re-defying the pshiinterface (once it is defined in the interface 
		//class change that also)
		// again we need to specify the deltaf = 0.2
			float deltaf = 0.2 ;
			for (int i = 0; i < 5; ++i)
			{
				if (Z_interface[i] >= deltaf)
				{
					pshiinterface[i] = fabs(Z_interface[i]) ;
				}
				else
				{
					pshiinterface[i] = 0.5*(pow(Z_interface[i],2) + pow(
						deltaf,2))/ deltaf ;
				}
			}

		// finally calculating the interface diffusionterm=diffusionfluxvetor
		// before that we can calculate tempvector as
		float tempvector[5];
		for (int i = 0; i < 5; ++i)
		{
			tempvector[i] = ( betainterface[i]*(gvector[i]+gvectorplus[i]) -
				pshiinterface[i]*right.alpha_interface[i] )/dt ;
		}

		// now the diffusionfluxvector
		for (int i = 0; i < 5; ++i)
		{
			diffusionfluxvector[i] = 0.0 ;
			for (int l = 0; l < 5; ++l)
			{
				diffusionfluxvector[i] += right.eigenvectormatrix[i][l] *
				tempvector[l] ;
			}
		}
	};
		// ~diffusion_flux();
	
};