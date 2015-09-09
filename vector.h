// #include <iostrem>
#include <iostream>
using namespace std ; 
class Point{
public:
	double xp,yp ;
	Point (){
		xp = 5.0; yp = 10.0; 
		cout << "Point instance created" << std::endl ;  
	}
	void offset(double offset_x, double offset_y){
		xp += offset_x ;
	    yp += offset_y ;
	}
	void print (){
		cout<< "(" << xp << "," << yp << " )" << std::endl ;
	}
};

class Vector{
public:
	Point start,end ;
	void offset(double offset_x, double offset_y);
	void print ();
};