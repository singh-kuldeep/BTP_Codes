#include "vector.h"
#include <iostream>
using namespace std ; 

// void Point::offset (double offset_x, double offset_y)
// {
// 	x += offset_x;
// 	y +=offset_y ;
// }

// void Point::print (){
// 	cout<< "(" << x << "," << y << " )" << endl ; 
// }

void Vector::offset(double offset_x, double offset_y){
	start.offset(offset_x,offset_y);
	end.offset(offset_x,offset_y);	
}

void Vector::print(){
	start.print();
	cout<<" --> ";
	end.print() ;
	cout << endl ;
}

int main (){
	Vector vec ;
	Point p ; 

	// vec.start.xp = 1.0 ; 
	// vec.start.yp = 2.0 ;
	vec.end.xp = 3.0 ;
	vec.end.yp = 4.0 ;

	// p.xp = vec.end.xp ; 
	// p.yp = vec.end.yp ;

	vec.print ();
	p.print();
	return 0;
}