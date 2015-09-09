class point
{
public:
	double x ;
	double y ;
	
};

class vector 
{
public:
	point start , end ;
};

#include <iostream>
using namespace std;

void printvetor ( vector v) 
{
	cout << "( " << v.start.x << " ," << v.start.y << " ) - > (" << v.end.x << " , " << v.end.y << " )" << endl ;	 
}

void offset (vector v, double offset_x , double offset_y )
{
	v.start.x += offset_x ;
	v.start.y += offset_y ;
	v.end.x += offset_x ;
	v.end.y += offset_y ; 
}

int main()
{
    vector  vec1 ;

    vec1.start.x = 1.2;
    vec1.start.y = 0.4 ;
    vec1.end.x = 2.0 ;
    vec1.end.y = 1.6 ;

 printvetor(vec1);
 cout << endl ;

 offset(vec1 , 1.0 , 1.5 ) ; 

 printvetor(vec1);
 cout << endl ;


cout << vec1.end.x << endl;
   
   return 0;
}