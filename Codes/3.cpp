/*
 * @Time   : 2018/09/19
 * @Author : Iydon
 * @File   : 3.cpp
*/
// include
#include <iostream>
#include <math.h>
// define
#define in(a,f,b)   a<=f && f<=b
#define len(A)      sizeof(A)/sizeof(*A)
#define lambda(x,f) [](double x)->double{return f;}
// namespace
using namespace std;

// Statement
double fixed_point(auto, double, int, int);
double Newton_method(auto, auto, double, int=32, int=6);
template<typename T>void print(const T&v){cout<<v<<"\n";}
template<typename T, typename... Args>void  print(const T&x,const Args&...as){
	cout << x << ", "; print(as...);
}

int main(){
	auto f  = lambda(x, cos(x));
	auto df = lambda(x, -sin(x));
	double zero = Newton_method(f, df, 0.7);
	print(f(zero),zero);
}



double fixed_point(auto f, double start, int max_step, int sign_dig){
	double eps = pow(10, -sign_dig);
	double new_val = f(start);
	double old_val = 0.0;
	for (int i=0; i<max_step; i++){
		old_val = new_val;
		new_val = f(old_val);
		if (in(-eps, old_val-new_val, eps))
			return new_val;
	}
	return 0.0;
}

double Newton_method(auto f, auto df, double start, int max_step, int sign_dig){
	auto fun = [f,df](double x)->double{return x-f(x)/df(x);};
	return fixed_point(fun, start, max_step, sign_dig);
}
