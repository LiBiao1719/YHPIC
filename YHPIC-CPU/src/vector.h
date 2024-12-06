#ifndef	__VECTOR_H
#define	__VECTOR_H

#include	<math.h>
#include "otypen.h"

class Vector3
{
	Scalar	x1, x2, x3;

public:
	Vector3(Scalar _x1=0.0, Scalar _x2=0.0, Scalar _x3=0.0) {x1=_x1; x2=_x2;
		x3=_x3;};
	Scalar	e1() const {return x1;};
	Scalar	e2() const {return x2;};
	Scalar	e3() const {return x3;};
	Scalar  e(int component){
		switch (component)
			{
			case 1:
				return x1;
			case 2:
				return x2;
			case 3:
				return x3;
			}
		return (0.0);
	}
	void	set_e1(Scalar _x1) {x1 = _x1;};
	void	set_e2(Scalar _x2) {x2 = _x2;};
	void	set_e3(Scalar _x3) {x3 = _x3;};
	void set(int component, Scalar value){
		switch (component)
			{
			case 1:
				x1 = value;
				break;
			case 2:
				x2 = value;
				break;
			case 3:
				x3 = value;
				break;
			}
	}
	Vector3 operator + (const Vector3& v) {return Vector3(x1+v.x1, x2+v.x2, x3+v.x3);};
	Vector3	operator - (const Vector3& v) {return Vector3(x1-v.x1, x2-v.x2, x3-v.x3);};
	Vector3	operator * (Scalar c) {return Vector3(c*x1, c*x2, c*x3);};
	Vector3	operator / (Scalar c) {Scalar ic = 1/c; return Vector3(ic*x1, ic*x2,
		ic*x3);};
	Scalar	operator * (const Vector3& v) {return x1*v.x1 + x2*v.x2 + x3*v.x3;};
	Vector3&	operator += (Scalar c) {x1 += c; x2 += c; x3 += c; return *this;};
	Vector3&	operator -= (Scalar c) {x1 -= c; x2 -= c; x3 -= c; return *this;};
	Vector3&	operator *= (Scalar c) {x1 *= c; x2 *= c; x3 *= c; return *this;};
	Vector3&	operator *= (const Vector3& v) {x1 *= v.x1; x2 *= v.x2; x3 *= v.x3;
		return *this;};
	Vector3&	operator += (const Vector3& v) {x1 += v.x1; x2 += v.x2; x3 += v.x3;
		return *this;};
	Vector3&	operator -= (const Vector3& v) {x1 -= v.x1; x2 -= v.x2; x3 -= v.x3;
		return *this;};
	Vector3& operator /= (Scalar c) {Scalar ic = 1/c; return *this *= ic;}
	Vector3	cross(const Vector3& v){return Vector3(x2*v.x3 - x3*v.x2, x3*v.x1
		- x1*v.x3, x1*v.x2 - x2*v.x1);};
	Vector3	jvMult(const Vector3& v) {return Vector3(x1*v.x1, x2*v.x2, x3*v.x3);};
	Vector3	scaleBy(const Vector3 v) {return jvMult(v);};
	Vector3	jvDivide(const Vector3& v) {return Vector3(x1/v.x1, x2/v.x2, x3/v.x3);}; //added by kc 8-30-94
	Scalar	magnitude() {return sqrt(x1*x1 + x2*x2 + x3*x3);};
	Scalar  energy() { return x1*x1+x2*x2+x3*x3;};
	Vector3	ComponentSqrt(){return Vector3(sqrt(x1), sqrt(x2),	sqrt(x3));}
	Vector3 ComponentSqr(){return Vector3(sqr(x1), sqr(x2), sqr(x3));}
	int	isNonZero() {if (x1 == 0.0 && x2 == 0.0 && x3 == 0) return 0;
		else return 1;}
};

inline Vector3	operator *(Scalar c,const Vector3& v)
{
	return Vector3(c*v.e1(), c*v.e2(), c*v.e3());
}


class Vector2
{
	Scalar	x1, x2;

public:
	Vector2(Scalar _x1=0.0, Scalar _x2=0.0) {x1=_x1; x2=_x2;};
	Vector2(const Vector2& v) {x1 = v.x1; x2 = v.x2;};
	Scalar e1() const {return x1;};
	Scalar e2() const  {return x2;};
	Scalar  e(int component){
		switch (component)
			{
			case 1:
				return x1;
			case 2:
				return x2;
			}
		return (0.0);
	}
	void	set_e1(Scalar _x1) {x1 = _x1;};
	void	set_e2(Scalar _x2) {x2 = _x2;};
	void set(int component, Scalar value){
		switch (component)
			{
			case 1:
				x1 = value;
				break;
			case 2:
				x2 = value;
				break;
			}
	}
	Vector2 operator + (const Vector2& v) {return Vector2(x1+v.x1, x2+v.x2);};
	Vector2	operator - (const Vector2& v) {return Vector2(x1-v.x1, x2-v.x2);};
	Vector2	operator * (Scalar c) {return Vector2(c*x1, c*x2);};
	Vector2	operator / (Scalar c) {Scalar ic = 1/c; return Vector2(ic*x1,
		ic*x2);};
	Scalar	operator * (const Vector2& v) {return x1*v.x1 + x2*v.x2;};
	Vector2&	operator += (Scalar c) {x1 += c; x2 += c; return *this;};
	Vector2&	operator -= (Scalar c) {x1 -= c; x2 -= c; return *this;};
	Vector2&	operator *= (Scalar c) {x1 *= c; x2 *= c; return *this;};
	Vector2&	operator *= (const Vector2& v) {x1 *= v.x1; x2 *= v.x2;
		return *this;};
	Vector2&	operator += (const Vector2& v) {x1 += v.x1; x2 += v.x2; return *this;};
	Vector2&	operator -= (const Vector2& v) {x1 -= v.x1; x2 -= v.x2; return *this;};
	Vector2	jvMult(const Vector2& v) {return Vector2(x1*v.x1, x2*v.x2);};
	Vector2	jvDivide(const Vector2& v) {return Vector2(x1/v.x1, x2/v.x2);}; // added by kc 8-30-94
	Scalar	magnitude() {return sqrt(x1*x1 + x2*x2);};
	Vector2	ComponentSqrt(){return Vector2(sqrt(x1), sqrt(x2));}
	Vector2 ComponentSqr(){return Vector2(sqr(x1), sqr(x2));}
	int	isNonZero() {if (x1 == 0.0 && x2 == 0.0) return 0;
		else return 1;}
};

inline Vector2	operator *(Scalar c, const Vector2& v)
{
	return Vector2(c*v.e1(), c*v.e2());
}

#endif	//	ifndef __OVECTOR_H
