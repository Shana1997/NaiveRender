#ifndef __MATH_H__
#define __MATH_H__
#include<iostream>
#define PI 3.14159265358
/*class vec2*/
class vec2 {
public:
	vec2();
	vec2(float x, float y);

	float x() const;
	float y() const;
	float& operator[](int i);
	float operator[](int i) const;

	vec2 operator-() const;
	vec2& operator*=(const float t);
	vec2& operator/=(const float t);

	float norm() const;
public:
	float vec[2];
};
std::ostream& operator<<(std::ostream& cout, const vec2 v);
vec2 operator+(const vec2& left, const vec2& right);
vec2 operator-(const vec2& left, const vec2& right);
vec2 operator*(const vec2& left, float t);
vec2 operator/(const vec2 left, float t);
/*class vec3*/
class vec3
{
public:
	vec3();
	vec3(float x, float y, float z);

	float x() const;
	float y() const;
	float z() const;
	float& operator[](int i);
	float operator[](int i)const;

	vec3 operator-()const;
	vec3& operator+=(const vec3& v);
	vec3& operator*=(const float t);
	vec3& operator/=(const float t);

	float norm()const;
public:
	float vec[3];
};
std::ostream& operator<<(std::ostream& out, const vec3& v);
vec3 operator+(const vec3& u, const vec3& v);
vec3 operator-(const vec3& u, const vec3& v);
vec3 operator*(const vec3& u, const vec3& v);
vec3 operator*(double t, const vec3& v);
vec3 operator*(const vec3& v, double t);
vec3 operator/(vec3 v, double t);
double dot(const vec3& u, const vec3& v);
vec3 cross(const vec3& u, const vec3& v);
vec3 normalnize(const vec3& v);
vec3 cwise_product(const vec3& a, const vec3& b);
/*class vec4*/
class vec4
{
public:
	vec4();
	vec4(float x, float y, float z, float w);

	float x()const;
	float y()const;
	float z()const;
	float w()const;

	float& operator[](int i);
	float operator[](int i)const;

	vec4& operator*=(const float t);
	vec4& operator/=(const float t);

public:
	float vec[4];
};
std::ostream& operator<<(std::ostream& out, const vec4& v);
vec4 to_vec4(const vec3& u, float w);
vec4 operator-(const vec4& u, const vec4& v);
vec4 operator+(const vec4& u, const vec4& v);
vec4 operator*(double t, const vec4& v);
vec4 operator*(const vec4& v, double t);

//mat3
class mat3
{
public:
	mat3();

	vec3& operator[](int i);
	vec3 operator[](int i)const;

	mat3 transpose()const;
	mat3 inverse()const;
	mat3 inverse_transpose()const;
	static mat3 identity();

public:
	vec3 rows[3];

};
std::ostream& operator<<(std::ostream& out, const mat3& m);

//mat4
class mat4
{
public:
	mat4();

	vec4& operator[](int i);
	vec4 operator[](int i)const;

	mat4 transpose()const;
	mat4 inverse()const;
	mat4 inverse_transpose()const;
	static mat4 identity();

public:
	vec4 rows[4];
};
std::ostream& operator<<(std::ostream& out, const mat4& m);
vec4 operator*(const mat4& m, const vec4 v);
mat4 operator*(const mat4& m1, const mat4& m2);

//transformation related functions
mat4 translate(float tx, float ty, float tz);
mat4 scale(float sx, float sy, float sz);
mat4 rotate_x(float angle);
mat4 rotate_y(float angle);
mat4 rotate_z(float angle);
mat4 lookAt(vec3 eyePos, vec3 target, vec3 up);
mat4 ortho(float left, float right, float bottom, float top, float near, float far);
mat4 perspective(float fov, float aspect, float near, float far);

//other related functions
static std::tuple<float, float, float> getBarycentric(const vec3 A, const vec3 B, const vec3 C, const vec3 P)
{
	float c = ((A.y() - B.y()) * P.x() + (B.x() - A.x()) * P.y() + A.x() * B.y() - B.x() * A.y()) /
		((A.y() - B.y()) * C.x() + (B.x() - A.x()) * C.y() + A.x() * B.y() - B.x() * A.y());
	float b = ((A.y() - C.y()) * P.x() + (C.x() - A.x()) * P.y() + A.x() * C.y() - C.x() * A.y()) /
		((A.y() - C.y()) * B.x() + (C.x() - A.x()) * B.y() + A.x() * C.y() - C.x() * A.y());
	float a = 1 - b - c;
	return { a,b,c };
}//static必须在文件范围内定义函数
#endif // !__MATH_H__
