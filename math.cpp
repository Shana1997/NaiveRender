#include"math.h"

vec2::vec2()
{
	this->vec[0] = 0.0f;
	this->vec[1] = 0.0f;
}
vec2::vec2(float x, float y)
{
	this->vec[0] = x;
	this->vec[1] = y;
}
float vec2::x() const
{
	return this->vec[0];
}
float vec2::y() const //如果没有const，重载+号中的v.y()会报错，x()同理
{
	return this->vec[1];
}

vec2 vec2::operator-() const
{
	return vec2(-vec[0], -vec[1]);
}
vec2& vec2::operator*=(const float t)//此处加const下面的this会变成非左值
{
	vec[0] *= t;
	vec[1] *= t;
	return *this;
}
vec2& vec2::operator/=(const float t)
{
	vec[0] /= t;
	vec[1] /= t;
	return *this;
}
float vec2::norm() const
{
	return std::sqrtf(vec[0] * vec[0] + vec[1] * vec[1]);
}
std::ostream& operator<<(std::ostream& cout, const vec2 v)
{
	std::cout << v.x() << "  " << v.y();
	return cout;
}
vec2 operator+(const vec2& left, const vec2& right)
{
	return vec2(left.x() + right.x(), left.y() + right.y());
}
vec2 operator-(const vec2& left, const vec2& right)
{
	return vec2(left.x() - right.x(), left.y() - right.y());
}
float& vec2::operator[](int i)
{
	return vec[i];
}
float vec2::operator[](int i) const
{
	return vec[i];
}
vec2 operator*(const vec2& left, float t)
{
	return vec2(left.x() * t, left.y() * t);
}
vec2 operator/(const vec2 left, float t)
{
	return vec2(left.x() / t, left.y() / t);
}
//vec3
vec3::vec3() :vec{ 0,0,0 } {}
vec3::vec3(float x, float y, float z) :vec{ x,y,z } {}

float vec3::x() const
{
	return vec[0];
}
float vec3::y() const
{
	return vec[1];
}
float vec3::z() const
{
	return vec[2];
}
float& vec3::operator[](int i)
{
	return vec[i];
}
float vec3::operator[](int i)const
{
	return vec[i];
}

vec3 vec3::operator-()const
{
	return vec3(-vec[0], -vec[1], -vec[2]);
}
vec3& vec3::operator+=(const vec3& v)
{
	//return vec3(vec[0] + v[0], vec[1] + v[1], vec[2] + v[2]);
	vec[0] += v[0]; vec[1] += v[1]; vec[2] += v[2];
	return *this;
}
vec3& vec3::operator*=(const float t)
{
	vec[0] *= t; vec[1] *= t; vec[2] *= t;
	return *this;
}
vec3& vec3::operator/=(const float t)
{
	vec[0] /= t; vec[1] /= t; vec[2] /= t;
	return *this;
}

float vec3::norm()const
{
	return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}
std::ostream& operator<<(std::ostream& out, const vec3& v)
{
	out << v.x() << ", " << v.y() << ", " << v.z();
	return out;
}
vec3 operator+(const vec3& u, const vec3& v)
{
	return vec3(u.x() + v.x(), u.y() + v.y(), u.z() + v.z());
}
vec3 operator-(const vec3& u, const vec3& v)
{
	return vec3(u.x() - v.x(), u.y() - v.y(), u.z() - v.z());
}
vec3 operator*(const vec3& u, const vec3& v)
{
	return vec3(u.x() * v.x(), u.y() * v.y(), u.z() * v.z());
}
vec3 operator*(double t, const vec3& v)
{
	return vec3(t * v.x(), t * v.y(), t * v.z());
}
vec3 operator*(const vec3& v, double t)
{
	return t * v;
}
vec3 operator/(vec3 v, double t)
{
	return (1 / t) * v;
}
double dot(const vec3& u, const vec3& v)
{
	return u.x() * v.x() + u.y() * v.y() + u.z() * v.z();
}
vec3 cross(const vec3& u, const vec3& v)
{
	return vec3(u.y() * v.z() - u.z() * v.y(),
		u.z() * v.x() - u.x() * v.z(),
		u.x() * v.y() - u.y() * v.x());
}
vec3 normalnize(const vec3& v)
{
	return v / v.norm();
}
vec3 cwise_product(const vec3& a, const vec3& b)
{
	return a * b;
}
//vec4
vec4::vec4() :vec{0,0,0,0 } {}
vec4::vec4(float x, float y, float z, float w) : vec{x, y, z, w} {}

float vec4::x()const
{
	return vec[0];
}
float vec4::y()const
{
	return vec[1];
}
float vec4::z()const
{
	return vec[2];
}
float vec4::w()const
{
	return vec[3];
}

float& vec4::operator[](int i)
{
	return vec[i];
}
float vec4::operator[](int i)const
{
	return vec[i];
}

vec4& vec4::operator*=(const float t)
{
	vec[0] *= t;
	vec[1] *= t;
	vec[2] *= t;
	vec[3] *= t;
	return *this;
}
vec4& vec4::operator/=(const float t)
{
	vec[0] /= t;
	vec[1] /= t;
	vec[2] /= t;
	vec[3] /= t;
	return *this;
}

std::ostream& operator<<(std::ostream& out, const vec4& v)
{
	out << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3];
	return out;
}
vec4 to_vec4(const vec3& u, float w)
{
	return vec4(u[0], u[1], u[2], w);
}
vec4 operator-(const vec4& u, const vec4& v)
{
	return vec4(u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3]);
}
vec4 operator+(const vec4& u, const vec4& v)
{
	return vec4(u[0] + v[0], u[1] + v[1], u[2] + v[2], u[3] + v[3]);
}
vec4 operator*(double t, const vec4& v)
{
	return vec4(t * v[0], t * v[1], t * v[2], t * v[3]);
}
vec4 operator*(const vec4& v, double t)
{
	return vec4(t * v[0], t * v[1], t * v[2], t * v[3]);
}

//mat3
mat3::mat3() {}
vec3& mat3::operator[](int i)
{
	return rows[i];
}
vec3 mat3::operator[](int i)const
{
	return rows[i];
}
/*for matrix's determinant adjoint and inverse*/
static float mat3_determinant(const mat3& m)
{
	float a = +m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]);
	float b = -m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]);
	float c = +m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
	return a + b + c;
}

static mat3 mat3_adjoint(const mat3& m)
{
	mat3 adjoint;
	adjoint[0][0] = +(m[1][1] * m[2][2] - m[2][1] * m[1][2]);
	adjoint[0][1] = -(m[1][0] * m[2][2] - m[2][0] * m[1][2]);
	adjoint[0][2] = +(m[1][0] * m[2][1] - m[2][0] * m[1][1]);
	adjoint[1][0] = -(m[0][1] * m[2][2] - m[2][1] * m[0][2]);
	adjoint[1][1] = +(m[0][0] * m[2][2] - m[2][0] * m[0][2]);
	adjoint[1][2] = -(m[0][0] * m[2][1] - m[2][0] * m[0][1]);
	adjoint[2][0] = +(m[0][1] * m[1][2] - m[1][1] * m[0][2]);
	adjoint[2][1] = -(m[0][0] * m[1][2] - m[1][0] * m[0][2]);
	adjoint[2][2] = +(m[0][0] * m[1][1] - m[1][0] * m[0][1]);
	return adjoint;
}
mat3 mat3::transpose()const
{
	mat3 temp;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			temp[i][j] = rows[j][i];
		}
	}
	return temp;
}
mat3 mat3::inverse_transpose()const
{
	int i, j;
	mat3 inverse_transpose, adjoint;
	float determinant, inv_determinant;

	adjoint = mat3_adjoint(*this);
	determinant = mat3_determinant(*this);
	inv_determinant = 1 / determinant;

	for (i = 0; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			inverse_transpose[i][j] = adjoint[i][j] * inv_determinant;

	return inverse_transpose;
}
mat3 mat3::inverse()const
{
	return ((*this).inverse_transpose()).transpose();
}
mat3 mat3::identity()
{
	mat3 identity;
	identity[0][0] = 1;
	identity[1][1] = 1;
	identity[2][2] = 1;
	return identity;
}
std::ostream& operator<<(std::ostream& out, const mat3& m)
{
	out << m[0] << std::endl << m[1] << std::endl << m[2];
	return out;
}

//mat4
mat4::mat4() : rows{vec4(), vec4(),vec4(),vec4()} {}

vec4& mat4::operator[](int i)
{
	return rows[i];
}
vec4 mat4::operator[](int i)const
{
	return rows[i];
}
vec4 operator*(const mat4& m, const vec4 v)
{
	vec4 result;
	for (int i = 0; i < 4; i++)
	{
		float a = m[i][0] * v[0];
		float b = m[i][1] * v[1];
		float c = m[i][2] * v[2];
		float d = m[i][3] * v[3];
		result[i] = a + b + c + d;
	}
	return result;
}
mat4 operator*(const mat4& m1, const mat4& m2)
{
	mat4 result;
	int i, j, k;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++)
				result[i][j] += m1[i][k] * m2[k][j];
	return result;
}

mat4 mat4::transpose()const
{
	mat4 temp;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			temp[i][j] = rows[j][i];
		}
	}
	return temp;
}

/*for matrix's minor, cofactor, and adjoint*/
static float mat4_minor(mat4 m, int r, int c)
{
	mat3 cut_down;
	int i, j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			int row = i < r ? i : i + 1;
			int col = j < c ? j : j + 1;
			cut_down[i][j] = m[row][col];
		}
	}
	return mat3_determinant(cut_down);
}

static float mat4_cofactor(mat4 m, int r, int c)
{
	float sign = (r + c) % 2 == 0 ? 1.0f : -1.0f;
	float minor = mat4_minor(m, r, c);
	return sign * minor;
}

static mat4 mat4_adjoint(mat4 m)
{
	mat4 adjoint;
	int i, j;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			adjoint[i][j] = mat4_cofactor(m, i, j);
	return adjoint;
}
mat4 mat4::inverse_transpose()const
{
	int i, j;
	float determinant, inv_determinant;
	mat4 adjoint, inverse_transpose;

	adjoint = mat4_adjoint(*this);
	determinant = 0;
	for (i = 0; i < 4; i++)
	{
		determinant += rows[0][i] * adjoint[0][i];
	}
	inv_determinant = 1 / determinant;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			inverse_transpose[i][j] = adjoint[i][j] * inv_determinant;

	return inverse_transpose;
}
mat4 mat4::inverse()const
{
	return ((*this).inverse_transpose()).transpose();
	//return this->inverse_transpose().transpose();
}
mat4 mat4::identity()
{
	mat4 identity;
	for (int i = 0; i < 4; i++) identity[i][i] = 1;
	return identity;
}
std::ostream& operator<<(std::ostream& out, const mat4& m)
{
	out << m[0] << std::endl << m[1] << std::endl << m[2] << std::endl << m[3];
	return out;
}
//transformation related functions
mat4 translate(float tx, float ty, float tz)
{
	mat4 trans = mat4::identity();
	trans[0][3] = tx;
	trans[1][3] = ty;
	trans[2][3] = tz;
	return trans;
}
mat4 scale(float sx, float sy, float sz)
{
	mat4 scale = mat4::identity();
	scale[0][0] = sx;
	scale[1][1] = sy;
	scale[2][2] = sz;
	return scale;
}
mat4 rotate_x(float angle)
{
	angle = angle * PI / 180.f;
	float c = cos(angle);
	float s = sin(angle);
	mat4 rotateX = mat4::identity();
	rotateX[1][1] = c;
	rotateX[1][2] = -s;
	rotateX[2][1] = s;
	rotateX[2][2] = c;
	return rotateX;
}
mat4 rotate_y(float angle)
{
	angle = angle * PI / 180.f;
	float c = cos(angle);
	float s = sin(angle);
	mat4 rotateY = mat4::identity();
	rotateY[0][0] = c;
	rotateY[0][2] = s;
	rotateY[2][0] = -s;
	rotateY[2][2] = c;
	return rotateY;
}
mat4 rotate_z(float angle)
{
	angle = angle * PI / 180.f;
	float c = cos(angle);
	float s = sin(angle);
	mat4 rotateZ = mat4::identity();
	rotateZ[0][0] = c;
	rotateZ[0][1] = -s;
	rotateZ[1][0] = s;
	rotateZ[1][1] = c;
	return rotateZ;
}
mat4 lookAt(vec3 eyePos, vec3 target, vec3 up)
{
	mat4 lookat = mat4::identity();
	vec3 z = normalnize(eyePos - target);
	vec3 x = normalnize(cross(up, z));
	vec3 y = normalnize(cross(z, x));

	lookat[0][0] = x[0];
	lookat[0][1] = x[1];
	lookat[0][2] = x[2];

	lookat[1][0] = y[0];
	lookat[1][1] = y[1];
	lookat[1][2] = y[2];

	lookat[2][0] = z[0];
	lookat[2][1] = z[1];
	lookat[2][2] = z[2];

	lookat[0][3] = -dot(x, eyePos);
	lookat[1][3] = -dot(y, eyePos);
	lookat[2][3] = -dot(z, eyePos);
	return lookat;
}
mat4 ortho(float left, float right, float bottom, float top, float near, float far)
{
	float x_range = right - left;
	float y_range = top - bottom;
	float z_range = near - far;//看向-z
	mat4 ortho = mat4::identity();
	ortho[0][0] = 2 / x_range;
	ortho[1][1] = 2 / y_range;
	ortho[2][2] = 2 / z_range;
	ortho[0][3] = -(left + right) / x_range;
	ortho[1][3] = -(bottom + top) / y_range;
	ortho[2][3] = -(near + far) / z_range;
	return ortho;
}
mat4 perspective(float fov, float aspect, float near, float far)
{
	fov = fov * PI / 180.f;
	float t = fabs(near) * tan(fov / 2);
	float r = aspect * t;

	mat4 perspect = mat4::identity();
	perspect[0][0] = near / r;
	perspect[1][1] = near / t;
	perspect[2][2] = (near + far) / (near - far);
	perspect[2][3] = 2 * near * far / (far - near);
	perspect[3][2] = 1;
	perspect[3][3] = 0;
	return perspect;
}
//other related functions
/*std::tuple<float, float, float> getBarycentric(const vec3 A, const vec3 B, const vec3 C, const vec3 P)
{
	float c = ((A.y() - B.y()) * P.x() + (B.x() - A.x()) * P.y() + A.x() * B.y() - B.x() * A.y()) /
		((A.y() - B.y()) * C.x() + (B.x() - A.x()) * C.y() + A.x() * B.y() - B.x() * A.y());
	float b = ((A.y() - C.y()) * P.x() + (C.x() - A.x()) * P.y() + A.x() * C.y() - C.x() * A.y()) /
		((A.y() - C.y()) * B.x() + (C.x() - A.x()) * B.y() + A.x() * C.y() - C.x() * A.y());
	float a = 1 - b - c;
	return { a,b,c };
}*/