#ifndef __LIGHT_H__
#define __LIGHT_H__
#include"math.h"
struct Light {
public:
	Light() = default;
	Light(vec3 pos, vec3 intensity) : position(pos), intensity(intensity) {}
	vec3 position;
	vec3 intensity;
};
#endif // !__LIGHT_H__
