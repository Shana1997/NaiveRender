#ifndef __CAMERA_H__
#define __CAMERA_H__
#include"math.h"
class Camera
{
public:
	Camera(vec3 pos, vec3 target, vec3 up, float aspect);
	~Camera();

	vec3 position;
	vec3 target;
	vec3 up;
	vec3 x;
	vec3 y;
	vec3 z;
	float aspect;
};

#endif // !__CAMERA_H__
