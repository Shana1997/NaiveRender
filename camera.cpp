#include"camera.h"
Camera::Camera(vec3 pos, vec3 target, vec3 up, float aspect)
	:position(pos),target(target),up(up),aspect(aspect) {}
Camera::~Camera(){}

