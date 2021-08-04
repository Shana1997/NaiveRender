#ifndef __PIPELINE_H__
#define __PIPELINE_H__
#include"model.h"
#include"TGAimage.h"
#include"texture.h"
#include"Light.hpp"
struct Payload
{
public:
	Payload() = default;
	//下面的函数如果再结构体外实现编译会通不过
	Payload(mat4 m, mat4 v, mat4 p, Model model, int index)
	{
		this->model = m;
		view = v;
		projection = p;
		for (int i = 0; i < 3; i++)
		{
			pos_modelSpace[i] = model.vert(model.fs(index).vert[i] - 1);
			normals[i] = model.getNormal(model.fs(index).normal[i] - 1);
			uv[i] = model.getTexcoord(model.fs(index).tex[i] - 1);
		}
	}
	mat4 model;
	mat4 view;
	mat4 projection;

	//Camera camera;
	//Light light;

	vec3 normals[3];
	vec3 pos_modelSpace[3];
	vec4 pos_worldSpace[3];
	vec4 pos_viewSpace[3];
	vec4 pos_clipSpace[3];
	vec3 pos_screenSpace[3];
	vec2 uv[3];

};

class pipeline
{
public:
	pipeline() = default;
	pipeline(int w, int h, Payload & payload);
	vec3 clip2screen(vec4 clipPos)const;
	void vertex_shader();
	void rasterize_triangle(TGAImage& image, float* zBuffer, texture& texImage, Light& light);
	vec3 fragment_shader(texture& texImage,Light& light, vec3& shadingP, vec3& Normal, vec2& uv);

private:
	int width;
	int height;
	Payload payload;
};
#endif // !__PIPELINE_H__
