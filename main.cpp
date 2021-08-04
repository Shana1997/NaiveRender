#include<iostream>
#include"model.h"
#include"draw.h"
#include"Light.hpp"
#include"camera.h"
#include"texture.h"
#include"pipeline.h"
#include<time.h>
#define WIDTH 1920
#define HEIGHT 1920
//模型加载->纹理加载->vertex shader->fragment shader->输出
Camera camera(vec3(0, 5, 10), vec3(0, 0, 0), vec3(0, 1, 0), 1.f);
Light light(vec3(0, 0, 1), vec3(1.0, 1.0, 1.0));
TGAImage image(WIDTH, HEIGHT, TGAImage::RGB);
int main(int argc, char** argv)
{	
	clock_t start, end;
	mat4 model = mat4::identity();
	mat4 view = mat4::identity();
	mat4 projection = mat4::identity();
	view = lookAt(camera.position, camera.target, camera.up);
	projection = perspective(50.f, 1.f, 0.1f, 50.f);
	//Model faceModel("../models/african_head.obj");
	//texture tex("../textures/african_head_diffuse.tga");
	Model faceModel("../models/mary.obj");
	texture tex("../textures/MC003_Kozakura_Mari.tga");
	float* zBuffer = new float[WIDTH * HEIGHT];
	for (int i = 0; i < WIDTH * HEIGHT; i++)zBuffer[i] = -std::numeric_limits<float>::max();
	start = clock();
	/*render one frame*/
	for (int i = 0; i < faceModel.facesNumber; i++)
	{
		Payload payload(model, view, projection, faceModel, i);
		pipeline renderLine(WIDTH, HEIGHT, payload);
		renderLine.vertex_shader();
		
		renderLine.rasterize_triangle(image, zBuffer, tex, light);
		//renderLine.fragment_shader(image);
	}
	image.write_tga_file("test.tga");
	end = clock();
	std::cout << double(end - start) / CLOCKS_PER_SEC;
	return 0;
}