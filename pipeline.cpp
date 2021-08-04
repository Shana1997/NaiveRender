#include"pipeline.h"
pipeline::pipeline(int w, int h, Payload& payload)
{
	width = w; height = h;
	this->payload = payload;
}
vec3 pipeline::clip2screen(vec4 clipPos)const
{
	clipPos /= clipPos.w();//to NDC space
	clipPos[0] = 0.5f * width * (clipPos.x() + 1.0f);
	clipPos[1] = 0.5f * height * (clipPos.y() + 1.0f);
	return vec3(int(clipPos.x()), int(clipPos.y()), 1.0f);
}
void pipeline::vertex_shader()
{
	for (int i = 0; i < 3; i++)
	{
		//to world space
		payload.pos_worldSpace[i] =
			payload.model * to_vec4(payload.pos_modelSpace[i], 1.f);
		//std::cout << payload.pos_modelSpace[i] << std::endl;
		//to view space
		payload.pos_viewSpace[i] =
			payload.view * payload.model * to_vec4(payload.pos_modelSpace[i], 1.f);
		//to clip space
		payload.pos_clipSpace[i] =
			payload.projection * payload.view * payload.model *
			to_vec4(payload.pos_modelSpace[i], 1.f);
		//to screen space
		payload.pos_screenSpace[i] =
			clip2screen(payload.pos_clipSpace[i]);
	}
}
vec3 pipeline::fragment_shader(texture& texImage, Light& light, vec3& shadingP, vec3& Normal, vec2& uv)
{
	//将法线 灯位置转换到视角空间
	vec4 normal_view = payload.view * payload.model * to_vec4(Normal, 1.f);
	vec4 lightPos_view = payload.view * payload.model * to_vec4(light.position, 1.f);
	vec3 normal = normalnize(vec3(normal_view.x(), normal_view.y(), normal_view.z()));
	vec3 lightPos = vec3(lightPos_view.x(), lightPos_view.y(), lightPos_view.z());
	vec3 lightDir = normalnize(lightPos - shadingP);
	vec3 viewDir = normalnize(-shadingP);//视角空间的摄像机位置为原点

	TGAColor texColor = texImage.getTexColor(uv.x(), uv.y());
	vec3 albedo = vec3(texColor.r, texColor.g, texColor.b) / 255.0f;
	//get ambient
	vec3 ambient = vec3(0.05, 0.05, 0.05);
	//get diffuse
	vec3 diffuse = albedo * light.intensity * max(0.0, dot(normal, lightDir));
	//get specular
	vec3 halfVector = normalnize(lightDir + viewDir);
	vec3 specularColor = vec3(1.0f, 1.0f, 1.0f);
	vec3 specular = specularColor * light.intensity * pow(max(0.0, dot(normal, halfVector)), 50.0);
	
	vec3 color = (ambient + diffuse + specular) * 255;
	//return normalnize(Normal) * 255.f;
	/*test*/
	vec3 ab = payload.pos_modelSpace[1] - payload.pos_modelSpace[0];
	vec3 ac = payload.pos_modelSpace[2] - payload.pos_modelSpace[0];
	vec3 n = normalnize(cross(ab, ac));
	vec3 lightDirection = vec3(0, 0, 1);
	float intensity = max((dot(n, lightDirection) / (n.norm() * lightDirection.norm())), 0.0);
	return intensity * albedo * 255.f;
	/*test end*/
}
void pipeline::rasterize_triangle(TGAImage& image, float* zBuffer, texture& texImage, Light& light)
{
	int xMin, xMax, yMin, yMax;
	xMin = std::min(payload.pos_screenSpace[0].x(),
		std::min(payload.pos_screenSpace[1].x(), payload.pos_screenSpace[2].x()));
	xMin = std::max(0, xMin);
	xMax = std::max(payload.pos_screenSpace[0].x(),
		std::max(payload.pos_screenSpace[1].x(), payload.pos_screenSpace[2].x()));
	xMax = std::min(width - 1, xMax);
	yMin = std::min(payload.pos_screenSpace[0].y(),
		std::min(payload.pos_screenSpace[1].y(), payload.pos_screenSpace[2].y()));
	yMin = std::max(0, yMin);
	yMax = std::max(payload.pos_screenSpace[0].y(),
		std::max(payload.pos_screenSpace[1].y(), payload.pos_screenSpace[2].y()));
	yMax = std::min(height - 1, yMax);
	vec3 p;
	for (p[0] = xMin; p[0] <= xMax; p[0]++)
	{
		for (p[1] = yMin; p[1] <= yMax; p[1]++)
		{
			auto [alpha, beta, gamma] = getBarycentric(
				payload.pos_screenSpace[0] + vec3(0.5, 0.5, 0.5), payload.pos_screenSpace[1] + vec3(0.5, 0.5, 0.5), payload.pos_screenSpace[2] + vec3(0.5, 0.5, 0.5), p
			);
			if (alpha < 0 || beta < 0 || gamma < 0)continue;
			p[2] = 0;
			float Zt = 1.0 /
				(alpha / payload.pos_viewSpace[0].z() + beta / payload.pos_viewSpace[1].z() + gamma / payload.pos_viewSpace[2].z());
			//深度插值
			p[2] +=
				(payload.pos_viewSpace[0].z() * alpha / payload.pos_viewSpace[0].z() +
					payload.pos_viewSpace[1].z() * beta / payload.pos_viewSpace[1].z() +
					payload.pos_viewSpace[2].z() * gamma / payload.pos_viewSpace[2].z()) / (1.0 / Zt);
			//法线插值
			/*vec3 Normal =
				(payload.normals[0] * alpha / payload.pos_viewSpace[0].z() +
					payload.normals[1] * beta / payload.pos_viewSpace[1].z() +
					payload.normals[2] * gamma / payload.pos_viewSpace[2].z()) / (1.0 / Zt);*/
			vec3 Normal = payload.normals[0] * alpha + payload.normals[1] * beta * payload.normals[2] * gamma;
			//纹理坐标插值
			/*vec2 UV =
				(payload.uv[0] * alpha / payload.pos_viewSpace[0].z() +
					payload.uv[1] * beta / payload.pos_viewSpace[1].z() +
					payload.uv[2] * gamma / payload.pos_viewSpace[2].z()) / (1.0 / Zt);*/
			vec2 UV = payload.uv[0] * alpha + payload.uv[1] * beta + payload.uv[2] * gamma;
			//着色点插值
			vec3 pos_view0 = vec3(payload.pos_viewSpace[0].x(), payload.pos_viewSpace[0].y(), payload.pos_viewSpace[0].z());
			vec3 pos_view1 = vec3(payload.pos_viewSpace[1].x(), payload.pos_viewSpace[1].y(), payload.pos_viewSpace[1].z());
			vec3 pos_view2 = vec3(payload.pos_viewSpace[2].x(), payload.pos_viewSpace[2].y(), payload.pos_viewSpace[2].z());
			vec3 shadingP =
				(pos_view0 * alpha + pos_view1 * beta + pos_view2 * gamma) / (1.0 / Zt);
			int idx = p.y() * image.get_width() + p.x();
			if (zBuffer[idx] < p.z())
			{
				zBuffer[idx] = p.z();
				vec3 phongColor = fragment_shader(texImage, light, shadingP, Normal, UV);
				TGAColor color = TGAColor(phongColor.x(), phongColor.y(), phongColor.z(), 255);
				image.set(p.x(), p.y(), color);
			}
		}
	}

}