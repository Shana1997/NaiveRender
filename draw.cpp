#include"draw.h"
void Line(vec3 p0, vec3 p1, TGAImage& image, TGAColor color)
{
	bool steep = false;
	if (fabs(p0.x() - p1.x()) < fabs(p0.y() - p1.y()))
	{
		std::swap(p0[0], p0[1]);
		std::swap(p0[0], p1[1]);
		steep = true;
	}
	if (p0.x() > p1.x())
	{
		std::swap(p0[0], p1[0]);
		std::swap(p0[1], p1[1]);
	}
	float k = (p1.y() - p0.y()) / (p1.x() - p0.x());
	float b = p0.y() - k * p0.x();
	for (float t = p0.x(); t <= p1.x(); t++)
	{
		float x = t;
		float y = k * x + b;
		if (steep == false)
		{
			image.set(x, y, color);
		}
		else
		{
			image.set(y, x, color);
		}
	}
}