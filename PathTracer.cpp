
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>

#include "TGA.h"
#include "Framebuffer.h"


struct float3
{
	float x, y, z;
	inline float3() {}
	inline float3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
	inline float operator[](int i) const { return (&x)[i]; }
	inline float3 operator-() const { return float3(-x, -y, -z); }
	inline float3 operator*(const float rhs) const { return float3(x * rhs, y * rhs, z * rhs); }
	inline float3& operator*=(const float rhs) { x *= rhs; y *= rhs; z *= rhs; return *this; }
	inline float3 operator+(const float3& rhs) const { return float3(x + rhs.x, y + rhs.y, z + rhs.z); }
	inline float3 operator-(const float3& rhs) const { return float3(x - rhs.x, y - rhs.y, z - rhs.z); }
	inline float3& operator+=(const float3& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
	inline float sqrLength() const { return x * x + y * y + z * z; }
	inline float length() const { return sqrtf(x * x + y * y + z * z); }
	inline float3& normalize() { operator*(1.f / length()); return *this; }
	inline float dot(const float3& rhs) const { return x * rhs.x + y * rhs.y + z * rhs.z; }
	inline float3 cross(const float3& rhs) const { return float3(y*rhs.z - z*rhs.y, z*rhs.x - x*rhs.z, x*rhs.y - y*rhs.x); }	
};
inline float3 operator*(const float3& lhs, const float3& rhs) { return float3(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z); }
typedef float3 vec3;

struct Ray
{
	static const float epsilon;
	float3 origin;
	float3 direction;
	float mint, maxt;
	inline float3 operator()(float t) const { return origin + direction*t; }
};
const float Ray::epsilon = 1e-3f;

struct HitRecord
{
	vec3 point; float t;
	vec3 normal; int material;
};

struct Sphere
{
	float3 center; 
	float radius;

	// in world space
	bool hit(const Ray& ray, HitRecord &rec) const
	{
		// to is in object space
		float3 to = ray.origin - center;
		float A = ray.direction.dot(ray.direction);
		float B = 2.f * to.dot(ray.direction);
		float C = to.dot(to) - radius*radius;
#if 0
		// pbrt
		float discriminant = B*B - 4.f*A*C;
		if (discriminant > 0.f) {
			float rootDiscriminant = sqrt(discriminant);
			float q = B < 0.f ? -0.5f * (B - rootDiscriminant) : -0.5f * (B + rootDiscriminant);
			float t0 = q / A;
			float t1 = C / q;
			if (t0 > t1) { float tmp = t1; t1 = t0; t0 = tmp; }
			if (t0 > ray.maxt || t1 < ray.mint)
				return false;
			float thit = t0;
			if (t0 < ray.mint) {
				thit = t1;
				if (thit > ray.maxt)
					return false;
			}
			return true;
			/*if (t0 <= ray.maxt && t1 >= ray.mint) {
				float thit = (t0 < ray.mint) ? t1 : t0;
				p = ray(thit);
				return thit < ray.maxt;
			}*/

		}
#else
		// glassner:
		B /= 2.0f;
		float discriminant = B*B - A*C;
		if (discriminant > 0.f) {
			float rootDiscriminant = sqrtf(discriminant);
			float temp = (-B - rootDiscriminant) / A;
			if (temp < ray.maxt && temp > ray.mint) { rec.point = ray(temp); rec.normal = (rec.point - center) * (1.f/radius); rec.t = temp; return true; }
			temp = (-B + rootDiscriminant) / A;
			if (temp < ray.maxt && temp > ray.mint) { rec.point = ray(temp); rec.normal = (rec.point - center) * (1.f/radius); rec.t = temp; return true; }
		}
#endif
		return false;
	}
};

struct Disk
{
	float height;
	float radius;

	bool hit(const Ray& ray, float3& point) const
	{
		// to is in object space
		float3 to = ray.origin - float3(0.f, 0.f, height);
		float thit = (height - to.y) / (to.y + 0.00001f);
		if (thit < ray.mint || thit > ray.maxt)
			return false;
		point = ray(thit);
		float distSqr = point.x*point.x + point.z*point.z;
		if (distSqr > radius*radius)
			return false;
		return true;
	}
};

struct Plane
{
	union {
		struct { float a, b, c, d; };
		struct { float3 n; float d; };
	};

	bool hit(const Ray& ray, float3& point) const
	{

	}
};

struct Camera
{
	float3 position;
	float3 forward;

	float width, height, nearPlane, farPlane, fov, aspectRatio;
	float invDenom, invTanAngle;
	
	void initialize(float _fov, float _near, float _far, float _width, float _height)
	{
		width = _width; height = _height;
		fov = _fov;  nearPlane = _near; farPlane = _far; aspectRatio = height / width;
		invDenom = 1.f / (farPlane - nearPlane);
		invTanAngle = 1.f / tan(fov * M_PI / 360.0f);
	}
	
	void project(uint16_t& x, uint16_t& y, const float3& p) const
	{
		float xp = p.x*invTanAngle*aspectRatio;
		float yp = p.y*invTanAngle;
		float zp = farPlane*invDenom * p.z - farPlane*nearPlane*invDenom;
		float wp = 1.f / p.z;
		xp *= wp; yp *= wp; zp *= wp;
		// project
		xp /= zp; yp /= zp;
		// convert to window 
		x = (uint16_t)((xp * 0.5f + 0.5f) * width);
		y = (uint16_t)((yp * 0.5f + 0.5f) * height);
	}

	void unproject(float3& p, float x, float y) const
	{
		float xndc = (x / width) * 2.f - 1.f;
		float yndc = (y / height) * 2.f - 1.f;

		p.x = 1.0f * (xndc / (invTanAngle*aspectRatio));
		p.y = 1.0f * (yndc / invTanAngle);
		p.z = nearPlane;
	}

	Ray generateRay(float x, float y) const
	{
		float3 offset;
		unproject(offset, x, y);
		// using same value for position & direction -> perspective ray converging to camera center
		return Ray{ position + offset, (position + offset).normalize(), 0.f, farPlane - nearPlane };
	}
};

float nextFloat(float min, float max)
{
	return min + ((float)rand() / (float)RAND_MAX) * (max - min);
}



struct Material
{
	float r, g, b;
	int type;
	float param; // metal(1)=fuzz, dieletric(2)=ior
};

struct Scene
{
	std::vector<Sphere> spheres;
	std::vector<Material> materials;
	Camera camera;

	vec3 color(Ray& ray, Scene& scene, int depth);
};

vec3 randomDirectionInUnitSphere()
{
	vec3 p;
	do {
		p = vec3(nextFloat(0.0f, 1.f), nextFloat(0.0f, 1.0f), nextFloat(0.0f, 1.0f)) * 2.0f - vec3(1.f, 1.f, 1.f);
	} 
	while (p.sqrLength() >= 1.0f);
	return p;
}

vec3 randomDirectionInHemisphere(const vec3& N)
{
	vec3 p;
	do {
		p = vec3(nextFloat(0.0f, 1.f), nextFloat(0.0f, 1.0f), nextFloat(0.0f, 1.0f)) * 2.0f - vec3(1.f, 1.f, 1.f);
	} while (p.sqrLength() >= 1.0f && N.dot(p) < 0.f);
	return p;
}

bool scatter(const Material& mat, const vec3& position, const vec3& direction, const vec3& N, Ray& scattered, vec3& attenuation)
{
	scattered.origin = position;
	scattered.mint = 0.f; 
	scattered.maxt = FLT_MAX;
	if (mat.type == 0) {
		scattered.direction = (position + N + randomDirectionInHemisphere(N)) - position;
		// Lambert
		//float NdotL = -N.z;
		attenuation = vec3(mat.r, mat.g, mat.b);// * NdotL;
		return true;
	}
	else if (mat.type == 1) {
		vec3 V = direction; 
		//V.normalize();
		float NdotV = V.dot(N);
		vec3 R = V - N*(2.f*NdotV);
		scattered.direction = R;
		// pseudo Phong
		//float VdotR = R.dot(-V);
		//float specular = saturate(pow(VdotR, 8.f));		
		attenuation = vec3(mat.r, mat.g, mat.b);// *specular;
		return N.dot(R) > 0.f;
	}
	else {
		vec3 V = direction;
		//V.normalize();
		float NdotV = V.dot(N);
		vec3 R = V - N*(2.f*NdotV);
		scattered.direction = R;
		// dielectric / refractive
		attenuation = vec3(1.f, 1.f, 1.f);

		// refract
		vec3 outN = N;
		if (NdotV > 0.0f) {
			outN = -N;
		}
		else {
		
		}

		float NdotOV = V.dot(outN);
		float discriminant = 1.0f - mat.param*mat.param*(1.0f - NdotOV*NdotOV);
		vec3 refracted{ 0.f, 0.f, 0.f };

		if (discriminant > 0.0f) {
			refracted = (V - N*NdotV)*mat.param - N*sqrtf(discriminant);
			//scattered.direction = refracted;
		}

		return true;
	}
	
}

vec3 Scene::color(Ray& ray, Scene& scene, int depth)
{
	HitRecord rec;
	int hitIndex = -1;
	int index = 0;

	// get closest hit
	for (auto& sphere : scene.spheres)
	{
		if (sphere.hit(ray, rec)) 
		{
			//ray.maxt = rec.point.z;
			ray.maxt = rec.t;
			hitIndex = index;
		}
		++index;
	}

	if (hitIndex >= 0)
	{
		vec3 N = rec.normal;
		Material& mat = scene.materials[hitIndex];
		vec3 attenuation;
		Ray scattered;
		if (depth < 10 && scatter(mat, rec.point, ray.direction, rec.normal, scattered, attenuation))
		{
			// recurse depending on material
			return attenuation * color(scattered, scene, depth + 1);
		}
		//else // fake emissive
		//	return vec3(1.0f, 1.0f, 1.0f);
	}
	else {
		// ambient sky light
		vec3 dir = ray.direction;
		float t = 0.5f * (dir.y + 1.0f);
		
		return vec3(1.0f, 1.0f, 1.0f)*t;// +vec3(0.0f, 1.0f, 0.0f) * (1.f - t);
	}

	return vec3(0.0f, 0.0f, 0.0f);
}

// ---

int main(void)
{
	// initialize ---

	Framebuffer fb;
	fb.create(1280, 720);

	Scene scene;
	
	scene.camera.position = float3(0.f, 0.f, 0.f);
	scene.camera.initialize(45.f, 0.5f, 100.f, fb.width, fb.height);

	scene.spheres.reserve(10);
	for (int i = 0; i < 10; i++) {
		scene.spheres.emplace_back(Sphere{ { nextFloat(-3.f, +3.f), nextFloat(-3.f, +3.f), nextFloat(0.5f, +10.f) }, nextFloat(0.1f, 1.f) });
	}

	scene.materials.reserve(10);
	for (int i = 0; i < 10; i++) {
		int id = 0;
		float param = 0.0f;
		if (i == 7) {
			id = 2;
			param = 1.5f;
		}
		else
		if (i < 7) {
			param = 0.0f;
		}
		else { 
			id = 1;
			param = 0.0f;
		}
		scene.materials.emplace_back(Material{ nextFloat(0.0f, 1.f), nextFloat(0.0f, 1.f), nextFloat(0.f, 1.f), id, param });
	}

	// main loop ---

	srand(2 ^ 24 - 1);

	Color* colorBuffer = (Color*)fb.colorBuffer;
	//for (int j = fb.height - 1; j >= 0; j--)
	for (int j = 0; j < fb.height; j++)
	{
		for (int i = 0; i < fb.width; i++)
		{
			Pixel pixel{ i, j, {0, 0, 0, 255} };
			vec3 color(0.0f, 0.0f, 0.0f);
#if 1
			const int numSamples = 10;
			for (int s = 0; s < numSamples; s++)
			{
				Ray ray = scene.camera.generateRay(float(i) + nextFloat(0.f, 1.f), float(j) + nextFloat(0.f, 1.f));

				/*HitRecord rec;
				int index = 0;
				for (auto& sphere : scene.spheres)
				{
					if (sphere.hit(ray, rec))
					{
						//ray.maxt = rec.point.z;
						ray.maxt = rec.t;
						Material& mat = scene.materials[index];
						vec3 N = rec.normal;
						if (mat.type == 0) {
							float NdotL = -N.z;
							color += vec3(mat.r * NdotL, mat.g * NdotL, mat.b * NdotL);
						}
						else {
							vec3 V = ray.direction; V.normalize();
							vec3 R = V - N*(2.f*V.dot(N));
							// pseudo Phong
							float VdotR = R.dot(-V);
							float specular = saturate(pow(VdotR, 8.f));
							//pixel.color.fromLinear(specular, specular, specular);
							pixel.color.fromLinear(specular*mat.r, specular*mat.g, specular*mat.b);
						}
					}
					++index;
				}*/
				color += scene.color(ray, scene, 0);
			}
			color *= (1.f/float(numSamples));
#else
			Ray ray = scene.camera.generateRay(i, j);
			color += scene.color(ray, scene, 0);
#endif

			pixel.color.fromLinear(color.x, color.y, color.z);
			*colorBuffer++ = pixel.color;
		}
	}

	for (int b = 0; b < 0; b++) {
		for (int j = 1; j < fb.height; j++)
		{
			for (int i = 1; i < fb.width; i++)
			{
				Color* A = (Color*)&fb.colorBuffer[(j - 1)*fb.width * 4 + (i - 1) * 4];
				Color* B = (Color*)&fb.colorBuffer[(j - 1)*fb.width * 4 + (i) * 4];
				Color* C = (Color*)&fb.colorBuffer[(j)*fb.width * 4 + (i - 1) * 4];
				Color* D = (Color*)&fb.colorBuffer[(j)*fb.width * 4 + (i) * 4];
				vec3 H0, H1;
				H0.x = (A->R / 255.f)*0.5f + (B->R / 255.f)*0.5f;
				H0.y = (A->G / 255.f)*0.5f + (B->G / 255.f)*0.5f;
				H0.z = (A->B / 255.f)*0.5f + (B->B / 255.f)*0.5f;
				H1.x = (C->R / 255.f)*0.5f + (D->R / 255.f)*0.5f;
				H1.y = (C->G / 255.f)*0.5f + (D->G / 255.f)*0.5f;
				H1.z = (C->B / 255.f)*0.5f + (D->B / 255.f)*0.5f;
				H0.x = H0.x*0.5f + H1.x*0.5f;
				H0.y = H0.y*0.5f + H1.y*0.5f;
				H0.z = H0.z*0.5f + H1.z*0.5f;
				D->R = (uint8_t)(H0.x*255.99f);
				D->G = (uint8_t)(H0.y*255.99f);
				D->B = (uint8_t)(H0.z*255.99f);
			}
		}
	}

	// terminate ---
	
	// output image
		
	TgaFileHeader header;
	memset(&header, 0, sizeof(header));
	header.imageType = 2;
	header.width = fb.width;
	header.height = fb.height;
	header.pixelSize = 4 * 8;
	
	char filename[256];
	strcpy(filename, "output");	
	strcat(filename, ".tga");

	FILE* tgaFile = fopen(filename, "wb");	
	fwrite(&header, sizeof(header), 1, tgaFile);
	fwrite(fb.colorBuffer, fb.width*fb.height*4, 1, tgaFile);
	fclose(tgaFile);

	// free resources
	fb.destroy();

	return 0;
}