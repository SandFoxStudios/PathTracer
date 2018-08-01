
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>

#include "TGA.h"
#include "Framebuffer.h"

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#define NO_MIN_MAX
#include <windows.h>
#elif defined(__APPLE__)
#include <CoreServices/CoreServices.h>
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

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
	inline float3& normalize() { operator*=(1.f / length()); return *this; }
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

struct Material
{
	enum { LAMBERT, METAL, DIELECTRIC };
	float r, g, b;
	int type;
	float param; // metal(1)=fuzz, dieletric(2)=ior
};

struct HitRecord
{
	vec3 point; float t;
	vec3 normal; int material;

	bool hit(const Ray& ray, const vec3& a, const vec3& b, const vec3& c)
	{
		vec3 e1 = b - a;
		vec3 e2 = c - a;
        // Left Handed, reverse determinant
		vec3 n = e1.cross(e2);
		float det = -n.dot(ray.direction);
        vec3 q = ray.direction.cross(e2);
        //float det = e1.dot(q);
        
        //if (det < 0.0001f)
		if (abs(det) < 0.0001f)
			return false;
        
		vec3 ap = (ray.origin - a);
        // backface culling
		float temp = n.dot(ap);
		if (temp < 0.f)
			return false;
		
		/*vec3 e = ray.direction.cross(ap);
		float v = e.dot(e2);
		if (v < 0.f || v > det)
			return false;
		float w = -e.dot(e1);
		if (w < 0.f || (v+w) > det)
			return false;*/
        
        //vec3 q = ray.direction.cross(e2);
        float v = ap.dot(q);
        if (v < 0.f || v > det)
            return false;
        vec3 r = ap.cross(e1);
        float w = ray.direction.dot(r);
        if (w < 0.f || (v+w) > det)
            return false;
        
		float ood = 1.f / det;
		
		t = temp * ood;
		point = ray(t);
		normal = n.normalize();
		return true;
	}
};

struct AABB
{
	float3 min;
	float3 max;

	bool test(const Ray& ray) const
	{
		float tmin = ray.mint, tmax = ray.maxt;
		for (int i = 0; i < 3; i++) {
			if (abs(ray.direction[i]) < 0.001f) {
				if (ray.origin[i] < min[i] || ray.origin[i] > max[i])
					return false;
			}
			else {
				float ood = 1.f / ray.direction[i];
				float t1 = (min[i] - ray.origin[i]) * ood;
				float t2 = (max[i] - ray.origin[i]) * ood;
				if (t1 > t2) std::swap(t1, t2);
				if (t1 > tmin) tmin = t1;
				if (t2 < tmax) tmax = t2;
				if (tmin > tmax)
					return false;
			}
		}
		return true;
	}

	// in world space
	bool hit(const Ray& ray, HitRecord &rec) const
	{
		float tmin = ray.mint, tmax = ray.maxt;
		for (int i = 0; i < 3; i++) {
			if (abs(ray.direction[i]) < 0.0001f) {
				if (ray.origin[i] < min[i] || ray.origin[i] > max[i])
					return false;
			}
			else {
				float ood = 1.f / ray.direction[i];
				float t1 = (min[i] - ray.origin[i]) * ood;
				float t2 = (max[i] - ray.origin[i]) * ood;
				if (t1 > t2) std::swap(t1, t2);
				if (t1 > tmin) tmin = t1;
				if (t2 < tmax) tmax = t2;
				if (tmin > tmax)
					return false;
			}
		}
		rec.point = ray(tmin); rec.normal = (rec.point - (min+max) * 0.5f).normalize(); rec.t = tmin; return true;
	}
};

struct Mesh
{
	std::vector<float3> vertices;
	std::vector<int> indices;
	AABB bounds;

	// in world space
	bool hit(Ray ray, HitRecord &rec) const
	{
		bool hasHit = false;
		//return bounds.hit(ray, rec);
		if (bounds.test(ray))
		{
			int numTri = indices.size() / 3;
			for (int i = 0; i < numTri; i++)
			{
				bool ok = rec.hit(ray, vertices[indices[i * 3]], vertices[indices[i * 3 + 1]], vertices[indices[i * 3 + 2]]);
				if (ok) {
					ray.maxt = rec.t;
				}
				hasHit |= ok;
			}
		}
		return hasHit;
	}
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
		float B = to.dot(ray.direction);
		float C = to.dot(to) - radius*radius;
#if 0
		B *= 2.f;
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
			//return true;
			/*if (t0 <= ray.maxt && t1 >= ray.mint) {
				float thit = (t0 < ray.mint) ? t1 : t0;
				p = ray(thit);
				return thit < ray.maxt;
			}*/

		}
#else
		// shirley (optimized):
		float discriminant = B*B - A*C;
		if (discriminant > 0.f) {
			float rootDiscriminant = sqrtf(discriminant);
			float temp = (-B - rootDiscriminant) / A;
			if (temp < ray.maxt && temp > ray.mint) { rec.point = ray(temp); rec.normal = (rec.point - center) * (1.f/radius); rec.t = temp; return true; }
			float temp2 = (-B + rootDiscriminant) / A;
			if (temp2 < ray.maxt && temp2 > ray.mint) { rec.point = ray(temp2); rec.normal = (rec.point - center) * (1.f/radius); rec.t = temp2; return true; }
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
		struct { float a, b, c; };
        struct { float3 n; };
        float d;
	};

	bool hit(const Ray& ray, float3& point) const
	{
        return false;
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

struct Scene
{
	std::vector<Sphere> spheres;
	std::vector<Mesh> meshes;
	std::vector<Material> materials;
	Camera camera;

	vec3 color(Ray& ray, Scene& scene, int depth, uint64_t& rayCount);
};

// ---

float nextFloat(float min, float max)
{
	return min + ((float)rand() / (float)RAND_MAX) * (max - min);
}

float mix(float a, float b, float t)
{
	return a*(1.f - t) + b*t;
}

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
vec3 randomDirectionOnHemisphere(/*const vec3& N*/)
{
	vec3 p;
	do {
		p = vec3(nextFloat(0.0f, 1.f), nextFloat(0.0f, 1.0f), nextFloat(0.0f, 1.0f)) * 2.0f - vec3(1.f, 1.f, 1.f);
	} while (p.sqrLength() >= 1.0f/* && N.dot(p) < 0.f*/);
	return p.normalize();
}
vec3 randomUnitVector()
{
	float z = nextFloat(0.0f, 1.0f) * 2.f - 1.f;
	float a = nextFloat(0.0f, 1.0f) * 2.f * M_PI;
	float r = sqrtf(1.f - z*z);
	float x = cos(a) * r;
	float y = sin(a) * r;
	return vec3(x, y, z);
}

vec3 ortho(const vec3& v) { return abs(v.x) > abs(v.z) ? vec3(-v.y, v.x, 0.f) : vec3(0.f, -v.z, v.y); }

vec3 getConeSample(const vec3& dir, float extent)
{
	vec3 d = dir;
	d.normalize();
	vec3 o1 = ortho(d).normalize();
	vec3 o2 = d.cross(o1).normalize();
	float rx = nextFloat(0.f, 1.f);
	float ry = nextFloat(0.f, 1.f);
	rx = rx * M_PI*2.f;
	ry = 1.f - ry*extent;

	float oneminus = sqrtf(1.0f - ry*ry);
	return o1*cos(rx)*oneminus + o2*sin(rx)*oneminus + d*ry;
}

// ---

bool scatter(const Material& mat, const vec3& position, const vec3& direction, const vec3& N, Ray& scattered, vec3& attenuation, float& pdf)
{
	scattered.origin = position;
	scattered.mint = 0.001f; 
	scattered.maxt = FLT_MAX;

	float NdotL = -N.z;

	//pdf = (M_1_PI * NdotL);// 1.f / (0.5f * M_1_PI);

	if (mat.type == Material::LAMBERT)
	{
		scattered.direction = N + randomUnitVector();
		scattered.direction.normalize();

		// BRDF
		NdotL = N.dot(scattered.direction);
		
		attenuation = vec3(mat.r, mat.g, mat.b) *M_1_PI;// *pdf);
		
		//scattered.origin += N*0.001f;
		return true;
	}
	else if (mat.type == Material::METAL) 
	{	
		vec3 V = direction;
		//V.normalize();
		
		vec3 outN = N;

		vec3 R = V - outN*(2.f*V.dot(outN));
		//R.normalize();
        scattered.direction = (R);// + randomDirectionInUnitSphere()).normalize();

		attenuation = vec3(mat.r, mat.g, mat.b);// * NdotL;
		
		//scattered.origin += outN*0.001f;
		return N.dot(scattered.direction) > 0.f;
	}
	else {
		
		float peta = mat.param;
		float eta = peta;

		vec3 V = direction;

		float NdotV = V.dot(N);

		vec3 outN = N;
		float cosine = -NdotV;// / V.length();
		if (NdotV > 0.0f) {
			outN = -N;
			cosine *= -peta;
		}
		else {
			eta = 1.0f / peta;
		}

		// Reflect
		vec3 reflected = V - outN*(2.f*V.dot(outN));
		//reflected.normalize();
		
		// dielectric / refractive
		//attenuation = vec3(mat.r, mat.g, mat.b);
		attenuation = vec3(1.f, 1.f, 1.f);

		float reflect_prob = 1.0f;

		// refract
		float NdotOV = V.dot(outN);
		vec3 refracted(0.f, 0.f, 0.f);;
		float discriminant = 1.0f - eta*eta*(1.0f - NdotOV*NdotOV);

		if (discriminant > 0.0f) 
		{
			// schlick-fresnel
			float r0 = (1.f - peta) / (1.f + peta);
			r0 *= r0;
			reflect_prob = r0 + (1.f - r0)*powf(1.0f - cosine, 5.0f);

			//scattered.direction = refracted;
			//scattered.origin += outN*-0.001f;
			//return true;
			refracted = (V - outN*NdotOV)*eta - outN*sqrtf(discriminant);
			refracted.normalize();
		}

		if (nextFloat(0.f, 1.f) < reflect_prob) {
			//scattered.origin += outN*0.001f;
			scattered.direction = reflected;
		}
		else {

			scattered.direction = refracted;
			//scattered.origin += outN*-0.001f;
		}

		return true;
	}
	return false;
}

vec3 Scene::color(Ray& ray, Scene& scene, int depth, uint64_t &rayCount)
{
	HitRecord rec;
	int hitIndex = -1;
	int index = 0;
	rayCount++;
	// get closest hit
	for (auto& sphere : scene.spheres)
	{
		//Material& mat = scene.materials[index];
		//if (mat.type != Material::METAL) {
			if (sphere.hit(ray, rec))
			{
				//ray.maxt = rec.point.z;
				ray.maxt = rec.t;
				hitIndex = index;
			}
		//}
		++index;
	}
	
	for (auto& mesh : scene.meshes)
	{
		Material& mat = scene.materials[0];
		//if (mat.type != Material::METAL) {
		if (mesh.hit(ray, rec))
		{
			//ray.maxt = rec.point.z;
			ray.maxt = rec.t;
			hitIndex = 0;
		}
		//}
		++index;
	}

	if (hitIndex >= 0)
	{
		vec3 N = rec.normal;
		Material& mat = scene.materials[hitIndex];
		vec3 attenuation;
		vec3 emission(0.f, 0.f, 0.f);
		Ray scattered;
		
		
		bool inShadow = false;
		if (mat.type == Material::LAMBERT && -N.z > 0.f)
		{
			/*rayCount++;
			HitRecord shadowRec;
			Ray shadowRay = {rec.point, vec3(0.f,0.f,-1.f), 0.001f, FLT_MAX };
			for (auto& sphere : scene.spheres)
			{
				if (sphere.hit(shadowRay, shadowRec))
				{
					inShadow = true;
					break;
				}
			}*/
			
			const float sunSize = 0.53f;
			const float sunAngularDiameterCos = cos(sunSize*M_PI / 180.0);
			const float solidAngle = 1e-5f*1000.f *19000.f*0.01f;
			
			vec3 sunSampleDir = getConeSample(vec3(0.f, 0.f, -1.f), 1.f - sunAngularDiameterCos);

			const float sunLight = N.dot(sunSampleDir);
			if (sunLight > 0.f && !inShadow) {
				emission = vec3(mat.r, mat.g, mat.b)*M_1_PI * sunLight*solidAngle;
			}
		}

		float pdf;
		if (depth < 10 && scatter(mat, rec.point, ray.direction, rec.normal, scattered, attenuation, pdf))
		{
			// recurse depending on material
			return emission + attenuation * color(scattered, scene, depth + 1, rayCount);
		}
		else // fake emissive
			return emission;
	}
	else {
		// terminal ray: ambient sky light
		vec3 dir = ray.direction;
		float t = 0.5f * (dir.y + 1.0f);

		return vec3(1.0f, 1.0f, 1.0f)*(1.0 - t) +vec3(1.f, 0.7f, 0.5f) * (t);
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
	scene.camera.initialize(60.f, 0.5f, 100.f, fb.width, fb.height);

    srand(2 ^ 17 - 1);
    
	scene.spheres.reserve(10);
	for (int i = 0; i < 10; i++) {
		scene.spheres.emplace_back(Sphere{ { nextFloat(-3.f, +3.f), nextFloat(-3.f, +3.f), nextFloat(0.5f, +10.f) }, nextFloat(0.1f, 1.f) });
	}

	scene.materials.reserve(10);
	for (int i = 0; i < 10; i++) {
		int id = Material::LAMBERT;
		float param = 0.0f;
		if (i == 7) {
			id = Material::DIELECTRIC;
			param = 1.55f;
		}
		else {
			if (i < 7) {
				param = 0.0f;
			}
			else {
				id = Material::METAL;
				param = 0.0f;
			}
		}
		scene.materials.emplace_back(Material{ nextFloat(0.0f, 1.f), nextFloat(0.0f, 1.f), nextFloat(0.f, 1.f), id, param });
	}

	scene.meshes.reserve(1);
	scene.meshes.emplace_back(Mesh{ { vec3{-1.f, -1.f, 5.f}, vec3{ 0.f, 1.f, 5.f }, vec3{ 1.f, -1.f, 5.f } }, { 0, 1, 2 }, { vec3{ -1.f, -1.f, 5.f }, vec3{ 1.f, 1.f, 5.f } } });

	// main loop ---

	srand(2 ^ 24 - 1);
    
#ifdef _WIN32
	LARGE_INTEGER t1, t2;
	QueryPerformanceCounter(&t1);
#else
    uint64_t t1 = mach_absolute_time();
#endif
    
	uint64_t totalRayCount = 0;
	Color* colorBuffer = (Color*)fb.colorBuffer;
	//for (int j = fb.height - 1; j >= 0; j--)

	#define TILESIZE 16
	const int numTilesW = fb.width / TILESIZE;
	const int numTilesH = fb.height / TILESIZE;

	for (int th = 0; th < numTilesH; th++)
	{
		const int ja = TILESIZE * th;
		for (int tw = 0; tw < numTilesW; tw++)
		{
			const int ia = TILESIZE * tw;

			for (int j = 0; j < TILESIZE; j++)
			{
				Color* buffer = &colorBuffer[fb.width*(ja + j) + (ia)];
				for (int i = 0; i < TILESIZE; i++)
				{
                    uint16_t pi = (i+ia), pj = (j+ja);
					Pixel pixel{ pi, pj, {0, 0, 0, 255} };
					vec3 color(0.0f, 0.0f, 0.0f);
#if 1
					const int numSamples = 1;// 256;
					for (int s = 0; s < numSamples; s++)
					{
						// stratify / jitter
						float sx = 1.f, sy = 1.f;
						//float sx = 0.25f*(s % 4), sy = 0.25f*(s / 4);
						Ray ray = scene.camera.generateRay(float(pixel.x) + sx*nextFloat(0.f, 1.f), float(pixel.y) + sy*nextFloat(0.f, 1.f));

						color += scene.color(ray, scene, 0, totalRayCount);
					}
					color *= (1.f / float(numSamples));
#else
					Ray ray = scene.camera.generateRay(i, j);
					color += scene.color(ray, scene, 0);
#endif

					pixel.color.fromLinear(color.x, color.y, color.z);
					*buffer++ = pixel.color;
				}
			}
		}
	}

#ifdef _WIN32
	QueryPerformanceCounter(&t2);
	uint64_t dt = t2.QuadPart - t1.QuadPart;
	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);
	double seconds = double(dt) / double(freq.QuadPart);
#else
    uint64_t t2 = mach_absolute_time();
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);
    double seconds = (double)((t2-t1) * timebase.numer / timebase.denom) / 1000000000.0;
#endif
    char buffer[128];
	sprintf(buffer, "%.2fms (%.1f FPS) %.1fMrays/sec %.2fMRays/frame\n", seconds*1000.f, 1.0f / seconds, totalRayCount / seconds * 1.0e-6f, totalRayCount * 1.0e-6f);
#ifdef _WIN32
	OutputDebugStringA(buffer);
#else
    printf("%s", buffer);
#endif
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
