
#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#elif defined(__APPLE__)
#include <CoreServices/CoreServices.h>
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <array>

#define TINYOBJLOADER_IMPLEMENTATION
#include "libs/tiny_obj_loader.h"

#if ONE_SHOT
#include "libs/TGA.h"
#else
#if defined(_WIN32)
#include <GL/glew.h>
#include <GL/GL.h>
#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glfw3.lib")
#pragma comment(lib, "glew32.lib")
#elif defined(__APPLE__)
#include <OpenGL/OpenGL.h>
#include <OpenGL/gl3.h>
#include <OpenGL/gl3ext.h>
#endif
#include "GLFW/glfw3.h"
#endif

#include "Framebuffer.h"
#include "Shader.h"

struct float3
{
	float x, y, z;
	inline float3() {}
    inline float3(const float scalar) : x(scalar), y(scalar), z(scalar) {}
	inline float3(const float _x, const float _y, const float _z) : x(_x), y(_y), z(_z) {}
	inline float operator[](const int i) const { return (&x)[i]; }
    inline float& operator[](const int i) { return (&x)[i]; }
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
    static float3 abs(const float3& vec) { return float3(fabs(vec.x), fabs(vec.y), fabs(vec.z)); }
};
inline float3 operator*(const float3& lhs, const float3& rhs) { return float3(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z); }
inline float3 cross(const float3& lhs, const float3& rhs) { return lhs.cross(rhs); }
inline float dot(const float3& lhs, const float3& rhs) { return lhs.dot(rhs); }
inline float3 normalize(const float3& v) { return v*(1.f/v.length()); }
typedef float3 vec3;

struct Vertex
{
    float3 position;
    float3 normal;
};

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
	float type;
	//float param; // metal(1)=fuzz, dieletric(2)=ior
};

struct HitRecord
{
	vec3 point; float t;
	vec3 normal; int material;
	vec3 uv;

	bool hit(Ray& ray, const Vertex& a, const Vertex& b, const Vertex& c)//, const vec3& na, const vec3& nb, const vec3& nc)
	{
		vec3 e1 = b.position - a.position;
		vec3 e2 = c.position - a.position;

		// Left Handed, reverse determinant
		//float det = -n.dot(ray.direction);
        vec3 q = ray.direction.cross(e2);
        float det = e1.dot(q);
        
        if (det < 0.0001f)
		//if (abs(det) < 0.0001f)
			return false;
        
		vec3 ap = (ray.origin - a.position);
        // backface culling
		//float temp = n.dot(ap);
		//if (temp < 0.f)
		//	return false;
		
		/*vec3 e = ray.direction.cross(ap);
		float v = e.dot(e2);
		if (v < 0.f || v > det)
			return false;
		float w = -e.dot(e1);
		if (w < 0.f || (v+w) > det)
			return false;*/
        
        //vec3 q = ray.direction.cross(e2);
        float u = ap.dot(q);
        if (u < 0.f || u > det)
            return false;
        vec3 r = ap.cross(e1);
        float v = ray.direction.dot(r);
        if (v < 0.f || (u + v) > det)
            return false;
        
		float ood = 1.f / det;
		float temp = e2.dot(r) * ood;
		if (temp > ray.mint && temp < ray.maxt) {
            ray.maxt = temp;
            t = temp;
			point = ray(temp);
            uv.y = u*ood; uv.z = v*ood; uv.x = 1.f - uv.y - uv.z;
            normal = a.normal * uv.x + b.normal * uv.y + c.normal * uv.z;
            //normal = e1.cross(e2);
            normal.normalize();

			return true;
		}
		return false;
	}
};

struct AABB
{
	float3 min;
	float3 max;

    void insert(const vec3& p) {
        if (p.x < min.x) min.x = p.x;
        if (p.y < min.y) min.y = p.y;
        if (p.z < min.z) min.z = p.z;
        if (p.x > max.x) max.x = p.x;
        if (p.y > max.y) max.y = p.y;
        if (p.z > max.z) max.z = p.z;
    }
    
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

/*======================== X-tests ========================*/

#define AXISTEST_X01(a, b, fa, fb)               \
p0 = a*v0.y - b*v0.z;                          \
p2 = a*v2.y - b*v2.z;                          \
if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
rad = fa * halfExtent.y + fb * halfExtent.z;   \
if(min>rad || max<-rad) return false;
#define AXISTEST_X2(a, b, fa, fb)               \
p0 = a*v0.z - b*v0.z;                       \
p1 = a*v1.z - b*v1.z;                          \
if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
rad = fa * halfExtent.y + fb * halfExtent.z;   \
if(min>rad || max<-rad) return false;
/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)               \
p0 = -a*v0.x + b*v0.z;                     \
p2 = -a*v2.x + b*v2.z;                             \
if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
rad = fa * halfExtent.x + fb * halfExtent.z;   \
if(min>rad || max<-rad) return false;
#define AXISTEST_Y1(a, b, fa, fb)               \
p0 = -a*v0.x + b*v0.z;                     \
p1 = -a*v1.x + b*v1.z;                           \
if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
rad = fa * halfExtent.x + fb * halfExtent.z;   \
if(min>rad || max<-rad) return false;
/*======================== Z-tests ========================*/
#define AXISTEST_Z12(a, b, fa, fb)               \
p1 = a*v1.x - b*v1.y;                       \
p2 = a*v2.x - b*v2.y;                          \
if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
rad = fa * halfExtent.x + fb * halfExtent.y;   \
if(min>rad || max<-rad) return false;
#define AXISTEST_Z0(a, b, fa, fb)               \
p0 = a*v0.x - b*v0.y;                   \
p1 = a*v1.x - b*v1.y;                       \
if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
rad = fa * halfExtent.x + fb * halfExtent.y;   \
if(min>rad || max<-rad) return false;

#define FINDMINMAX(x0,x1,x2,min,max) \
min = max = x0;   \
if(x1<min) min=x1;\
if(x1>max) max=x1;\
if(x2<min) min=x2;\
if(x2>max) max=x2;

struct RegularGrid
{
    struct Cell {
        AABB bounds;
        int left, right;
        std::vector<const int*> content;
        
        // adapted from http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox3.txt
        bool overlap(const int* indices, const std::vector<float3>& vertices)
        {
            float min, max, p0, p1, p2, rad;
            vec3 halfExtent = (bounds.max - bounds.min)*0.5f;
            vec3 center = (bounds.max + bounds.min)*0.5f;
            vec3 v0 = vertices[indices[0]] - center;
            vec3 v1 = vertices[indices[1]] - center;
            vec3 v2 = vertices[indices[2]] - center;
            
            vec3 e0 = v1 - v0;
            vec3 e1 = v2 - v1;
            vec3 e2 = v0 - v2;
            vec3 fe = vec3::abs(e0);
            AXISTEST_X01(e0.z, e0.y, fe.z, fe.y);
            AXISTEST_Y02(e0.z, e0.x, fe.z, fe.x);
            AXISTEST_Z12(e0.y, e0.x, fe.y, fe.x);
            fe = float3::abs(e1);
            AXISTEST_X01(e1.z, e1.y, fe.z, fe.y);
            AXISTEST_Y02(e1.z, e1.x, fe.z, fe.x);
            AXISTEST_Z0(e1.y, e1.x, fe.y, fe.x);
            fe = float3::abs(e2);
            AXISTEST_X2(e2.z, e2.y, fe.z, fe.y);
            AXISTEST_Y1(e2.z, e2.x, fe.z, fe.x);
            AXISTEST_Z12(e2.y, e2.x, fe.y, fe.x);
            
            /*  1. first test overlap in the {x,y,z}-directions */
            /*  find min, max of the triangle each direction, and test for overlap in */
            /*  that direction -- this is equivalent to testing a minimal AABB around */
            /*  the triangle against the AABB */
            
            /* test in X-direction */
            FINDMINMAX(v0.x,v1.x,v2.x,min,max);
            if(min>halfExtent.x || max<-halfExtent.x) return false;
            /* test in Y-direction */
            FINDMINMAX(v0.y,v1.y,v2.y,min,max);
            if(min>halfExtent.y || max<-halfExtent.y) return 0;
            /* test in Z-direction */
            FINDMINMAX(v0.z,v1.z,v2.z,min,max);
            if(min>halfExtent.z || max<-halfExtent.z) return 0;
            
            // 2. Plane-AABB test
            vec3 n = cross(e0, e1);
            //float e = halfExtent.x*abs(center.x) + halfExtent.y*abs(center.y) + halfExtent.z*abs(center.z);
            //float s = dot(center, n);
            vec3 vmin, vmax;
            for (int q = 0; q < 3; q++) {
                float v = v0[q];
                if(n[q] > 0.0f) {
                    vmin[q] = -halfExtent[q] - v;    // -NJMP-
                    vmax[q] = halfExtent[q] - v;    // -NJMP-
                }
                else {
                    vmin[q] = halfExtent[q] - v;    // -NJMP-
                    vmax[q] = -halfExtent[q] - v;    // -NJMP-
                }
            }
            if(dot(n, vmin) > 0.0f) return false;    // -NJMP-
            if(dot(n, vmax) < 0.0f) return false;    // -NJMP-

            content.push_back(indices);
            return true;
        }
    };
    
    AABB root;
    vec3 cellExtents;
    Cell* cells;
    int subdivisions;
    
    RegularGrid() {}
    ~RegularGrid() {
        delete[] cells;
    }
    RegularGrid(const AABB& bound, int subdiv = 1)
    {
        root = bound;
        subdivisions = subdiv;
        cellExtents = (root.max - root.min)*(1.f/(float)subdiv);
        
        cells = new Cell[subdivisions*subdivisions*subdivisions];
        // TODO: morton sequence
        for (int w = 0; w < subdiv; w++) {
            vec3 slice = root.min;
            slice.z += cellExtents.z*w;
            for (int v = 0; v < subdiv; v++) {
                vec3 row = slice;
                row.y += cellExtents.y*v;
                for (int u = 0; u < subdiv; u++) {
                    int index = (w*subdiv*subdiv) + v*subdiv + u;
                    cells[index].bounds.min = row;
                    cells[index].bounds.max = row + cellExtents;
                    row.x += cellExtents.y;
                }
            }
        }
    }

    void construct(const std::vector<float3>& vertices, const std::vector<int>& indices)
    {
        vec3 extents(1.f / cellExtents.x, 1.f / cellExtents.y, 1.f / cellExtents.z);
        int numTris = indices.size() / 3;
        int misses = 0;
        for (int tri = 0; tri < numTris; tri++)
        {
            // 1. calc bounding box and calc bounding box overlap as a crude test
            AABB triBox;
            const vec3& p0 = vertices[indices[tri*3]];
            const vec3& p1 = vertices[indices[tri*3+1]];
            const vec3& p2 = vertices[indices[tri*3+2]];
            triBox.min = p0; triBox.max = p0;
            triBox.insert(p1);
            triBox.insert(p2);
            triBox.min = (triBox.min - root.min) * extents;
            triBox.max = (triBox.max - root.min) * extents;
            int minX = (int)floor(triBox.min.x), minY = (int)floor(triBox.min.y), minZ = (int)floor(triBox.min.z);
            int maxX = (int)ceil(triBox.max.x), maxY = (int)ceil(triBox.max.y), maxZ = (int)ceil(triBox.max.z);
            // 2. refine overlapping cells
            for (int w = minZ; w < maxZ; w++) {
                for (int v = minY; v < maxY; v++) {
                    for (int u = minX; u < maxX; u++) {
                        int index = (w*subdivisions*subdivisions) + v*subdivisions + u;
                        //if (cells[index].overlap(&indices[tri*3], vertices) == false)
                        //    misses++;
                        cells[index].content.push_back(&indices[tri*3]);
                    }
                }
            }
        }
    }
    
    // in world space
    bool hit(Ray& ray, HitRecord &rec, const std::vector<float3>& vertices, const std::vector<float3>& normals) const
    {
        if (root.hit(ray, rec))
        {
            Ray endRay = ray;
            endRay.origin = rec.point + ray.direction * 50.f;
            endRay.direction = -ray.direction;
            HitRecord endRec;
            root.hit(endRay, endRec);
            
            vec3 delta = endRec.point - rec.point;
            
            vec3 extents(1.f / cellExtents.x, 1.f / cellExtents.y, 1.f / cellExtents.z);
            // uses 3DDDA
            vec3 min = (rec.point - root.min) * extents;
            int i = (int)(min.x), j = (int)(min.y), k = (int)(min.z);
            vec3 max = (endRec.point - root.min) * extents;
            int iend = (int)(max.x), jend = (int)(max.y), kend = (int)(max.z);
            // this shouldn't be necessary if we have a better max
            if (iend >= subdivisions) iend = subdivisions - 1;
            if (jend >= subdivisions) jend = subdivisions - 1;
            if (kend >= subdivisions) kend = subdivisions - 1;
            
            // todo: early exit if (ijk)end == (ijk)
            //if (abs(delta.x)<0.0001f && abs(delta.y)<0.0001f && abs(delta.z)<0.0001f)
            //    return false;
            
            // calc ray direction signs
            //int di = ray.direction.x > 0 ? 1 : (ray.direction.x < 0 ? -1 : 0);
            //int dj = ray.direction.y > 0 ? 1 : (ray.direction.y < 0 ? -1 : 0);
            //int dk = ray.direction.z > 0 ? 1 : (ray.direction.z < 0 ? -1 : 0);
            int di = (i < iend) ? 1 : ((i > iend) ? -1 : 0);
            int dj = (j < jend) ? 1 : ((j > jend) ? -1 : 0);
            int dk = (k < kend) ? 1 : ((k > kend) ? -1 : 0);
            
            float minx = floorf(min.x) * cellExtents.x, maxx = minx + cellExtents.x;
            float miny = floorf(min.y) * cellExtents.y, maxy = miny + cellExtents.y;
            float minz = floorf(min.z) * cellExtents.z, maxz = minz + cellExtents.z;
            
            vec3 diff = vec3::abs(max - min);
            float tx = ((min.x < max.x) ? (min.x - minx) : (maxx - min.x)) / diff.x;
            float ty = ((min.y < max.y) ? (min.y - miny) : (maxy - min.y)) / diff.y;
            float tz = ((min.z < max.z) ? (min.z - minz) : (maxz - min.z)) / diff.z;
            
            float deltatx = extents.x / diff.x;
            float deltaty = extents.y / diff.y;
            float deltatz = extents.z / diff.z;
            
            /*for (;;)
            {
                int index = (k*subdivisions*subdivisions) + j*subdivisions + i;
                if (cells[index].content.size()) {
                    bool hasHit = false;
                    int hitCount = 0;
                    for (const int* indices : cells[index].content) {
                        bool didHit = rec.hit(ray, positions[indices[0]], positions[indices[1]], positions[indices[2]]
                                              , normals[indices[0]], normals[indices[1]], normals[indices[2]]);
                        if (didHit) {
                            hasHit |= didHit;
                            hitCount++;
                        }
                    }
                    if (hasHit)
                        return true;
                }
                
                if (tx <= ty && tx <= tz)
                {
                    if (i == iend) break;
                    tx += deltatx;
                    i += di;
                } else if (ty <= tx && ty <= tz)
                {
                    if (j == jend) break;
                    ty += deltaty;
                    j += dj;
                } else
                {
                    if (k == kend) break;
                    tz += deltatz;
                    k += dk;
                }
            }*/
        }
        return false;
    }
};

struct Mesh
{
    std::vector<Vertex> vertices;
	//std::vector<float3> positions;
	//std::vector<float3> normals;
	std::vector<int> indices;
	AABB bounds;
    RegularGrid* grid;

    Mesh() : grid(nullptr) {}
    Mesh(Mesh&& m) : vertices(std::move(m.vertices))/*, normals(std::move(m.normals))*/, indices(std::move(m.indices)), bounds(std::move(m.bounds)) {
        grid = nullptr;
        std::swap(grid, m.grid);
    }
    ~Mesh() {
        if (grid)
            delete grid;
    }
    
	bool hit(const Ray& ray, HitRecord &rec, const int a, const int b, const int c) const
	{
		const vec3& va = vertices[a].position;
		const vec3 e1 = vertices[b].position - va;
		const vec3 e2 = vertices[c].position - va;
		
		// Left Handed, reverse determinant
		//float det = -n.dot(ray.direction);
		const vec3 q = ray.direction.cross(e2);
		float det = e1.dot(q);

		if (det < 0.0001f)
			//if (abs(det) < 0.0001f)
			return false;

		const vec3 ap = (ray.origin - va);
		// backface culling
		//float temp = n.dot(ap);
		//if (temp < 0.f)
		//	return false;

		/*vec3 e = ray.direction.cross(ap);
		float v = e.dot(e2);
		if (v < 0.f || v > det)
		return false;
		float w = -e.dot(e1);
		if (w < 0.f || (v+w) > det)
		return false;*/

		//vec3 q = ray.direction.cross(e2);
		float u = ap.dot(q);
		if (u < 0.f || u > det)
			return false;
		const vec3 r = ap.cross(e1);
		float v = ray.direction.dot(r);
		if (v < 0.f || (u + v) > det)
			return false;

		
		float ood = 1.f / det;
		float temp = e2.dot(r) * ood;
		if (temp > ray.mint && temp < ray.maxt) {
			rec.t = temp;
			rec.point = ray(temp);
			rec.uv.y = u*ood; rec.uv.z = v*ood; rec.uv.x = 1.f - rec.uv.y - rec.uv.z;
			rec.normal = vertices[a].normal * rec.uv.x + vertices[b].normal * rec.uv.y + vertices[c].normal * rec.uv.z;
			rec.normal.normalize();
			return true;
		}
		return false;
	}

	// in world space
	bool hit(Ray ray, HitRecord &rec) const
	{
		bool hasHit = false;
		//return bounds.hit(ray, rec);
#if 1//DONT_USE_GRID
        if (bounds.test(ray))
		{
			int numTri = indices.size() / 3;
			for (int i = 0; i < numTri; i++)
			{
				//bool ok = hit(ray, rec, indices[i * 3], indices[i * 3 + 1], indices[i * 3 + 2]);
                bool ok = rec.hit(ray, vertices[indices[i * 3]], vertices[indices[i * 3 + 1]], vertices[indices[i * 3 + 2]]);
                //, normals[indices[i * 3]], normals[indices[i * 3 + 1]], normals[indices[i * 3 + 2]]);
				if (ok) {
					ray.maxt = rec.t;
				}
				hasHit |= ok;
			}
        }
#else
        hasHit = grid->hit(ray, rec, vertices, normals);
#endif

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
    float3 lowerLeftCorner;
    float3 horizontal;
    float3 vertical;

	float width, height, nearPlane, farPlane, fov, aspectRatio;
	float invDenom, invTanAngle;
	
    void lookAt(const float3& eye, const float3& target, const float3& up)
    {
        position = eye;
        float3 w = normalize(target - eye);
        float3 u = normalize(cross(up, w));
        float3 v = cross(w, u);
        
        // projection
        float theta = tan(fov * M_PI / 360.0f);
        invTanAngle = 1.f / theta;
        float hh = theta;
        float hw = theta*aspectRatio;
        float focusDist = 1.f;
        lowerLeftCorner = position - hw*u*focusDist - hh*v*focusDist + w*focusDist;
        horizontal = (2.f*hw*focusDist) * u;
        vertical = (2.f*hh*focusDist) * v;
    }
    
	void initialize(float _fov, float _near, float _far, float _width, float _height)
	{
		width = _width; height = _height;
		fov = _fov;  nearPlane = _near; farPlane = _far; aspectRatio = width/height;
		invDenom = 1.f / (farPlane - nearPlane);
        
        //lookAt(eye, target, float3(0.0f, 1.0f, 0.0f));
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
        x/=width; y/=height;
        return Ray{position, (lowerLeftCorner + x*horizontal + y*vertical - position).normalize(), 0.f, farPlane - nearPlane };
		float3 offset;
		unproject(offset, x, y);
		// using same value for position & direction -> perspective ray converging to camera center
		return Ray{ position/* + offset*/, (position + offset).normalize(), 0.f, farPlane - nearPlane };
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
	scattered.origin = position + N*0.001f;
	scattered.mint = 0.00f; 
	scattered.maxt = FLT_MAX;

	float NdotL = -N.z;

	//pdf = (M_1_PI * NdotL);// 1.f / (0.5f * M_1_PI);

	if (mat.type == Material::LAMBERT)
	{
		scattered.direction = N + randomUnitVector();
		scattered.direction.normalize();

		// BRDF
		NdotL = N.dot(scattered.direction);
		
		attenuation = vec3(mat.r, mat.g, mat.b);// M_1_PI;// *pdf);
		
		//scattered.origin += N*0.001f;
		return true;
	}
	/*else if (mat.type == Material::METAL) 
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
		// TODO
        float peta = 1.55f;//mat.param;
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
	}*/
	return false;
}

vec3 Scene::color(Ray& ray, Scene& scene, int depth, uint64_t &rayCount)
{
    vec3 accumulator = vec3(0.0f);
    vec3 energy(1.0f);
    
    
    bool stop = false;
    
    for (int level = 0; level < depth; ++level)
    {

        // 1. intersect scene
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
            //Material& mat = scene.materials[0];
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
             // get position and normal at the intersection point
            vec3 P = rec.point;
            vec3 N = rec.normal;
            Material& mat = scene.materials[hitIndex];

            // compute surface lighting, and calc indirect lighting ray
            //Ray scattered;
            float pdf;
            vec3 attenuation(1.f);
#if 1
            if (scatter(mat, P, ray.direction, N, ray, attenuation, pdf))
            {
                // recurse depending on material
                //vec3 indirect = color(scattered, scene, depth + 1, rayCount);
                
                // final lighting
                //return emission + attenuation * indirect;
				//ray.origin = vec3(-1000.0, -1000.0, -1000.0);
				//ray.direction = normalize(vec3(-1000.0, -1000.0, -1000.0));
            }

			energy = energy * attenuation;
            
            bool inShadow = false;
            if (mat.type == Material::LAMBERT && -N.z > 0.f)
            {
                //rayCount++;
                //HitRecord shadowRec;
                //Ray shadowRay = {rec.point, vec3(0.f,0.f,-1.f), 0.001f, FLT_MAX };
                //for (auto& sphere : scene.spheres)
                //{
                //    if (sphere.hit(shadowRay, shadowRec))
                //    {
                //        inShadow = true;
                //        break;
                //    }
                //}
                
                // compute direct lighting
                /*const float sunSize = 0.53f;
                const float sunAngularDiameterCos = cos(sunSize*M_PI / 180.0);
                const float sunSolidAngle = 1e-5f*1000.f *19000.f*0.01f;
                
                vec3 sunSampleDir = getConeSample(vec3(0.f, 0.f, -1.f), 1.f - sunAngularDiameterCos);
                vec3 direct(0.0f);
                const float sunLight = N.dot(sunSampleDir);
                if (sunLight > 0.f && !inShadow) {
                    direct = vec3(mat.r, mat.g, mat.b)*M_1_PI * sunLight*sunSolidAngle;
                }
                //emission += direct;
                accumulator += emission * direct;*/
            }
#endif
        }
        // if nothing found, return background color
        else {
            // ambient sky light


            break;
        }
        
        if ((energy.x+energy.y+energy.z) == 0.f)
            break;
    }
    
    vec3 dir = ray.direction;
    float t = 0.5f * (dir.y + 1.0f);
    vec3 skylight = vec3(1.0 - t) +vec3(1.f, 0.7f, 0.5f) * (t);
    accumulator += energy * skylight;
    
    return accumulator;
}

// ---

void errorCallback(int error, const char* description)
{
    printf("GLFW Error %d : %s\n", error, description );
}

int main(void)
{
	// initialize ---
    const int width = 1280, height = 720;
    
    Framebuffer fb;
    fb.create(width, height);
    
    Scene scene;

#if !ONE_SHOT
    GLFWwindow* m_window;
    glfwSetErrorCallback(errorCallback);
    if (glfwInit() != GLFW_TRUE)
        throw std::runtime_error("Could not init GLFW");
    // note: MSAA is pointless here since we are not rasterizing polygons to the backbuffer
    //glfwWindowHint(GLFW_SAMPLES, MSAA);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    #ifdef _WIN32
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    #else
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    #endif
    // Why Etna : allessandro volta is italian, volta, tesla... -> Etna
    m_window = glfwCreateWindow(width, height, "Etna", NULL, NULL);
    if (m_window == NULL) {
        glfwTerminate();
        return -1;
    }
    /*glfwSetWindowUserPointer(m_window, this);
    glfwSetKeyCallback(m_window, key_callback);
    glfwSetMouseButtonCallback(m_window, mouse_button_callback);
    glfwSetCursorPosCallback(m_window, cursor_position_callback);
    glfwSetScrollCallback(m_window, wheel_callback);*/
    glfwMakeContextCurrent(m_window);
    #ifdef _WIN32
    // Initialise GLEW
    glewExperimental = true; // NÈcessaire dans le profil de base
    if (glewInit() != GLEW_OK) {
        glfwTerminate();
        throw std::runtime_error("Could not init GLEW");
    }
    glfwSwapInterval(0);
    #endif
    
	glClearColor(0.f, 0.f, 0.f, 1.f);
	glDisable(GL_BLEND);
	glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE);

    // Determine GLSL version of this context.
    //  We'll use this info to generate a GLSL shader source string
    //  with the proper version preprocessor string prepended.
    float  glLanguageVersion;
    sscanf((char *)glGetString(GL_SHADING_LANGUAGE_VERSION), "%f", &glLanguageVersion);
    GLuint version = (GLuint)(100.f * glLanguageVersion);
    Shader::SetVersionString(version);
    
    Shader m_copyShader;

    GLuint m_fakeVAO;
    GLuint m_fboID[2];
    GLuint m_depthID;
    GLuint m_textureID[2];                            // ID of the texture used for rendering
    int m_currentTexture = 0;
    GLuint m_whiteTextureID;
    
    // INIT COPY SHADER
    m_copyShader.loadVertexShader("data/shaders/postprocess.vs");
    m_copyShader.loadFragmentShader("data/shaders/copy.fs");
    m_copyShader.createProgram(Shader::ATTRIBUTE_LESS);
    
    glGenFramebuffers(2, m_fboID);
    glGenTextures(2, m_textureID);
    for (int i = 0; i < 2; i++) {
        glBindTexture(GL_TEXTURE_2D, m_textureID[i]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        //glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA32F, width, height);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, nullptr);
        glBindFramebuffer(GL_FRAMEBUFFER, m_fboID[i]);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_textureID[i], 0);
        GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
        assert(status == GL_FRAMEBUFFER_COMPLETE);
		
		glClear(GL_COLOR_BUFFER_BIT);
    }
    glGenVertexArrays(1, &m_fakeVAO);
    glBindVertexArray(m_fakeVAO);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glDisableVertexAttribArray(0);
    glBindVertexArray(0);
#endif
    
    scene.camera.position = float3(0.f, 0.f, 0.f);
    scene.camera.initialize(60.f, 0.1f, 100.f, fb.width, fb.height);
    scene.camera.lookAt(float3(0.0f, 0.0f, 0.0f), float3(0.0f, 0.0f, 1.f), float3(0.0f, 1.0f, 0.0f));
    
    // Load mesh
    
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string err;
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, "data/monkeyhead2.obj");
    scene.meshes.reserve(1);
    Mesh objMesh;
    objMesh.bounds.min = vec3(FLT_MAX);
    objMesh.bounds.max = vec3(-FLT_MAX);
    objMesh.vertices.reserve(attrib.vertices.size() / 3);
    //objMesh.positions.reserve(attrib.vertices.size() / 3);
	//objMesh.normals.reserve(attrib.vertices.size() / 3);
    const float scale = 1.0f;
    const vec3 translation(0.0f, 0.0f, 3.f);
    for (int i = 0; i < attrib.vertices.size(); i+=3) {
        const float x = attrib.vertices[i]*scale;
        const float y = attrib.vertices[i+1]*scale;
        const float z = attrib.vertices[i+2]*scale + translation.z;
        //objMesh.positions.emplace_back(vec3(x, y, z));
        if (x < objMesh.bounds.min.x) objMesh.bounds.min.x = x;
        if (y < objMesh.bounds.min.y) objMesh.bounds.min.y = y;
        if (z < objMesh.bounds.min.z) objMesh.bounds.min.z = z;
        if (x > objMesh.bounds.max.x) objMesh.bounds.max.x = x;
        if (y > objMesh.bounds.max.y) objMesh.bounds.max.y = y;
        if (z > objMesh.bounds.max.z) objMesh.bounds.max.z = z;
		const float nx = attrib.normals[i];
		const float ny = attrib.normals[i + 1];
		const float nz = attrib.normals[i + 2];
		//objMesh.normals.emplace_back(vec3(nx, ny, nz));
        objMesh.vertices.emplace_back(Vertex{{x, y, z}, {nx, ny, nz}});
    }
    for (auto &shape : shapes) {
        for (int indices = 0; indices < shape.mesh.indices.size(); indices+=3) {
            int a = shape.mesh.indices[indices].vertex_index;
            int b = shape.mesh.indices[indices+1].vertex_index;
            int c = shape.mesh.indices[indices+2].vertex_index;
            objMesh.indices.push_back(a);
            objMesh.indices.push_back(b);
            objMesh.indices.push_back(c);
        }
        //for (auto& index : shape.mesh.indices) {
        //    objMesh.indices.push_back(index.vertex_index);
        //}
    }
    //objMesh.grid = new RegularGrid(objMesh.bounds);
    //objMesh.grid->construct(objMesh.vertices, objMesh.indices);
    scene.meshes.emplace_back(std::move(objMesh));
    
    //scene.meshes.emplace_back(Mesh{ { vec3{-1.f, -1.f, 5.f}, vec3{ 0.f, 1.f, 5.f }, vec3{ 1.f, -1.f, 5.f } }, { 0, 1, 2 }, { vec3{ -1.f, -1.f, 5.f }, vec3{ 1.f, 1.f, 5.f } } });
    
#ifndef _WIN32
    srand(2 ^ 17 - 1);
#endif

	/*scene.spheres.reserve(10);
	for (int i = 0; i < 10; i++) {
		scene.spheres.emplace_back(Sphere{ { nextFloat(-3.f, +3.f), nextFloat(-3.f, +3.f), nextFloat(0.5f, +10.f) }, nextFloat(0.1f, 1.f) });
	}*/

	scene.materials.reserve(10);
    //scene.materials.emplace_back(Material{ 1.f, 1.f, 1.f, Material::LAMBERT });
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
		scene.materials.emplace_back(Material{ nextFloat(0.0f, 1.f), nextFloat(0.0f, 1.f), nextFloat(0.f, 1.f), (float)id/*, param*/ });
	}

#define USE_GLSL
#if defined(USE_GLSL)
    Shader m_pathTracingShader;
    m_pathTracingShader.loadVertexShader("data/shaders/postprocess.vs");
    m_pathTracingShader.loadFragmentShader("data/shaders/pathtracing.fs");
    m_pathTracingShader.createProgram(Shader::ATTRIBUTE_LESS);

    GLuint spheresTexture;
    glGenTextures(1, &spheresTexture);
    glBindTexture(GL_TEXTURE_2D, spheresTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, 10, 1, 0, GL_RGBA, GL_FLOAT, scene.spheres.data());
    
    GLuint materialsTexture;
    glGenTextures(1, &materialsTexture);
    glBindTexture(GL_TEXTURE_2D, materialsTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, 10, 1, 0, GL_RGBA, GL_FLOAT, scene.materials.data());
    
    GLuint indicesTexture;
    glGenTextures(1, &indicesTexture);
    glBindTexture(GL_TEXTURE_RECTANGLE, indicesTexture);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_RECTANGLE, 0, GL_RGB32UI, 1024, 1, 0, GL_RGB_INTEGER, GL_UNSIGNED_INT, scene.meshes[0].indices.data());
    //glTexSubImage2D(GL_TEXTURE_RECTANGLE, 0, 0, 0, scene.meshes[0].indices.size(), 1, GL_RGB_INTEGER, GL_UNSIGNED_INT, scene.meshes[0].indices.data());

    GLenum error = glGetError();
    assert(error == GL_NO_ERROR);
    
    GLuint verticesTexture;
    glGenTextures(1, &verticesTexture);
    glBindTexture(GL_TEXTURE_RECTANGLE, verticesTexture);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_RECTANGLE, 0, GL_RGB32F, 1024, 1, 0, GL_RGB, GL_FLOAT, scene.meshes[0].vertices.data());
    //glTexSubImage2D(GL_TEXTURE_RECTANGLE, 0, 0, 0, scene.meshes[0].vertices.size()*2, 1, GL_RGB, GL_FLOAT, scene.meshes[0].vertices.data());

    error = glGetError();
    assert(error == GL_NO_ERROR);
#else
#endif
    
	// main loop ---
    
    srand(2 ^ 24 - 1);
    
#if !ONE_SHOT
    int samples = 0;
    do {
        glfwPollEvents();
#endif
    
        samples++;
        int m_nextTexture = m_currentTexture ^ 1;
        
#ifdef _WIN32
	LARGE_INTEGER t1, t2;
	QueryPerformanceCounter(&t1);
#else
    uint64_t t1 = mach_absolute_time();
#endif
    
    uint64_t totalRayCount = 0;
	GLuint program;
#if !defined(USE_GLSL)
	
	Color* colorBuffer = (Color*)fb.colorBuffer;
	//for (int j = fb.height - 1; j >= 0; j--)

	#define TILESIZE 16
	const int numTilesW = fb.width / TILESIZE;
	const int numTilesH = fb.height / TILESIZE;

    float3* accum = (float3*)fb.accumulationBuffer;
    
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
					
#if 1
					const int numSamples = 1;// 256;
					for (int s = 0; s < numSamples; s++)
					{
						// stratify / jitter
						float sx = 1.f, sy = 1.f;
						//float sx = 0.25f*(s % 4), sy = 0.25f*(s / 4);
						Ray ray = scene.camera.generateRay(float(pixel.x) /*+ sx*nextFloat(0.f, 1.f)*/, float(pixel.y) /*+ sy*nextFloat(0.f, 1.f)*/);

						(*accum) += scene.color(ray, scene, 10, totalRayCount);
					}
					(*accum) *= (1.f / float(numSamples));
                    
                    vec3 color = (*accum)*(1.f / float(samples));
                    accum++;
#else
					Ray ray = scene.camera.generateRay(i, j);
					color += scene.color(ray, scene, 0);
#endif
					pixel.color.fromLinear(color.z, color.y, color.x);
                    
					*buffer++ = pixel.color;
				}
			}
		}
	}
#else
        program = m_pathTracingShader.getProgram();
        glUseProgram(program);
        glUniform3fv(glGetUniformLocation(program, "u_CameraPosition"), 1, &scene.camera.position.x);
        glUniform3fv(glGetUniformLocation(program, "u_CameraLowerLeftCorner"), 1, &scene.camera.lowerLeftCorner.x);
        glUniform3fv(glGetUniformLocation(program, "u_CameraHorizontal"), 1, &scene.camera.horizontal.x);
        glUniform3fv(glGetUniformLocation(program, "u_CameraVertical"), 1, &scene.camera.vertical.x);
        glUniform1i(glGetUniformLocation(program, "u_Depth"), 10);
		auto loc = glGetUniformLocation(program, "u_Samples");
		float invSamples = (float)samples;
		glUniform1f(loc, invSamples);
        
        loc = glGetUniformLocation(program, "u_NumTriangles");
        glUniform1i(loc, scene.meshes[0].indices.size()/3);
        
        glBindFramebuffer(GL_FRAMEBUFFER, m_fboID[m_currentTexture]);
        glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_2D, materialsTexture);
        glUniform1i(glGetUniformLocation(program, "u_MaterialsTexture"), 4);
        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_RECTANGLE, indicesTexture);
        glUniform1i(glGetUniformLocation(program, "u_IndicesTexture"), 3);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_RECTANGLE, verticesTexture);
        glUniform1i(glGetUniformLocation(program, "u_VerticesTexture"), 2);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, spheresTexture);
        glUniform1i(glGetUniformLocation(program, "u_SpheresTexture"), 1);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, m_textureID[m_nextTexture]);
		glUniform1i(glGetUniformLocation(program, "u_Texture"), 0);
        
        glBindVertexArray(m_fakeVAO);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glFinish();
#endif

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
	sprintf(buffer, "sample [%d] %.2fms (%.1f FPS) %.1fMrays/sec %.2fMRays/frame\n", samples, seconds*1000.f, 1.0f / seconds, totalRayCount / seconds * 1.0e-6f, totalRayCount * 1.0e-6f);
#ifdef _WIN32
	OutputDebugStringA(buffer);
#endif
    printf("%s", buffer);
    
    
#if ONE_SHOT
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
#else
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        
        program = m_copyShader.getProgram();
        glUseProgram(program);
		glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, m_textureID[m_currentTexture]);
		glUniform1i(glGetUniformLocation(program, "u_Texture"), 0);

#if defined(USE_GLSL)
        //glBindFramebuffer(GL_FRAMEBUFFER, m_nextTexture);
        loc = glGetUniformLocation(program, "u_InvNumSamples");
		invSamples = (float)samples;
        glUniform1f(loc, invSamples);
#else
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, fb.colorBuffer);
#endif

        glBindVertexArray(m_fakeVAO);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        
		m_currentTexture ^= 1;// m_nextTexture;

		//GLenum error = glGetError();
		//assert(error == GL_NO_ERROR);
		//::Sleep(1000);
		glfwSwapBuffers(m_window);
        
    } while (glfwGetKey(m_window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(m_window) == 0);
   
#if defined(USE_GLSL)
    glDeleteTextures(1, &verticesTexture);
    glDeleteTextures(1, &indicesTexture);
    glDeleteTextures(1, &materialsTexture);
    glDeleteTextures(1, &spheresTexture);
#endif
    
    m_copyShader.destroyProgram();
    
    glDeleteVertexArrays(1, &m_fakeVAO);
    //glDeleteBuffers(1, &m_fboID);
    //glDeleteTextures(1, &m_depthID);
    glDeleteFramebuffers(2, m_fboID);
    glDeleteTextures(2, m_textureID);
    
    glfwTerminate();
#endif
    
	// free resources
	fb.destroy();

	return 0;
}
