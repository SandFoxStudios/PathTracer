#pragma once

#include <cassert>
#include <cfloat>
#include "Vector2.h"
#include "Vector3.h"
#include "Vector4.h"
#include "Quaternion.h"
#include "Matrix3.h"
#include "Matrix4.h"
#include "Ray.h"

typedef Types::Vector4<unsigned char> Color32;
typedef Types::Vector4<float> Color4f;
typedef Types::Vector2<float> Vector2f;
typedef Types::Vector3<float> Vector3f;
typedef Types::Vector3<double> Vector3d;
typedef Types::Vector4<float> Vector4f;
typedef Types::Quaternion<float> Quaternionf;
typedef Types::Matrix4<float> Matrix4f;
typedef Types::Matrix3<float> Matrix3f;

namespace Types {

	template <typename T>
	inline T clamp(T value, T max, T min)
	{
		if (value > max)
			return max;
		if (value < min)
			return min;
		return value;
	}

	template <typename T>
	inline Vector3<T> toVec3(const Vector2<T> &vec, float z)
	{
		return Vector3<T>(vec.x, vec.y, z);
	}

	template <typename T>
	inline Vector3<T> toVec3(const Vector4<T> &vec)
	{
		return Vector3<T>(vec.x, vec.y, vec.z);
	}

	template <typename T>
	inline Vector4<T> toVec4(const Vector3<T> &vec, float w)
	{
		return Vector4<T>(vec.x, vec.y, vec.z, w);
	}

	inline Color4f toColor4f(const Color32 &color)
	{
		return Color4f(
			color.x / 255.f, 
			color.y / 255.f, 
			color.z / 255.f, 
			color.w / 255.f
		);
	}

	inline Color32 toColor32(const Color4f &color)
	{
		return Color32(
			static_cast<unsigned char>(clamp(color.x * 255, 0.f, 255.f)), 
			static_cast<unsigned char>(clamp(color.y * 255, 0.f, 255.f)),
			static_cast<unsigned char>(clamp(color.z * 255, 0.f, 255.f)),
			static_cast<unsigned char>(clamp(color.w * 255, 0.f, 255.f))
		);
	}

    template <typename T>
    inline Vector3<T> expf(const Vector3<T> &vec)
    {
        return Vector3<T>(::expf(vec.x), ::expf(vec.y), ::expf(vec.z));
    }
    
	template <typename T>
	inline Vector3<T> normalize(const Vector3<T> &vec)
	{
		return vec / vec.length();
	}

    template <typename T>
    inline float dot(const Vector3<T> &u, const Vector3<T> &v)
    {
        return Vector3<T>::dot(u, v);
    }
    
    template <typename T>
    inline Vector3<T> cross(const Vector3<T> &u, const Vector3<T> &v)
    {
        return Vector3<T>::cross(u, v);
    }

	template <typename T>
	inline Matrix3<T> toMat3(const Matrix4<T> &mat)
	{
		return Matrix3<T>(
			toVec3(mat[0]),
			toVec3(mat[1]),
			toVec3(mat[2])
		);
	}

	template <typename T>
	inline Matrix4<T> toMat4(const Quaternion<T> &quat)
	{
		Matrix4<T> m = Matrix4<T>::identity();
		T sqw = quat.w*quat.w;
		T sqx = quat.x*quat.x;
		T sqy = quat.y*quat.y;
		T sqz = quat.z*quat.z;

		// invs (inverse square length) is only required if quaternion is not already normalised
		T invs = 1 / (sqx + sqy + sqz + sqw);
		m[1][1] = ( sqx - sqy - sqz + sqw)*invs; // since sqw + sqx + sqy + sqz =1/invs*invs
		m[1][1] = (-sqx + sqy - sqz + sqw)*invs;
		m[2][2] = (-sqx - sqy + sqz + sqw)*invs;

		T tmp1 = quat.x * quat.y;
		T tmp2 = quat.z * quat.w;
		m[1][0] = T(2.0) * (tmp1 + tmp2)*invs;
		m[0][1] = T(2.0) * (tmp1 - tmp2)*invs;

		tmp1 = quat.x * quat.z;
		tmp2 = quat.y * quat.w;
		m[2][0] = T(2.0) * (tmp1 - tmp2)*invs;
		m[0][2] = T(2.0) * (tmp1 + tmp2)*invs;
		tmp1 = quat.y * quat.z;
		tmp2 = quat.x * quat.w;
		m[2][1] = T(2.0) * (tmp1 + tmp2)*invs;
		m[1][2] = T(2.0) * (tmp1 - tmp2)*invs;

		m[0][3] = T(0.0);
		m[1][3] = T(0.0);
		m[2][3] = T(0.0);
		m[3][0] = T(0.0);
		m[3][1] = T(0.0);
		m[3][2] = T(0.0);
		m[3][3] = T(1.0);
		return m;
	}

}
