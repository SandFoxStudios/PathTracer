#pragma once

#include "Vector3.h"
#include "Vector4.h"

namespace Types {

	template <typename T>
	class Quaternion : public Vector4<T>
	{
		static_assert(std::is_floating_point<T>::value, "Template must be floating point");
	public:
		Quaternion() {}
		Quaternion(T a1, T a2, T a3, T a4) : Vector4<T>(a1, a2, a3, a4) {}
		~Quaternion() {}

        static Quaternion<T> fromAxisAngle(const T angle, const Vector3<T>& axis);
		static Quaternion<T> fromEuler(const T x, const T y, const T z);

		Vector3<T> euler();

		static Quaternion<T> identity();
	};

    template<typename T>
    inline Quaternion<T> Quaternion<T>::fromAxisAngle(const T angle, const Vector3<T>& axis)
    {
        Quaternion<T> q;
        q.x = axis.x * sin(angle/2.f);
        q.y = axis.y * sin(angle/2.f);
        q.z = axis.z * sin(angle/2.f);
        q.w = cos(angle/2.f);
        return q;
    }
    
	template<typename T>
	inline Quaternion<T> Quaternion<T>::fromEuler(const T yaw, const T pitch, const T roll)
	{
		Quaternion<T> quat;
		yaw *= 0.5;
		pitch *= 0.5;
		roll *= 0.5;

		double c1 = cos(yaw);
		double c2 = cos(pitch);
		double c3 = cos(roll);
		double s1 = sin(yaw);
		double s2 = sin(pitch);
		double s3 = sin(roll);

		double c1c2 = c1*c2;
		double s1s2 = s1*s2;

		quat.w = c1c2*c3 - s1s2*s3;
		quat.x = c1c2*s3 + s1s2*c3;
		quat.y = s1*c2*c3 + c1*s2*s3;
		quat.z = c1*s2*c3 - s1*c2*s3;
		return quat;
	}

	template<typename T>
	inline Vector3<T> Quaternion<T>::euler()
	{
        // strange behavior of xcode, xyzw are not visible!
        T x = this->x, y = this->y, z = this->z, w = this->w;
		Vector3<T> vec;
		vec.x = atan2(2 * (w * x + y * z), 1 - 2 * (x*x + y*y));
		vec.y = asin(2 * (w * y - z * x));
		vec.z = atan2(2 * (w * z + x * y), 1 - 2 * (y*y + z*z));
		return vec;
	}

	template<typename T>
	inline Quaternion<T> Quaternion<T>::identity()
	{
		Quaternion<T> quat = Quaternion<T>(0.f, 0.f, 0.f, 1.f);
		return quat;
	}
}
