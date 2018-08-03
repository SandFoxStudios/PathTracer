#pragma once

#include <iostream>

namespace Types {

	template <typename T>
	class Vector4 {
	public:
		union {
			struct {
				T x, y, z, w;
			};
			T data[4];
		};

		Vector4();
		explicit Vector4(T a);
		Vector4(T a1, T a2, T a3, T a4);
		template<typename U>
		Vector4(const Vector4<U> &vec);

		// comparisons operator
		bool operator<(const Vector4<T> &v) const;
		bool operator==(const Vector4<T> &c) const;
		bool operator!=(const Vector4<T> &c) const;

		T& operator[](unsigned int idx);
		const T& operator[](unsigned int idx) const;
		// operations
		float length() const;
		Vector4<T> reverse() const;
		Vector4<T> normalize() const;
		static T dot(const Vector4<T> &lhs, const Vector4<T> &rhs);

		// arithmetic operator
		friend Vector4<T> operator+(const Vector4<T>& lhs, const Vector4<T>& rhs);
		friend Vector4<T> operator-(const Vector4<T>& lhs, const Vector4<T>& rhs);
        friend std::ostream& operator<<(std::ostream& os, const Vector4<T>& c);
	};

	//------------------------
	//--- Class Definition ---
	//------------------------
	template <typename T>
	inline Vector4<T>::Vector4() : Vector4(T(), T(), T(), T())
	{
	}
	template <typename T>
	inline Vector4<T>::Vector4(T a) : Vector4(a, a, a, a)
	{
	}
	template <typename T>
	inline Vector4<T>::Vector4(T a1, T a2, T a3, T a4) : x(a1), y(a2), z(a3), w(a4)
	{
	}
	template<typename T>
	template<typename U>
	inline Vector4<T>::Vector4(const Vector4<U> &vec) : 
		x(static_cast<T>(vec.x)), y(static_cast<T>(vec.y)), z(static_cast<T>(vec.z)), w(static_cast<T>(vec.w))
	{
	}

	template <typename T>
	inline bool Vector4<T>::operator<(const Vector4<T> &v) const
	{
		for (int i = 0; i < 4; i++)
		{
			if (data[i] < v.data[i])
				return true;
			if (data[i] > v.data[i])
				return false;
		}
		return false;
	}

	template <typename T>
	inline bool Vector4<T>::operator==(const Vector4<T> &c) const
	{
		return (x == c.x && y == c.y && z == c.z && w == c.w);
	}
	template <typename T>
	inline bool Vector4<T>::operator!=(const Vector4<T> &c) const
	{
		return (x != c.x || y != c.y || z != c.z || w != c.w);
	}

	template <typename T>
	inline T& Vector4<T>::operator[](unsigned int idx)
	{
		assert(idx < 4);
		return data[idx];
	}
	template <typename T>
	inline const T& Vector4<T>::operator[](unsigned int idx) const
	{
		assert(idx < 4);
		return data[idx];
	}
	template <typename T>
	inline float Vector4<T>::length() const
	{
		return sqrtf(x*x + y*y + z*z + w*w);
	}
	template <typename T>
	inline Vector4<T> Vector4<T>::reverse() const
	{
		return Vector4<T>(-x, -y, -z, -w);
	}
	template <typename T>
	inline Vector4<T> Vector4<T>::normalize() const
	{
		float l = length();
		return Vector4<T>(x / l, y / l, z / l, w / l);
	}
	template <typename T>
	T Vector4<T>::dot(const Vector4<T> &lhs, const Vector4<T> &rhs)
	{
		return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w;
	}
	template <typename T>
	std::ostream& operator<<(std::ostream& os, const Vector4<T>& c)
	{
		os << "(" << c.r << ", " << c.g << ", " << c.b << ", " << c.a << ")";
		return os;
	}
	template <typename T>
	Vector4<T> operator+(const Vector4<T>& lhs, const Vector4<T>& rhs)
	{
		return Vector4<T>(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
	}
	template <typename T>
	Vector4<T> operator-(const Vector4<T>& lhs, const Vector4<T>& rhs)
	{
		return Vector4<T>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
	}
}
