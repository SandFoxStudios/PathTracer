#pragma once

#include <iostream>
#include <cmath>

namespace Types {

	template <typename T>
	class Vector3 {
	public:
		union {
			struct {
				T x, y, z;
			};
			T data[3];
		};

		Vector3();
		explicit Vector3(T a);
		Vector3(T a1, T a2, T a3);
		template<typename U>
		Vector3(const Vector3<U> &vec);

		// comparison
		bool operator<(const Vector3<T> &v) const;
		bool operator==(const Vector3<T> &v) const;
		bool operator!=(const Vector3<T> &v) const;

		T& operator[](unsigned int idx);
		const T& operator[](unsigned int idx) const;

		// operations
		T length() const;
        Vector3<T> operator-() const;
		Vector3<T> reverse() const;
        Vector3<T> normalize() const;

		float dot(const Vector3<T> &vec) const;
        Vector3<T> cross(const Vector3<T> &vec) const;

		static T dot(const Vector3<T> &lhs, const Vector3<T> &rhs);
		static Vector3<T> cross(const Vector3<T> &lhs, const Vector3<T> &rhs);

		// arithmetic operator
        void operator*=(float scalar);
        
        template <typename U>
        friend Vector3<U> operator*(const Vector3<U> &u, const Vector3<U> &v);
        template <typename U>
        friend Vector3<U> operator*(float scalar, const Vector3<U> &vec);
        template <typename U>
        friend Vector3<U> operator*(const Vector3<U> &vec, float scalar);
        template <typename U>
		friend Vector3<U> operator/(const Vector3<U> &vec, float scalar);
        template <typename U>
		friend Vector3<U> operator+(const Vector3<U>& lhs, const Vector3<U>& rhs);
        template <typename U>
		friend Vector3<U> operator-(const Vector3<U>& lhs, const Vector3<U>& rhs);
        template <typename U>
		friend std::ostream& operator <<(std::ostream& os, const Vector3<U>& vec);
	};

	//------------------------
	//--- Class Definition ---
	//------------------------
	template<typename T>
	inline Vector3<T>::Vector3() : Vector3<T>(T(), T(), T())
	{
	}
	template<typename T>
	inline Vector3<T>::Vector3(T a) : Vector3<T>(a, a, a)
	{
	}
	template<typename T>
	inline Vector3<T>::Vector3(T a1, T a2, T a3) : x(a1), y(a2), z(a3)
	{
	}
	template<typename T>
	template<typename U>
	inline Vector3<T>::Vector3(const Vector3<U> &vec) : x(static_cast<T>(vec.x)), y(static_cast<T>(vec.y)), z(static_cast<T>(vec.z))
	{
	}

	template<typename T>
	inline bool Vector3<T>::operator<(const Vector3<T> &v) const {
		for (int i = 0; i < 3; i++)
		{
			if (data[i] < v.data[i])
				return true;
			if (data[i] > v.data[i])
				return false;
		}
		return false;
	}
	template<typename T>
	inline bool Vector3<T>::operator==(const Vector3<T> &v) const {
		return (x == v.x && y == v.y && z == v.z);
	}
	template<typename T>
	inline bool Vector3<T>::operator!=(const Vector3<T> &v) const {
		return (x != v.x || y != v.y || z != v.z);
	}
	template<typename T>
	inline T& Vector3<T>::operator[](unsigned int idx)
	{
		assert(idx < 3);
		return data[idx];
	}
	template<typename T>
	inline const T& Vector3<T>::operator[](unsigned int idx) const
	{
		assert(idx < 3);
		return data[idx];
	}
	template<typename T>
	inline T Vector3<T>::length() const
	{
		return sqrt(x*x + y*y + z*z);
	}
	template<typename T>
	inline Vector3<T> Vector3<T>::reverse() const
	{
		return Vector3<T>(-x, -y, -z);
	}
    template<typename T>
    inline Vector3<T> Vector3<T>::operator-() const
    {
        return reverse();
    }
	template<typename T>
	Vector3<T> Vector3<T>::normalize() const
	{
		return (*this) / length();
	}
	template<typename T>
	T Vector3<T>::dot(const Vector3<T> &lhs, const Vector3<T> &rhs)
	{
		return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
	}
	template<typename T>
	Vector3<T> Vector3<T>::cross(const Vector3<T> &lhs, const Vector3<T> &rhs)
	{
		return Vector3<T>(
			lhs.y * rhs.z - lhs.z * rhs.y,
			lhs.z * rhs.x - lhs.x * rhs.z,
			lhs.x * rhs.y - lhs.y * rhs.x
		);
	}
    template<typename T>
    inline Vector3<T> Vector3<T>::cross(const Vector3<T> &vec) const
    {
        return Vector3<T>(y * vec.z - z * vec.y, z * vec.x - x * vec.z, x * vec.y - y * vec.x);
    }

    template<typename T>
    inline void Vector3<T>::operator*=(float scalar)
    {
        x *= scalar; y *= scalar; z *= scalar;
    }
    
	template <typename T>
	Vector3<T> operator*(const Vector3<T> &vec, float scalar)
	{
		return Vector3<T>(vec.x * scalar, vec.y * scalar, vec.z * scalar);
	}
    
    template <typename T>
    Vector3<T> operator*(float scalar, const Vector3<T> &vec)
    {
        return vec * scalar;
    }
    
    template <typename T>
    Vector3<T> operator*(const Vector3<T> &u, const Vector3<T> &v)
    {
        return Vector3<T>(u.x * v.x, u.y * v.y, u.z * v.z);
    }

	template <typename T>
	Vector3<T> operator/(const Vector3<T> &vec, float scalar)
	{
		return Vector3<T>(vec.x / scalar, vec.y / scalar, vec.z / scalar);
	}

	template <typename T>
	std::ostream& operator<<(std::ostream& os, const Vector3<T>& vec)
	{
		os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
		return os;
	}

	template <typename T>
	Vector3<T> operator+(const Vector3<T>& lhs, const Vector3<T>& rhs)
	{
		return Vector3<T>(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
	}

    template <typename T>
	Vector3<T> operator-(const Vector3<T>& lhs, const Vector3<T>& rhs)
	{
		return Vector3<T>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
	}

}
