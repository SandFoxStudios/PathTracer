#pragma once

#include <iostream>

namespace Types {

	template <typename T>
	class Vector2 {
	public:
		union {
			struct {
				T x, y;
			};
			T data[2];
		};

		Vector2();
		explicit Vector2(T a);
		Vector2(T a1, T a2);
		template <typename U>
		Vector2(const Vector2<U> &vec);

		// comparison
		bool operator<(const Vector2<T> &v) const;
		bool operator==(const Vector2<T> &v) const;
		bool operator!=(const Vector2<T> &v) const;
		T& operator[](unsigned int idx);
		const T& operator[](unsigned int idx) const;


		// operations
		float length() const;
		Vector2<T> reverse() const;
		static float dot(const Vector2<T> &lhs, const Vector2<T> &rhs);

		// arithmetic operator
        template <typename U>
		friend Vector2<U> operator*(const Vector2<U> &vec, const float scalar);
        template <typename U>
		friend Vector2<U> operator+(const Vector2<U>& lhs, const Vector2<U>& rhs);
        template <typename U>
		friend Vector2<U> operator-(const Vector2<U>& lhs, const Vector2<U>& rhs);
        template <typename U>
        friend std::ostream& operator <<(std::ostream& os, const Vector2<U>& vec);
	};
    
	//------------------------
	//--- Class Definition ---
	//------------------------
	template <typename T>
	inline Vector2<T>::Vector2() : Vector2(T(), T())
	{
	}
	template <typename T>
	inline Vector2<T>::Vector2(T a) : Vector2(a, a)
	{
	}
	template <typename T>
	inline Vector2<T>::Vector2(T a1, T a2) : x(a1), y(a2)
	{
	}
	template<typename T>
	template<typename U>
	inline Vector2<T>::Vector2(const Vector2<U> &vec) : x(static_cast<T>(vec.x)), y(static_cast<T>(vec.y))
	{
	}

	template <typename T>
	inline bool Vector2<T>::operator<(const Vector2<T> &v) const {
		for (int i = 0; i < 2; i++)
		{
			if (data[i] < v.data[i])
				return true;
			if (data[i] > v.data[i])
				return false;
		}
		return false;
	}

	template <typename T>
	inline bool Vector2<T>::operator==(const Vector2<T> &v) const {
		return (x == v.x && y == v.y);
	}

	template <typename T>
	inline bool Vector2<T>::operator!=(const Vector2<T> &v) const {
		return (x != v.x || y != v.y);
	}
	template <typename T>
	inline T& Vector2<T>::operator[](unsigned int idx)
	{
		assert(idx < 2);
		return data[idx];
	}
	template <typename T>
	inline const T& Vector2<T>::operator[](unsigned int idx) const
	{
		assert(idx < 2);
		return data[idx];
	}
	template <typename T>
	inline float Vector2<T>::length() const {
		return sqrtf(x*x + y*y);
	}
	template <typename T>
	inline Vector2<T> Vector2<T>::reverse() const {
		return Vector2<T>(-x, -y);
	}
	template <typename T>
	float Vector2<T>::dot(const Vector2<T> &lhs, const Vector2<T> &rhs)
	{
		return lhs.x * rhs.x + lhs.y * rhs.y;
	}
	template <typename T>
	Vector2<T> operator*(const Vector2<T> &vec, const float scalar)
	{
		return Vector2<T>(vec.x * scalar, vec.y * scalar);
	}
	template <typename T>
	std::ostream& operator<<(std::ostream& os, const Vector2<T>& vec)
	{
		os << "(" << vec.x << ", " << vec.y << ")";
		return os;
	}
	template <typename T>
	Vector2<T> operator+(const Vector2<T>& lhs, const Vector2<T>& rhs)
	{
		return Vector2<T>(lhs.x + rhs.x, lhs.y + rhs.y);
	}
	template <typename T>
	Vector2<T> operator-(const Vector2<T>& lhs, const Vector2<T>& rhs)
	{
		return Vector2<T>(lhs.x - rhs.x, lhs.y - rhs.y);
	}

}

