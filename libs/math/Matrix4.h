#pragma once

#include "Vector4.h"

namespace Types {
    
	// Columns majors matrix 4x4
	template <typename T>
	class Matrix4 {
	public:
		typedef Vector4<T> Col;
		Matrix4();
		Matrix4(const Col &col1, const Col &col2, const Col &col3, const Col &col4);
		Matrix4(const T * data);
		~Matrix4();

		inline Matrix4<T> transpose() const;
		inline Matrix4<T> inverse() const;
		static Matrix4<T> identity();

        static Matrix4<T> translate(const Vector3<T>& translation);
        static Matrix4<T> rotate(float angle, const Vector3<T>& axis);

		inline float det() const;

		static Matrix4<T> TRS(const Vector3<T> &translation, const Quaternion<T> &rotation, const Vector3<T> &scale);

	public:
        template <typename U>
		friend Matrix4<U> operator*(const Matrix4<U> &lhs, const Matrix4<U> &rhs);
        template <typename U>
        friend Vector4<U> operator*(const Matrix4<U> &lhs, const Vector4<U> &rhs);

		inline bool operator==(const Matrix4<T> &c) const;
		inline bool operator!=(const Matrix4<T> &c) const;

	public:
		inline Col &operator[](const unsigned int col);
		inline const Col &operator[](const unsigned int col) const;

        inline void setCol(int col, const Col& c) { m_cols[col] = c; }
        inline Col &getCol(const unsigned int col) { return (*this)[col]; }
        inline const Col &getCol(const unsigned int col) const { return (*this)[col]; }
        
	private:
		Col m_cols[4];
	};
    
	template<typename T>
	Matrix4<T> operator*(const Matrix4<T>& lhs, const Matrix4<T>& rhs)
	{
		Matrix4<T> out = Matrix4<T>(Vector4<T>(0.f), Vector4<T>(0.f), Vector4<T>(0.f), Vector4<T>(0.f));
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				for (int k = 0; k < 4; k++)
					out[i][j] += lhs[i][k] * rhs[k][j];
		return out;
	}

	template<typename T>
	Vector4<T> operator*(const Matrix4<T> &lhs, const Vector4<T> &rhs)
	{
		return Vector4<T>(
			lhs.m_cols[0].x * rhs.x + lhs.m_cols[1].x * rhs.y + lhs.m_cols[2].x * rhs.z + lhs.m_cols[3].x * rhs.w,
			lhs.m_cols[0].y * rhs.x + lhs.m_cols[1].y * rhs.y + lhs.m_cols[2].y * rhs.z + lhs.m_cols[3].y * rhs.w,
			lhs.m_cols[0].z * rhs.x + lhs.m_cols[1].z * rhs.y + lhs.m_cols[2].z * rhs.z + lhs.m_cols[3].z * rhs.w,
			lhs.m_cols[0].w * rhs.x + lhs.m_cols[1].w * rhs.y + lhs.m_cols[2].w * rhs.z + lhs.m_cols[3].w * rhs.w
		);
	}

	template<typename T>
	inline Matrix4<T>::Matrix4()
	{
	}

	template<typename T>
	inline Matrix4<T>::Matrix4(const Col & col1, const Col & col2, const Col & col3, const Col & col4)
	{
		m_cols[0] = col1;
		m_cols[1] = col2;
		m_cols[2] = col3;
		m_cols[3] = col4;
	}

	template<typename T>
	inline Matrix4<T>::Matrix4(const T * data)
	{
		for (int col = 0; col < 4; col++)
			for (int row = 0; row < 4; row++)
				m_cols[col][row] = data[col * 4 + row];
	}

	template<typename T>
	inline Matrix4<T>::~Matrix4()
	{
	}
   
	template<typename T>
	inline Matrix4<T> Matrix4<T>::transpose() const
	{
		Matrix4<T> out;
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				out[j][i] = m_cols[i][j];
		return out;
	}
	template<typename T>
	inline Matrix4<T> Matrix4<T>::inverse() const
	{
		float d =
			m_cols[0][0] * m_cols[1][1] * m_cols[2][2] * m_cols[3][3] + m_cols[0][0] * m_cols[1][2] * m_cols[2][3] * m_cols[3][1] + m_cols[0][0] * m_cols[1][3] * m_cols[2][1] * m_cols[3][2] +
			m_cols[0][1] * m_cols[1][0] * m_cols[2][3] * m_cols[3][2] + m_cols[0][1] * m_cols[1][2] * m_cols[2][0] * m_cols[3][3] + m_cols[0][1] * m_cols[1][3] * m_cols[2][2] * m_cols[3][0] +
			m_cols[0][2] * m_cols[1][0] * m_cols[2][1] * m_cols[3][3] + m_cols[0][2] * m_cols[1][1] * m_cols[2][3] * m_cols[3][0] + m_cols[0][2] * m_cols[1][3] * m_cols[2][0] * m_cols[3][1] +
			m_cols[0][3] * m_cols[1][0] * m_cols[2][2] * m_cols[3][1] + m_cols[0][3] * m_cols[1][1] * m_cols[2][0] * m_cols[3][2] + m_cols[0][3] * m_cols[1][2] * m_cols[2][1] * m_cols[3][0] -
			m_cols[0][0] * m_cols[1][1] * m_cols[2][3] * m_cols[3][2] - m_cols[0][0] * m_cols[1][2] * m_cols[2][1] * m_cols[3][3] - m_cols[0][0] * m_cols[1][3] * m_cols[2][2] * m_cols[3][1] -
			m_cols[0][1] * m_cols[1][0] * m_cols[2][2] * m_cols[3][3] - m_cols[0][1] * m_cols[1][2] * m_cols[2][3] * m_cols[3][0] - m_cols[0][1] * m_cols[1][3] * m_cols[2][0] * m_cols[3][2] -
			m_cols[0][2] * m_cols[1][0] * m_cols[2][3] * m_cols[3][1] - m_cols[0][2] * m_cols[1][1] * m_cols[2][0] * m_cols[3][3] - m_cols[0][2] * m_cols[1][3] * m_cols[2][1] * m_cols[3][0] -
			m_cols[0][3] * m_cols[1][0] * m_cols[2][1] * m_cols[3][2] - m_cols[0][3] * m_cols[1][1] * m_cols[2][2] * m_cols[3][0] - m_cols[0][3] * m_cols[1][2] * m_cols[2][0] * m_cols[3][1];

		if (d == 0.f)
		{
			std::cerr << "Cannot inverse matrix" << std::endl;
			return (*this);
		}
		d = 1.f / d;

		Matrix4<T> m;
		m[0][0] = m_cols[1][1] * m_cols[2][2] * m_cols[3][3] * d + m_cols[1][2] * m_cols[2][3] * m_cols[3][1] * d + m_cols[1][3] * m_cols[2][1] * m_cols[3][2] * d - m_cols[1][1] * m_cols[2][3] * m_cols[3][2] * d - m_cols[1][2] * m_cols[2][1] * m_cols[3][3] * d - m_cols[1][3] * m_cols[2][2] * m_cols[3][1] * d;
		m[0][1] = m_cols[0][1] * m_cols[2][3] * m_cols[3][2] * d + m_cols[0][2] * m_cols[2][1] * m_cols[3][3] * d + m_cols[0][3] * m_cols[2][2] * m_cols[3][1] * d - m_cols[0][1] * m_cols[2][2] * m_cols[3][3] * d - m_cols[0][2] * m_cols[2][3] * m_cols[3][1] * d - m_cols[0][3] * m_cols[2][1] * m_cols[3][2] * d;
		m[0][2] = m_cols[0][1] * m_cols[1][2] * m_cols[3][3] * d + m_cols[0][2] * m_cols[1][3] * m_cols[3][1] * d + m_cols[0][3] * m_cols[1][1] * m_cols[3][2] * d - m_cols[0][1] * m_cols[1][3] * m_cols[3][2] * d - m_cols[0][2] * m_cols[1][1] * m_cols[3][3] * d - m_cols[0][3] * m_cols[1][2] * m_cols[3][1] * d;
		m[0][3] = m_cols[0][1] * m_cols[1][3] * m_cols[2][2] * d + m_cols[0][2] * m_cols[1][1] * m_cols[2][3] * d + m_cols[0][3] * m_cols[1][2] * m_cols[2][1] * d - m_cols[0][1] * m_cols[1][2] * m_cols[2][3] * d - m_cols[0][2] * m_cols[1][3] * m_cols[2][1] * d - m_cols[0][3] * m_cols[1][1] * m_cols[2][2] * d;
		m[1][0] = m_cols[1][0] * m_cols[2][3] * m_cols[3][2] * d + m_cols[1][2] * m_cols[2][0] * m_cols[3][3] * d + m_cols[1][3] * m_cols[2][2] * m_cols[3][0] * d - m_cols[1][0] * m_cols[2][2] * m_cols[3][3] * d - m_cols[1][2] * m_cols[2][3] * m_cols[3][0] * d - m_cols[1][3] * m_cols[2][0] * m_cols[3][2] * d;
		m[1][1] = m_cols[0][0] * m_cols[2][2] * m_cols[3][3] * d + m_cols[0][2] * m_cols[2][3] * m_cols[3][0] * d + m_cols[0][3] * m_cols[2][0] * m_cols[3][2] * d - m_cols[0][0] * m_cols[2][3] * m_cols[3][2] * d - m_cols[0][2] * m_cols[2][0] * m_cols[3][3] * d - m_cols[0][3] * m_cols[2][2] * m_cols[3][0] * d;
		m[1][2] = m_cols[0][0] * m_cols[1][3] * m_cols[3][2] * d + m_cols[0][2] * m_cols[1][0] * m_cols[3][3] * d + m_cols[0][3] * m_cols[1][2] * m_cols[3][0] * d - m_cols[0][0] * m_cols[1][2] * m_cols[3][3] * d - m_cols[0][2] * m_cols[1][3] * m_cols[3][0] * d - m_cols[0][3] * m_cols[1][0] * m_cols[3][2] * d;
		m[1][3] = m_cols[0][0] * m_cols[1][2] * m_cols[2][3] * d + m_cols[0][2] * m_cols[1][3] * m_cols[2][0] * d + m_cols[0][3] * m_cols[1][0] * m_cols[2][2] * d - m_cols[0][0] * m_cols[1][3] * m_cols[2][2] * d - m_cols[0][2] * m_cols[1][0] * m_cols[2][3] * d - m_cols[0][3] * m_cols[1][2] * m_cols[2][0] * d;
		m[2][0] = m_cols[1][0] * m_cols[2][1] * m_cols[3][3] * d + m_cols[1][1] * m_cols[2][3] * m_cols[3][0] * d + m_cols[1][3] * m_cols[2][0] * m_cols[3][1] * d - m_cols[1][0] * m_cols[2][3] * m_cols[3][1] * d - m_cols[1][1] * m_cols[2][0] * m_cols[3][3] * d - m_cols[1][3] * m_cols[2][1] * m_cols[3][0] * d;
		m[2][1] = m_cols[0][0] * m_cols[2][3] * m_cols[3][1] * d + m_cols[0][1] * m_cols[2][0] * m_cols[3][3] * d + m_cols[0][3] * m_cols[2][1] * m_cols[3][0] * d - m_cols[0][0] * m_cols[2][1] * m_cols[3][3] * d - m_cols[0][1] * m_cols[2][3] * m_cols[3][0] * d - m_cols[0][3] * m_cols[2][0] * m_cols[3][1] * d;
		m[2][2] = m_cols[0][0] * m_cols[1][1] * m_cols[3][3] * d + m_cols[0][1] * m_cols[1][3] * m_cols[3][0] * d + m_cols[0][3] * m_cols[1][0] * m_cols[3][1] * d - m_cols[0][0] * m_cols[1][3] * m_cols[3][1] * d - m_cols[0][1] * m_cols[1][0] * m_cols[3][3] * d - m_cols[0][3] * m_cols[1][1] * m_cols[3][0] * d;
		m[2][3] = m_cols[0][0] * m_cols[1][3] * m_cols[2][1] * d + m_cols[0][1] * m_cols[1][0] * m_cols[2][3] * d + m_cols[0][3] * m_cols[1][1] * m_cols[2][0] * d - m_cols[0][0] * m_cols[1][1] * m_cols[2][3] * d - m_cols[0][1] * m_cols[1][3] * m_cols[2][0] * d - m_cols[0][3] * m_cols[1][0] * m_cols[2][1] * d;
		m[3][0] = m_cols[1][0] * m_cols[2][2] * m_cols[3][1] * d + m_cols[1][1] * m_cols[2][0] * m_cols[3][2] * d + m_cols[1][2] * m_cols[2][1] * m_cols[3][0] * d - m_cols[1][0] * m_cols[2][1] * m_cols[3][2] * d - m_cols[1][1] * m_cols[2][2] * m_cols[3][0] * d - m_cols[1][2] * m_cols[2][0] * m_cols[3][1] * d;
		m[3][1] = m_cols[0][0] * m_cols[2][1] * m_cols[3][2] * d + m_cols[0][1] * m_cols[2][2] * m_cols[3][0] * d + m_cols[0][2] * m_cols[2][0] * m_cols[3][1] * d - m_cols[0][0] * m_cols[2][2] * m_cols[3][1] * d - m_cols[0][1] * m_cols[2][0] * m_cols[3][2] * d - m_cols[0][2] * m_cols[2][1] * m_cols[3][0] * d;
		m[3][2] = m_cols[0][0] * m_cols[1][2] * m_cols[3][1] * d + m_cols[0][1] * m_cols[1][0] * m_cols[3][2] * d + m_cols[0][2] * m_cols[1][1] * m_cols[3][0] * d - m_cols[0][0] * m_cols[1][1] * m_cols[3][2] * d - m_cols[0][1] * m_cols[1][2] * m_cols[3][0] * d - m_cols[0][2] * m_cols[1][0] * m_cols[3][1] * d;
		m[3][3] = m_cols[0][0] * m_cols[1][1] * m_cols[2][2] * d + m_cols[0][1] * m_cols[1][2] * m_cols[2][0] * d + m_cols[0][2] * m_cols[1][0] * m_cols[2][1] * d - m_cols[0][0] * m_cols[1][2] * m_cols[2][1] * d - m_cols[0][1] * m_cols[1][0] * m_cols[2][2] * d - m_cols[0][2] * m_cols[1][1] * m_cols[2][0] * d;

		return m;
	}
	template<typename T>
	inline Matrix4<T> Matrix4<T>::identity()
	{
		return Matrix4<T>(
			Col(1, 0, 0, 0),
			Col(0, 1, 0, 0),
			Col(0, 0, 1, 0),
			Col(0, 0, 0, 1)
			);
	}
	template<typename T>
	inline float Matrix4<T>::det() const
	{
		return
			m_cols[0][3] * m_cols[1][2] * m_cols[2][1] * m_cols[3][0] - m_cols[0][2] * m_cols[1][3] * m_cols[2][1] * m_cols[3][0] -
			m_cols[0][3] * m_cols[1][1] * m_cols[2][2] * m_cols[3][0] + m_cols[0][1] * m_cols[1][3] * m_cols[2][2] * m_cols[3][0] +
			m_cols[0][2] * m_cols[1][1] * m_cols[2][3] * m_cols[3][0] - m_cols[0][1] * m_cols[1][2] * m_cols[2][3] * m_cols[3][0] -
			m_cols[0][3] * m_cols[1][2] * m_cols[2][0] * m_cols[3][1] + m_cols[0][2] * m_cols[1][3] * m_cols[2][0] * m_cols[3][1] +
			m_cols[0][3] * m_cols[1][0] * m_cols[2][2] * m_cols[3][1] - m_cols[0][0] * m_cols[1][3] * m_cols[2][2] * m_cols[3][1] -
			m_cols[0][2] * m_cols[1][0] * m_cols[2][3] * m_cols[3][1] + m_cols[0][0] * m_cols[1][2] * m_cols[2][3] * m_cols[3][1] +
			m_cols[0][3] * m_cols[1][1] * m_cols[2][0] * m_cols[3][2] - m_cols[0][1] * m_cols[1][3] * m_cols[2][0] * m_cols[3][2] -
			m_cols[0][3] * m_cols[1][0] * m_cols[2][1] * m_cols[3][2] + m_cols[0][0] * m_cols[1][3] * m_cols[2][1] * m_cols[3][2] +
			m_cols[0][1] * m_cols[1][0] * m_cols[2][3] * m_cols[3][2] - m_cols[0][0] * m_cols[1][1] * m_cols[2][3] * m_cols[3][2] -
			m_cols[0][2] * m_cols[1][1] * m_cols[2][0] * m_cols[3][3] + m_cols[0][1] * m_cols[1][2] * m_cols[2][0] * m_cols[3][3] +
			m_cols[0][2] * m_cols[1][0] * m_cols[2][1] * m_cols[3][3] - m_cols[0][0] * m_cols[1][2] * m_cols[2][1] * m_cols[3][3] -
			m_cols[0][1] * m_cols[1][0] * m_cols[2][2] * m_cols[3][3] + m_cols[0][0] * m_cols[1][1] * m_cols[2][2] * m_cols[3][3];
	}
	template<typename T>
	inline Matrix4<T> Matrix4<T>::TRS(const Vector3<T>& translation, const Quaternion<T> & rotation, const Vector3<T>& scale)
	{
		Matrix4<T> t = Matrix4<T>::identity();
		t[3] = toVec4(translation, 1.f);
		Matrix4<T> r = toMat4(rotation);
		Matrix4<T> s = Matrix4<T>::identity();
		for (int i = 0; i < 3; i++)
			s[i][i] = scale[i];
		return t * (r * s);
	}
    template<typename T>
    inline Matrix4<T> Matrix4<T>::translate(const Vector3<T>& translation)
    {
        Matrix4<T> t = Matrix4<T>::identity();
        t[3] = toVec4(translation, 1.f);
        return t;
    }
    template<typename T>
    inline Matrix4<T> Matrix4<T>::rotate(float angle, const Vector3<T>& axis)
    {
        Quaternion<T> rotation;
        Matrix4<T> r = toMat4(rotation);
        return r;
    }
	template<typename T>
	inline bool Matrix4<T>::operator==(const Matrix4<T>& c) const
	{
		return (m_cols[0] == c.m_cols[0] && m_cols[1] == c.m_cols[1] && m_cols[2] == c.m_cols[2] && m_cols[3] == c.m_cols[3]);
	}
	template<typename T>
	inline bool Matrix4<T>::operator!=(const Matrix4<T>& c) const
	{
		return (m_cols[0] != c.m_cols[0] || m_cols[1] != c.m_cols[1] || m_cols[2] != c.m_cols[2] || m_cols[3] != c.m_cols[3]);
	}

	template<typename T>
	inline Vector4<T>& Matrix4<T>::operator[](const unsigned int col)
	{
		return m_cols[col];
	}

	template<typename T>
	inline const Vector4<T>& Matrix4<T>::operator[](const unsigned int col) const
	{
		return m_cols[col];
	}

}
