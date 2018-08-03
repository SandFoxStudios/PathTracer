#pragma once

#include "Vector3.h"

namespace Types {

	// Columns majors matrix 4x4
	template <typename T>
	class Matrix3 {
	public:
		typedef Vector3<T> Col;
		Matrix3();
		Matrix3(const Col &col1, const Col &col2, const Col &col3);
		Matrix3(const T * data);
		~Matrix3();

		inline Matrix3<T> transpose() const;
		inline Matrix3<T> inverse() const;
		static Matrix3<T> identity();

	public:
        template <typename U>
		friend Matrix3<U> operator*(const Matrix3<U> &lhs, const Matrix3<U> &rhs);
        template <typename U>
        friend Vector3<U> operator*(const Matrix3<U> &lhs, const Vector3<U> &rhs);

		inline bool operator==(const Matrix3<T> &c) const;
		inline bool operator!=(const Matrix3<T> &c) const;

	public:
		inline Col &operator[](const unsigned int col);
		inline const Col &operator[](const unsigned int col) const;

	private:
		Col m_cols[3];
	};

	template<typename T>
	inline Matrix3<T> operator*(const Matrix3<T>& lhs, const Matrix3<T>& rhs)
	{
		Matrix3<T> out;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					out[i][j] += lhs[i][k] * rhs[k][j];
		return out;
	}

	template<typename T>
	inline Vector3<T> operator*(const Matrix3<T> &lhs, const Vector3<T> &rhs)
	{
		return Vector3<T>(
			lhs.m_cols[0].x * rhs.x + lhs.m_cols[1].x * rhs.y + lhs.m_cols[2].x * rhs.z,
			lhs.m_cols[0].y * rhs.x + lhs.m_cols[1].y * rhs.y + lhs.m_cols[2].y * rhs.z,
			lhs.m_cols[0].z * rhs.x + lhs.m_cols[1].z * rhs.y + lhs.m_cols[2].z * rhs.z
		);
	}

	template<typename T>
	inline Matrix3<T>::Matrix3()
	{
	}

	template<typename T>
	inline Matrix3<T>::Matrix3(const Col & col1, const Col & col2, const Col & col3)
	{
		m_cols[0] = col1;
		m_cols[1] = col2;
		m_cols[2] = col3;
	}

	template<typename T>
	inline Matrix3<T>::Matrix3(const T * data)
	{
		for (int col = 0; col < 3; col++)
			for (int row = 0; row < 3; row++)
				m_cols[col][row] = data[col * 3 + row];
	}

	template<typename T>
	inline Matrix3<T>::~Matrix3()
	{
	}

	template<typename T>
	inline Matrix3<T> Matrix3<T>::transpose() const
	{
		Matrix3<T> out;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				out[j][i] = m_cols[i][j];
		return out;
	}
	template<typename T>
	inline Matrix3<T> Matrix3<T>::inverse() const
	{
		double det = m_cols[0][0] * (m_cols[1][1] * m_cols[2][2] - m_cols[2][1] * m_cols[1][2]) -
			m_cols[0][1] * (m_cols[1][0] * m_cols[2][2] - m_cols[1][2] * m_cols[2][0]) +
			m_cols[0][2] * (m_cols[1][0] * m_cols[2][1] - m_cols[1][1] * m_cols[2][0]);
		
		if (det == 0.f)
		{
			std::cerr << "Cannot inverse matrix" << std::endl;
			return (*this);
		}
		det = 1.f / det;

		Matrix3<T> m; // inverse of matrix m
		m[0][0] = (m_cols[1][1] * m_cols[2][2] - m_cols[2][1] * m_cols[1][2]) * det;
		m[0][1] = (m_cols[0][2] * m_cols[2][1] - m_cols[0][1] * m_cols[2][2]) * det;
		m[0][2] = (m_cols[0][1] * m_cols[1][2] - m_cols[0][2] * m_cols[1][1]) * det;
		m[1][0] = (m_cols[1][2] * m_cols[2][0] - m_cols[1][0] * m_cols[2][2]) * det;
		m[1][1] = (m_cols[0][0] * m_cols[2][2] - m_cols[0][2] * m_cols[2][0]) * det;
		m[1][2] = (m_cols[1][0] * m_cols[0][2] - m_cols[0][0] * m_cols[1][2]) * det;
		m[2][0] = (m_cols[1][0] * m_cols[2][1] - m_cols[2][0] * m_cols[1][1]) * det;
		m[2][1] = (m_cols[2][0] * m_cols[0][1] - m_cols[0][0] * m_cols[2][1]) * det;
		m[2][2] = (m_cols[0][0] * m_cols[1][1] - m_cols[1][0] * m_cols[0][1]) * det;
		return m;
	}
	template<typename T>
	inline Matrix3<T> Matrix3<T>::identity()
	{
		return Matrix3<T>(
			Col(1, 0, 0),
			Col(0, 1, 0),
			Col(0, 0, 1)
			);
	}
	template<typename T>
	inline bool Matrix3<T>::operator==(const Matrix3<T>& c) const
	{
		return (m_cols[0] == c.m_cols[0] && m_cols[1] == c.m_cols[1] && m_cols[2] == c.m_cols[2]);
	}
	template<typename T>
	inline bool Matrix3<T>::operator!=(const Matrix3<T>& c) const
	{
		return (m_cols[0] != c.m_cols[0] || m_cols[1] != c.m_cols[1] || m_cols[2] != c.m_cols[2]);
	}

	template<typename T>
	inline Vector3<T>& Matrix3<T>::operator[](const unsigned int col)
	{
		return m_cols[col];
	}

	template<typename T>
	inline const Vector3<T>& Matrix3<T>::operator[](const unsigned int col) const
	{
		return m_cols[col];
	}

}
