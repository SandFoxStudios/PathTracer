#pragma once

#include "Vector3.h"

namespace Types {
	template <typename T>
	class Ray
	{
		static_assert(std::is_floating_point<T>::value, "Template must be floating point");
	public:
        Vector3<T> origin;
        Vector3<T> direction;
        float minT;
        float maxT;
        
		Ray() {}
        Ray(Vector3<T>& o, Vector3<T>& d, float mint = 0.f, float maxt = std::numeric_limits<T>()) : origin(o), direction(d), minT(mint), maxT(maxt) {}
		~Ray() {}

		inline Vector3<T> eval(float t) const { return origin + direction*t; }

        inline const T epsilon() { return (T)1e-3f; }
	};
}
