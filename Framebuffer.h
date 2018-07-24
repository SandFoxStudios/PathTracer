#pragma once

inline float saturate(float f)
{
	if (f < 0.f)
		f = 0.f;
	else if (f > 1.f)
		f = 1.f;
	return f;
}

struct Color
{
	uint8_t R, G, B, A;

	/*void fromLinear(float r, float g, float b) {
		float sRGB = 1.0f / 2.2f;
		R = uint8_t(255.99f * pow(saturate(r), sRGB));
		G = uint8_t(255.99f * pow(saturate(g), sRGB));
		B = uint8_t(255.99f * pow(saturate(b), sRGB));
		A = 255;
	}*/
	inline float toSRGB(float x) { return  1.055f * pow(x, 0.416666667) - 0.055f; }
	void fromLinear(float r, float g, float b) {
		float x = toSRGB(r > 0.0f ? r : 0.0f);
		float y = toSRGB(g > 0.0f ? g : 0.0f);
		float z = toSRGB(b > 0.0f ? b : 0.0f);
		R = uint8_t(255.99f * saturate(x));
		G = uint8_t(255.99f * saturate(y));
		B = uint8_t(255.99f * saturate(z));
		A = 255;
	}
};

struct Pixel
{
	uint16_t x, y;
	Color color;
};

struct Framebuffer
{
	uint16_t width;
	uint16_t height;
	uint8_t* colorBuffer;

	void create(int w, int h) {
		width = w; height = h;
		colorBuffer = new uint8_t[width * height * sizeof(Color)];
		memset(colorBuffer, 0, width * height * sizeof(Color));
	}
	
	void destroy() {
		delete[] colorBuffer;
	}
};
