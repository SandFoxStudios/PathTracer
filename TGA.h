//
// TGA  
//
// original code : Malek Bengougam
//
// (c) Rewind Engine 2010-2016, All rights reserved.
//

#ifndef REWIND_RESOURCES_TGA_H
#define REWIND_RESOURCES_TGA_H

// ___ Includes _______________________________________________________________

enum eTgaAlpha {
	eTgaFullAlpha = 0,
	eTgaOneBitAlpha = 1,
};

enum {
	eTgaCMapType_Paletted = 1,
	eTgaDescriptorFlipBit = 1 << 5
};

enum eTgaImageType {
	eTgaImageType_Empty = 0,
	eTgaImageType_Indexed = 1,
	eTgaImageType_Rgb = 2,
	eTgaImageType_Grey = 3,
	eTgaImageType_Rle = 8		// +8 = rle compressed
};

#define TGAFILEFOOTER "TRUEVISION-XFILE"

// ___ Structs ________________________________________________________________

#pragma pack(push, 1)
struct TgaFileHeader
{
	uint8_t		idFieldSize;			// Size of ID field that follows header (0)
	uint8_t		colorMapType;			// 0 = none, 1 = paletted
	uint8_t		imageType;              // eTgaImageType
	uint16_t	colorMapStart;          // First colour map entry
	uint16_t	colorMapLength;         // Number of colors
	uint8_t 	colorMapBits;			// bits per palette entry
	uint16_t	xOrigin;
	uint16_t	yOrigin;
	uint16_t	width;
	uint16_t	height;
	uint8_t		pixelSize;
	uint8_t		descriptor;
};
#pragma pack(pop)

#endif // REWIND_RESOURCES_TGA_H
