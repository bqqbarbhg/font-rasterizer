#pragma once

#include "common.h"
#include <vector>

struct Line
{
	Vec2f a, b;
};

struct Bezier
{
	Vec2f a, b, c;
};

struct Font
{
	std::vector<Line> lines;
	std::vector<Bezier> beziers;
};

struct RasterizeOptions
{
	Vec2f offset;
	Vec2f scale;
};

void rasterizeFont(const Font &font, const RasterizeOptions &options, uint32_t width, uint32_t height, float *distances);
