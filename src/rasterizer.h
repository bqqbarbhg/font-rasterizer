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

struct Bounds
{
	Vec2f boundsMin, boundsMax;

	template <typename T>
	T distSq(const Vec2<T> &point) const {
		Vec2<T> delta = point - clamp(point, boundsMin.cast<T>(), boundsMax.cast<T>());
		return dot2(delta);
	}
};

struct BvhLeaf
{
	enum {
		Line,
		Bezier,
	};

	Bounds bounds;
	uint16_t itemIndex;
	uint8_t itemType;
};

struct BvhNode
{
	Bounds bounds;
	uint32_t childIndex = 0;
	uint16_t childCount = 0;
	bool childLeaf = false;
};

struct FontBvh
{
	std::vector<BvhLeaf> leafs;
	std::vector<BvhNode> nodes;
};

struct RasterizeOptions
{
	Vec2f offset;
	Vec2f scale;
};

void rasterizeFont(const Font &font, const RasterizeOptions &options, uint32_t width, uint32_t height, float *distances);
