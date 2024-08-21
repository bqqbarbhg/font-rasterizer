#include "rasterizer.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>

// -- Implementation

uint64_t bezierCalls;

template <typename T, typename U>
static bool intersectLine(T y, const Vec2<U> &A, const Vec2<U> &B, float* hits)
{
	T p0 = T(A.y) - y, p1 = T(B.y) - y;
	T a = p1 - p0;
	T b = p0;

	T ax = T(B.x - A.x);
	T bx = T(A.x);

	T t = -b / a;
	auto hit = abs(a) >= 0.0001f & (t >= 0.0f) & (t <= 1.0f);
	vectorStore(hits, maskSelect(hit, ax*t + bx, T(Inf)));
	return anyTrue(hit);
}

template <typename T, typename U>
static bool intersectBezier(T y, const Vec2<U> &A, const Vec2<U> &B, const Vec2<U> &C, float *hits)
{
	T p0 = T(A.y) - y, p1 = T(B.y) - y, p2 = T(C.y) - y;
	T a = p0 - p1 * 2.0f + p2;
	T b = (p1 - p0) * 2.0f;
	T c = p0;

	T ax = T(A.x - B.x * 2.0f + C.x);
	T bx = T((B.x - A.x) * 2.0f);
	T cx = T(A.x);

	bool anyHit = false;
	T miss = Inf;
	T hit0 = miss, hit1 = miss;

	auto valid = abs(a) > 0.0f;
	if (anyTrue(valid)) {
		T disc2 = b*b - a*c*4.0f;

		auto exists = valid & (disc2 >= 0.0f);
		if (anyTrue(exists)) {
			T disc = sqrt(disc2);

			T rcpA = T(-0.5f) / a;
			T t0 = (b + disc) * rcpA;
			T t1 = (b - disc) * rcpA;

			auto h0 = exists & (t0 >= 0.0f) & (t0 <= 1.0f);
			auto h1 = exists & (t1 >= 0.0f) & (t1 <= 1.0f);

			hit0 = maskSelect(h0, (ax*t0+bx)*t0+cx, miss);
			hit1 = maskSelect(h1, (ax*t1+bx)*t1+cx, miss);
			anyHit = anyTrue(h0 | h1);
		}
	}
	if (anyFalse(valid)) {
		T t = c / -b;
		auto h = !valid & (t >= 0.0f) & (t <= 1.0f);
		hit0 = maskSelect(h, (ax*t+bx)*t+cx, hit0);
		anyHit |= anyTrue(h);
	}

	if (anyHit) {
		vectorStore(hits + 0, hit0);
		vectorStore(hits + vectorWidth<T>(), hit1);
		return true;
	}
	return false;
}

// https://www.shadertoy.com/view/WltSD7
template <typename T>
T cosAcos3(T x) {
	x = sqrt(max(x * 0.5f + 0.5f, 0.0f));
	return x*(x*(x*(x*-0.008972f+0.039071f)-0.107074f)+0.576975f)+0.5f;
}

template <typename T>
T distanceToLine(const Vec2<T> &pos, const Vec2<T> &A, const Vec2<T> &B)
{
    Vec2<T> pa = pos - A, ba = B - A;
    T h = saturate(dot(pa,ba)/dot(ba,ba));
    return dot2(pa - ba*h);
}

// https://www.shadertoy.com/view/MlKcDD
template <typename T>
T distanceToBezier(const Vec2<T> &pos, const Vec2<T> &A, const Vec2<T> &B, const Vec2<T> &C)
{
	bezierCalls++;

	Vec2<T> a = B - A;
	Vec2<T> b = A - B*2.0f + C;
	Vec2<T> c = a * 2.0f;
	Vec2<T> d = A - pos;

	T kk = T(1.0f) / (T(1e-7f) + dot(b, b));
	T kx = kk * dot(a, b);
	T ky = kk * (dot(a, a)*2.0f + dot(d, b)) * (1.0f/3.0f);
	T kz = kk * dot(d, a);

	T p = ky - kx*kx;
	T q = kx*(kx*kx*2.0f - ky*3.0f) + kz;
	T p3 = p*p*p;
	T q2 = q*q;
	T h = q2 + p3*4.0f;

	T resTrue = Inf, resFalse = Inf;
	auto mask = h >= 0.0f;
	if (anyTrue(mask)) {
		h = sqrt(h);

		Vec2<T> x = vec2(h - q, -h - q) * 0.5f;
		T uvx = approxCbrt(x.x);
		T uvy = approxCbrt(x.y);
		T t = uvx + uvy;

		t = t - (t*(t*t+p*3.0f)+q)/(t*t*3.0f+p*3.0f);

		t = saturate(t - kx);
		Vec2<T> w = d + (c + b*t)*t;
		resTrue = dot2(w);
	}
	if (anyFalse(mask)) {
		T z = sqrt(max(-p, 0.0f));
		T m = cosAcos3(q / (p*z*2.0f));
		T n = sqrt(T(1.0f) - m*m);
		n = n * sqrt(3.0f);
		T tx = saturate((m+m)*z-kx);
		T ty = saturate((-n-m)*z-kx);
		Vec2<T> qx = d + (c + b*tx)*tx;
		Vec2<T> qy = d + (c + b*ty)*ty;
		resFalse = min(dot2(qx), dot2(qy));
	}

	return maskSelect(mask, resTrue, resFalse);
}

template <typename T>
struct PrecomputedBezier
{
	Vec2<T> a, b, c;
	Vec2<T> dOffset;
	T kk, kx;
	Vec2<T> kyDot;
	T kyBias;
	Vec2<T> kzDot;
};

template <typename T>
void precomputeBezier(PrecomputedBezier<T> &r, const Vec2<T> &A, const Vec2<T> &B, const Vec2<T> &C)
{
	r.dOffset = A.cast<T>();
	r.a = (B - A).cast<T>();
	r.b = (A - B*2.0f + C).cast<T>();
	r.c = r.a * 2.0f;
	r.kk = T(1.0f) / (T(1e-7f) + dot(r.b, r.b));
	r.kx = r.kk * dot(r.a, r.b);
	r.kyDot = r.b * r.kk * (1.0f / 3.0f);
	r.kyBias = r.kk * dot(r.a, r.a) * (2.0f / 3.0f);
	r.kzDot = r.a * r.kk;
}

// https://www.shadertoy.com/view/MlKcDD
template <typename T>
T distanceToBezier(const PrecomputedBezier<T> &bezier, const Vec2<T> &pos)
{
	bezierCalls++;

	Vec2<T> d = bezier.dOffset - pos;

	T kk = bezier.kk;
	T kx = bezier.kx;
	T ky = dot(bezier.kyDot, d) + bezier.kyBias;
	T kz = dot(bezier.kzDot, d);

	T p = ky - kx*kx;
	T q = kx*(kx*kx*2.0f - ky*3.0f) + kz;
	T p3 = p*p*p;
	T q2 = q*q;
	T h = q2 + p3*4.0f;

	T resTrue = Inf, resFalse = Inf;
	auto mask = h >= 0.0f;
	if (anyTrue(mask)) {
		h = sqrt(h);

		Vec2<T> x = vec2(h - q, -h - q) * 0.5f;
		T uvx = approxCbrt(x.x);
		T uvy = approxCbrt(x.y);
		T t = uvx + uvy;

		// t = t - (t*(t*t+p*3.0f)+q)/(t*t*3.0f+p*3.0f);

		t = saturate(t - kx);
		Vec2<T> w = d + (bezier.c + bezier.b*t)*t;
		resTrue = dot2(w);
	}
	if (anyFalse(mask)) {
		T z = sqrt(max(-p, 0.0f));
		T m = cosAcos3(q / (p*z*2.0f));
		T n = sqrt(T(1.0f) - m*m);
		n = n * sqrt(3.0f);
		T tx = saturate((m+m)*z-kx);
		T ty = saturate((-n-m)*z-kx);
		Vec2<T> qx = d + (bezier.c + bezier.b*tx)*tx;
		Vec2<T> qy = d + (bezier.c + bezier.b*ty)*ty;
		resFalse = min(dot2(qx), dot2(qy));
	}

	return maskSelect(mask, resTrue, resFalse);
}

struct Rasterizer
{
	const Font &font;
	const RasterizeOptions &options;
	uint32_t width;
	uint32_t height;
	float *distances;
	FontBvh bvh;
};

template <typename T>
struct Precomputed
{
	std::vector<PrecomputedBezier<T>> beziers;
};

template <typename T>
inline T transformX(Rasterizer &r, T x) {
	return x*r.options.scale.x + r.options.offset.x;
}

template <typename T>
inline T transformY(Rasterizer &r, T y) {
	return y*r.options.scale.y + r.options.offset.y;
}

template <typename T>
void setupWinding(Rasterizer &r)
{
	constexpr uint32_t step = vectorWidth<T>();
	const float epsilonUp =   1.000001f;
	const float epsilonDown = 0.999999f;

	uint32_t countX[step];
	float hitX[step][256];

	for (uint32_t py = 0; py < r.height; py += step) {
		float hits[2 * step];
		float init[step];
		for (uint32_t dy = 0; dy < step; dy++) {
			init[dy] = float(py + dy);
			countX[dy] = 0;
		}

		T y = transformY(r, vectorLoad<T>(init));

		// HACK: Adjust y not to lie at control points
		T base = floor(y);
		y = clamp(y, base * epsilonUp, (base + 1.0f) * epsilonDown);

		for (const Line &line : r.font.lines) {
			float hits[step];
			if (intersectLine(y, line.a, line.b, hits)) {
				for (uint32_t i = 0; i < step; i++) {
					hitX[i][countX[i]] = hits[i];
					countX[i] += hits[i] < Inf ? 1 : 0;
				}
			}
		}
		for (const Bezier &bezier : r.font.beziers) {
			if (intersectBezier(y, bezier.a, bezier.b, bezier.c, hits)) {
				for (uint32_t i = 0; i < step; i++) {
					hitX[i][countX[i]] = hits[i];
					countX[i] += hits[i] < Inf ? 1 : 0;
					hitX[i][countX[i]] = hits[i + step];
					countX[i] += hits[i + step] < Inf ? 1 : 0;
				}
			}
		}

		for (uint32_t dy = 0; dy < step; dy++) {
			uint32_t hitCount = countX[dy];
			float *hits = hitX[dy];

			std::sort(hits, hits + hitCount);
			hits[hitCount] = Inf;

			uint32_t dstY = py + dy;

			uint32_t index = 0;
			bool inside = false;
			for (uint32_t px = 0; px < r.width; px++) {
				float x = transformX(r, float(px));
				while (x >= hits[index])
				{
					inside = !inside;
					index++;
				}

				r.distances[dstY * r.width + px] = inside ? 1.0f : -1.0f;
			}
		}
	}
}

template <typename T>
void calculateDistance(Rasterizer &r)
{
	constexpr uint32_t step = vectorWidth<T>();

	for (uint32_t py = 0; py < r.height; py += step) {
		float init[step];
		for (uint32_t dy = 0; dy < step; dy++) {
			init[dy] = float(py + dy);
		}
		T y = transformY(r, vectorLoad<T>(init));

		for (uint32_t px = 0; px < r.width; px++) {
			float x = transformX(r, float(px));
			Vec2<T> pos = vec2(T(x), y);

			T dist2 = Inf;
			for (const Line &line : r.font.lines) {
				T d = distanceToLine(pos, line.a.cast<T>(), line.b.cast<T>());
				dist2 = min(dist2, d);
			}
			for (const Bezier &bezier : r.font.beziers) {
				T d = distanceToBezier(pos, bezier.a.cast<T>(), bezier.b.cast<T>(), bezier.c.cast<T>());
				dist2 = min(dist2, d);
			}

			vectorStore(init, sqrt(dist2));
			for (uint32_t dy = 0; dy < step; dy++) {
				uint32_t dstY = py + dy;
				r.distances[dstY * r.width + px] *= init[dy];
			}
		}
	}
}

template <typename T>
void precompute(Precomputed<T> &pre, const Font &font)
{
	pre.beziers.resize(font.beziers.size());
	for (size_t i = 0; i < font.beziers.size(); i++) {
		const Bezier &bezier = font.beziers[i];
		Vec2<T> a = bezier.a.cast<T>();
		Vec2<T> b = bezier.b.cast<T>();
		Vec2<T> c = bezier.c.cast<T>();
		precomputeBezier(pre.beziers[i], a, b, c);
	}
}

template <typename T>
void calculateDistancePre(Rasterizer &r, Precomputed<T> &pre)
{
	constexpr uint32_t step = vectorWidth<T>();

	for (uint32_t py = 0; py < r.height; py += step) {
		float init[step];
		for (uint32_t dy = 0; dy < step; dy++) {
			init[dy] = float(py + dy);
		}
		T y = transformY(r, vectorLoad<T>(init));

		for (uint32_t px = 0; px < r.width; px++) {
			float x = transformX(r, float(px));
			Vec2<T> pos = vec2(T(x), y);

			T dist2 = Inf;
			for (const Line &line : r.font.lines) {
				T d = distanceToLine(pos, line.a.cast<T>(), line.b.cast<T>());
				dist2 = min(dist2, d);
			}
			for (const PrecomputedBezier<T> &bezier : pre.beziers) {
				T d = distanceToBezier(bezier, pos);
				dist2 = min(dist2, d);
			}

			vectorStore(init, sqrt(dist2));
			for (uint32_t dy = 0; dy < step; dy++) {
				uint32_t dstY = py + dy;
				r.distances[dstY * r.width + px] *= init[dy];
			}
		}
	}
}

inline float project(const BvhLeaf &node, uint32_t axis) {
	return ((float*)&node.bounds.boundsMin)[axis] + ((float*)&node.bounds.boundsMax)[axis];
}

void buildBvhNode(FontBvh &bvh, uint32_t dstIndex, uint32_t startIndex, uint32_t count)
{
	BvhLeaf *leafs = bvh.leafs.data() + startIndex;

	Vec2f boundsMin = vec2(Inf, Inf);
	Vec2f boundsMax = vec2(-Inf, -Inf);
	for (uint32_t i = 0; i < count; i++) {
		boundsMin = min(boundsMin, leafs[i].bounds.boundsMin);
		boundsMax = max(boundsMax, leafs[i].bounds.boundsMax);
	}

	BvhNode &dst = bvh.nodes[dstIndex];
	dst.bounds.boundsMin = boundsMin;
	dst.bounds.boundsMax = boundsMax;
	if (count <= 8) {
		dst.childIndex = startIndex;
		dst.childCount = count;
		dst.childLeaf = true;
		return;
	}

	uint32_t childIndex = (uint32_t)bvh.nodes.size();
	dst.childIndex = childIndex;
	dst.childCount = 2;
	dst.childLeaf = false;

	bvh.nodes.emplace_back();
	bvh.nodes.emplace_back();

	Vec2f boundsExtent = boundsMax - boundsMin;
	uint32_t axis = boundsExtent.x > boundsExtent.y ? 0 : 1;

	std::sort(leafs, leafs + count, [=](const BvhLeaf &a, const BvhLeaf &b) {
		return project(a, axis) < project(b, axis);
	});

	uint32_t midpoint = count / 2;
	buildBvhNode(bvh, childIndex + 0, startIndex, midpoint);
	buildBvhNode(bvh, childIndex + 1, startIndex + midpoint, count - midpoint);
}

void buildBvh(FontBvh &bvh, const Font &font)
{
	uint32_t lineIndex = 0;
	for (const Line &line : font.lines) {
		bvh.leafs.push_back(BvhLeaf{
			min(line.a, line.b),
			max(line.a, line.b),
			uint16_t(lineIndex),
			BvhLeaf::Line,
		});
		lineIndex++;
	}

	uint32_t bezierIndex = 0;
	for (const Bezier &bezier : font.beziers) {
		bvh.leafs.push_back(BvhLeaf{
			min(min(bezier.a, bezier.b), bezier.c),
			max(max(bezier.a, bezier.b), bezier.c),
			uint16_t(bezierIndex),
			BvhLeaf::Bezier,
		});
		bezierIndex++;
	}

	bvh.nodes.emplace_back();
	buildBvhNode(bvh, 0, 0, uint32_t(bvh.leafs.size()));
}

template <typename T>
void calculateDistanceBvh(Rasterizer &r, Precomputed<T> &pre)
{
	constexpr uint32_t step = vectorWidth<T>();

	struct StackEntry {
		const BvhNode *node;
		T dist2;
	};

	StackEntry stack[64];

	for (uint32_t py = 0; py < r.height; py += step) {
		float init[step];
		for (uint32_t dy = 0; dy < step; dy++) {
			init[dy] = float(py + dy);
		}
		T y = transformY(r, vectorLoad<T>(init));

		T prevDist2 = Inf;

		for (uint32_t px = 0; px < r.width; px++) {
			float x = transformX(r, float(px));
			Vec2<T> pos = vec2(T(x), y);

			T dist2 = prevDist2;

			uint32_t stackCount = 1;
			stack[0].node = r.bvh.nodes.data();
			stack[0].dist2 = 0.0f;
			
			while (stackCount > 0) {
				StackEntry entry = stack[--stackCount];
				if (!anyTrue(entry.dist2 < dist2)) continue;

				const BvhNode &node = *entry.node;
				if (node.childLeaf) {
					for (uint32_t i = 0; i < node.childCount; i++) {
						const BvhLeaf &leaf = r.bvh.leafs[node.childIndex + i];
						T leafDist2 = leaf.bounds.distSq(pos);
						if (!anyTrue(leafDist2 < dist2)) continue;

						if (leaf.itemType == BvhLeaf::Line) {
							const Line &line = r.font.lines[leaf.itemIndex];
							T d = distanceToLine(pos, line.a.cast<T>(), line.b.cast<T>());
							dist2 = min(dist2, d);
						} else if (leaf.itemType == BvhLeaf::Bezier) {
#if 0
							const PrecomputedBezier<T> &bezier = pre.beziers[leaf.itemIndex];
							T d = distanceToBezier(bezier, pos);
#else
							const Bezier &bezier = r.font.beziers[leaf.itemIndex];
							T d = distanceToBezier(pos, bezier.a.cast<T>(), bezier.b.cast<T>(), bezier.c.cast<T>());
#endif
							dist2 = min(dist2, d);
						}
					}
					continue;
				}

				const BvhNode *childA = &r.bvh.nodes[node.childIndex + 0];
				const BvhNode *childB = &r.bvh.nodes[node.childIndex + 1];
				T distA = childA->bounds.distSq(pos);
				T distB = childB->bounds.distSq(pos);
				if (vectorSum(distB) > vectorSum(distA)) {
					std::swap(childA, childB);
					std::swap(distA, distB);
				}

				if (anyTrue(distA < dist2)) {
					stack[stackCount++] = { childA, distA };
					if (anyTrue(distB < dist2)) {
						stack[stackCount++] = { childB, distB };
					}
				}
			}

			T dist = sqrt(dist2);

			// TODO: Why does this not work??
			T nextDist = dist + T(1.0f);
			// prevDist2 = nextDist*nextDist;

			vectorStore(init, dist);
			for (uint32_t dy = 0; dy < step; dy++) {
				uint32_t dstY = py + dy;
				r.distances[dstY * r.width + px] *= init[dy];
			}
		}
	}
}

template <typename T>
void renderSdf(Rasterizer &r)
{
	bezierCalls = 0;

	buildBvh(r.bvh, r.font);
	setupWinding<T>(r);
	Precomputed<T> pre;
	precompute(pre, r.font);
	// calculateDistance<T>(r);
	calculateDistanceBvh<T>(r, pre);

	//calculateDistancePre(r, pre);

	printf(">>> %llu (%.2f/px)\n", bezierCalls, (double)bezierCalls / double(r.width*r.height));
}

void rasterizeFont(const Font &font, const RasterizeOptions &options, uint32_t width, uint32_t height, float *distances)
{
	Rasterizer r = { font, options, width, height, distances };
	renderSdf<float>(r);
	// renderSdf<SseFloat4>(r);
}
