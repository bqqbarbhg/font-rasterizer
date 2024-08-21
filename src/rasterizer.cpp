#include "rasterizer.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>

// -- Implementation

uint64_t bezierCalls;
uint64_t bezierCallsNR;

template <typename T, typename U>
static bool intersectLine(T y, const Vec2<U> &A, const Vec2<U> &B, float* hits)
{
	T p0 = T(A.y) - y, p1 = T(B.y) - y;
	T a = p1 - p0;
	T b = p0;

	T ax = T(B.x - A.x);
	T bx = T(A.x);

	T t = -b / a;
	auto hit = (abs(a) >= 0.0001f) & (t >= 0.0f) & (t <= 1.0f);
	vectorStore(hits, maskSelect(hit, ax*t + bx, T(Inf)));
	return anyTrue(hit);
}

template <typename T, typename U>
static bool intersectLineEx2(T y, const Vec2<U> &A, const Vec2<U> &B, T xs[1])
{
	T p0 = T(A.y), p1 = T(B.y);
	T a = p1 - p0;
	T b = p0 - y;

	T ax = T(B.x - A.x);
	T bx = T(A.x);

	T t = -b / a;
	auto hit = (abs(a) >= 0.0001f) & (t >= 0.0f) & (t <= 1.0f);
	xs[0] = maskSelect(hit, ax*t + bx, T(Inf));
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

template <typename T>
struct BezierEx
{
	Vec2<T> a, b, c;

	void inPlaceSwapXY() {
		a = a.yx();
		b = b.yx();
		c = c.yx();
	}
};

template <typename T, typename U>
static BezierEx<T> setupBezierEx(const Vec2<U> &A, const Vec2<U> &B, const Vec2<U> &C)
{
	BezierEx<T> b;
	b.a = (A - B * 2.0f + C).cast<T>();
	b.b = ((B - A) * 2.0f).cast<T>();
	b.c = A.cast<T>();
	return b;
}

template <typename T>
static inline bool intersectBezierEx(T y, const BezierEx<T> &bezier, float *hits, float *ts)
{
	T a = bezier.a.y;
	T b = bezier.b.y;
	T c = bezier.c.y - y;

	T ax = bezier.a.x;
	T bx = bezier.b.x;
	T cx = bezier.c.x;

	bool anyHit = false;
	T miss = Inf;
	T hit0 = miss, hit1 = miss;
	T t0 = miss, t1 = miss;

	auto valid = abs(a) > 0.0f;
	if (anyTrue(valid)) {
		T disc2 = b*b - a*c*4.0f;

		auto exists = valid & (disc2 >= 0.0f);
		if (anyTrue(exists)) {
			T disc = sqrt(disc2);

			T rcpA = T(-0.5f) / a;
			t0 = (b + disc) * rcpA;
			t1 = (b - disc) * rcpA;

			auto h0 = exists & (t0 >= 0.0f) & (t0 <= 1.0f);
			auto h1 = exists & (t1 >= 0.0f) & (t1 <= 1.0f);

			hit0 = maskSelect(h0, (ax*t0+bx)*t0+cx, miss);
			hit1 = maskSelect(h1, (ax*t1+bx)*t1+cx, miss);
			anyHit = anyTrue(h0 | h1);
		}
	}
	if (anyFalse(valid)) {
		t0 = c / -b;
		auto h = !valid & (t0 >= 0.0f) & (t0 <= 1.0f);
		hit0 = maskSelect(h, (ax*t0+bx)*t0+cx, hit0);
		anyHit |= anyTrue(h);
	}

	if (anyHit) {
		vectorStore(hits + 0, hit0);
		vectorStore(hits + vectorWidth<T>(), hit1);
		vectorStore(ts + 0, t0);
		vectorStore(ts + vectorWidth<T>(), t1);
		return true;
	}
	return false;
}

template <typename T>
static inline bool intersectBezierEx2(T y, const BezierEx<T> &bezier, T xs[2], T ts[2], T windings[2])
{
	T a = bezier.a.y * 2.0f;
	T b = bezier.b.y;
	T c = bezier.c.y - y;

	T ax = bezier.a.x;
	T bx = bezier.b.x;
	T cx = bezier.c.x;

	bool anyHit = false;
	T miss = Inf;
	T hit0 = miss, hit1 = miss;
	T t0 = miss, t1 = miss;

	auto valid = abs(a) > 0.0f;
	if (anyTrue(valid)) {
		T disc2 = b*b - a*c*2.0f;

		auto exists = valid & (disc2 >= 0.0f);
		if (anyTrue(exists)) {
			T disc = sqrt(disc2);

			T rcpA = T(-1.0f) / a;
			t0 = (b + disc) * rcpA;
			t1 = (b - disc) * rcpA;

			auto rh0 = exists & (t0 >= 0.0f) & (t0 <= 1.0f);
			auto rh1 = exists & (t1 >= 0.0f) & (t1 <= 1.0f);

			// Compact hits
			auto h0 = rh0 | rh1;
			auto h1 = rh0 & rh1;
			t0 = maskSelect(rh0, t0, t1);

			t0 = maskSelect(h0, t0, Inf);
			t1 = maskSelect(h1, t1, Inf);

			hit0 = maskSelect(h0, (ax*t0+bx)*t0+cx, miss);
			hit1 = maskSelect(h1, (ax*t1+bx)*t1+cx, miss);
			anyHit = anyTrue(h0);
		}
	}
	if (anyFalse(valid)) {
		t0 = c / -b;
		auto h = !valid & (t0 >= 0.0f) & (t0 <= 1.0f);
		hit0 = maskSelect(h, (ax*t0+bx)*t0+cx, hit0);
		t0 = maskSelect(h, t0, Inf);
		anyHit |= anyTrue(h);
	}

	if (anyHit) {
		xs[0] = hit0;
		xs[1] = hit1;
		ts[0] = t0;
		ts[1] = t1;
		windings[0] = a*t0 + b;
		windings[1] = a*t1 + b;
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

struct BezierJob
{
	uint32_t bezierIndex;
	uint32_t px;
	uint32_t py;
};

struct Rasterizer
{
	const Font &font;
	const RasterizeOptions &options;
	uint32_t width;
	uint32_t height;
	float *distances;

	uint32_t safeWidth;
	uint32_t safeHeight;
	FontBvh bvh;

	std::vector<float> distSq;
	std::vector<float> offsetX;
	std::vector<float> offsetY;
	std::vector<int32_t> winding;

	std::vector<uint32_t> scratch;
	std::vector<BezierJob> bezierJobs;
	std::vector<float> rasterYs;
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
inline T inverseTransformX(Rasterizer &r, T x) {
	return (x - r.options.offset.x) / r.options.scale.x;
}

template <typename T>
inline T inverseTransformY(Rasterizer &r, T y) {
	return (y - r.options.offset.y) / r.options.scale.y;
}

template <typename T>
T nudgeOffInteger(T y) {
	const float epsilonUp =   1.000001f;
	const float epsilonDown = 0.999999f;
	T base = floor(y);
	return clamp(y, base * epsilonUp, (base + 1.0f) * epsilonDown);
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
void setupJobs(Rasterizer &r)
{
	r.scratch.resize(r.width * r.height);

	constexpr uint32_t step = vectorWidth<T>();
	float hits[step * 2];
	float ts[step * 2];

	uint32_t token = 0;

	int32_t padding = 2;

	uint32_t bezierIndex = 0;
	for (const Bezier &bezier : r.font.beziers) {
		Vec2f boundsMin = min(min(bezier.a, bezier.b), bezier.c);
		Vec2f boundsMax = max(max(bezier.a, bezier.b), bezier.c);

		float boundX = (float)(r.width - 1);
		float boundY = (float)(r.height - 1);

		uint32_t px0 = (uint32_t)clamp(inverseTransformX(r, boundsMin.x), 0.0f, boundX);
		uint32_t px1 = (uint32_t)clamp(inverseTransformX(r, boundsMax.x) + 1, 0.0f, boundX);
		uint32_t py0 = (uint32_t)clamp(inverseTransformY(r, boundsMax.y), 0.0f, boundY);
		uint32_t py1 = (uint32_t)clamp(inverseTransformY(r, boundsMin.y) + 1, 0.0f, boundY);

		token++;

		BezierEx<float> bex = setupBezierEx<float>(bezier.a, bezier.b, bezier.c);
		for (uint32_t py = py0; py < py1; py++) {
			float y = transformY(r, (float)py);
			if (!intersectBezierEx(y, bex, hits, ts)) continue;

			for (uint32_t i = 0; i < 2; i++) {
				if (hits[i] < Inf) {
					float x = hits[i];
					int32_t pxb = (int32_t)clamp(inverseTransformX(r, x), 0.0f, boundX);

					for (int32_t dx = -padding; dx <= padding; dx++) {
						int32_t pxi = pxb + dx;
						if (pxi < 0 || pxi >= (int32_t)r.width) continue;
						uint32_t px = (uint32_t)pxi;

						uint32_t &refToken = r.scratch[py * r.width + px];
						if (refToken == token) continue;
						refToken = token;
						r.bezierJobs.push_back({ bezierIndex, px, py });
					}
				}
			}
		}

		bex.inPlaceSwapXY();
		for (uint32_t px = px0; px < px1; px++) {
			float x = transformX(r, (float)px);
			if (!intersectBezierEx(x, bex, hits, ts)) continue;

			for (uint32_t i = 0; i < 2; i++) {
				if (hits[i] < Inf) {
					float y = hits[i];
					int32_t pyb = (int32_t)clamp(inverseTransformY(r, y), 0.0f, boundY);

					for (int32_t dy = -padding; dy <= padding; dy++) {
						int32_t pyi = pyb + dy;
						if (pyi < 0 || pyi >= (int32_t)r.height) continue;
						uint32_t py = (uint32_t)pyi;

						uint32_t &refToken = r.scratch[py * r.width + px];
						if (refToken == token) continue;
						refToken = token;
						r.bezierJobs.push_back({ bezierIndex, px, py });
					}
				}
			}
		}

		bezierIndex++;
	}

	while (r.bezierJobs.size() % step != 0) {
		BezierJob last = r.bezierJobs.back();
		r.bezierJobs.push_back(last);
	}
}


template <typename T>
void calculateDistanceJobs(Rasterizer &r)
{
	constexpr uint32_t step = vectorWidth<T>();

	float initX[step];
	float initY[step];
	float initAX[step];
	float initAY[step];
	float initBX[step];
	float initBY[step];
	float initCX[step];
	float initCY[step];
	float dist[step];

	uint32_t pixelCount = r.width * r.height;
	for (uint32_t i = 0; i < pixelCount; i++) {
		r.distances[i] = Inf;
	}

	uint32_t bezierJobCount = (uint32_t)r.bezierJobs.size();
	for (uint32_t i = 0; i < bezierJobCount; i += step) {
		for (uint32_t j = 0; j < step; j++) {
			BezierJob &job = r.bezierJobs[i + j];
			const Bezier &bezier = r.font.beziers[job.bezierIndex];
			initX[j] = (float)job.px;
			initY[j] = (float)job.py;
			initAX[j] = bezier.a.x;
			initAY[j] = bezier.a.y;
			initBX[j] = bezier.b.x;
			initBY[j] = bezier.b.y;
			initCX[j] = bezier.c.x;
			initCY[j] = bezier.c.y;
		}

		Vec2<T> pos = vec2(vectorLoad<T>(initX), vectorLoad<T>(initY));
		pos.x = transformX(r, pos.x);
		pos.y = transformY(r, pos.y);
		Vec2<T> a = vec2(vectorLoad<T>(initAX), vectorLoad<T>(initAY));
		Vec2<T> b = vec2(vectorLoad<T>(initBX), vectorLoad<T>(initBY));
		Vec2<T> c = vec2(vectorLoad<T>(initCX), vectorLoad<T>(initCY));

		T d = distanceToBezier(pos, a, b, c);
		
		vectorStore(dist, d);
		for (uint32_t j = 0; j < step; j++) {
			BezierJob &job = r.bezierJobs[i + j];
			float &dr = r.distances[job.py * r.width + job.px];
			dr = min(dr, dist[j]);
		}
	}

	for (uint32_t i = 0; i < pixelCount; i += step) {
		T t = vectorLoad<T>(r.distances + i);
		t = sqrt(t);
		vectorStore(r.distances + i, t);
	}
}

#define STATS 0

template <typename T>
struct BezierNR
{
	Vec2<T> A, B, C;
	Vec2<T> a, b;
	T aa, ab, bb;
};

template <typename T>
BezierNR<T> setupBezierNR(const Vec2<T> &A, const Vec2<T> &B, const Vec2<T> &C)
{
	BezierNR<T> bz;
	bz.A = A;
	bz.B = B;
	bz.C = C;
	bz.a = A - B*2.0f + C;
	bz.b = (B - A) * 2.0f;
	bz.aa = dot(bz.a, bz.a);
	bz.ab = dot(bz.a, bz.b);
	bz.bb = dot(bz.b, bz.b);
	return bz;
}

template <typename T>
T distanceToBezierNR(const Vec2<T> &pos, const Vec2<T> &A, const Vec2<T> &B, const Vec2<T> &C, T t)
{
	bezierCallsNR++;

	// TODO: Precomputed
	Vec2<T> a = A - B*2.0f + C;
	Vec2<T> b = (B - A) * 2.0f;
	Vec2<T> c = A - pos;

	T aa = dot(a, a);
	T ab = dot(a, b);
	T bb = dot(b, b);
	T ac = dot(a, c);
	T bc = dot(b, c);

	Vec2<T> wr = c + (b + a*t)*t;
	T dr = dot2(wr);

	T prevT = t;

	{
		T vy = bc + (bb+ac*2.0f)*t + (ab*3.0f)*t*t + aa*2.0f*t*t*t;
		T dy = (bb + ac*2.0f) + (ab*6.0f)*t + (aa*6.0f)*t*t;
		t = t - vy / dy;
	}

	{
		T vy = bc + (bb+ac*2.0f)*t + (ab*3.0f)*t*t + aa*2.0f*t*t*t;
		T dy = (bb + ac*2.0f) + (ab*6.0f)*t + (aa*6.0f)*t*t;
		t = t - vy / dy;
	}

	// printf("%f -> %f\n", prevT, t);

	if (abs(t - prevT) > 0.5f) {
		return distanceToBezier(pos, A, B, C);
	}
	t = saturate(t);

	// t = saturate(t - kx);
	Vec2<T> w = c + (b + a*t)*t;
	return min(dr, dot2(w));
}

template <typename T>
__forceinline T distanceToBezierNR(const Vec2<T> &pos, const BezierNR<T> &bz, T t)
{
	bezierCallsNR++;

	Vec2<T> c = bz.A - pos;
	T aa = bz.aa;
	T ab = bz.ab;
	T bb = bz.bb;
	T ac = dot(bz.a, c);
	T bc = dot(bz.b, c);

	T p = bb + ac * 2.0f;
	T q = ab * 3.0f;
	T r = aa * 2.0f;

	Vec2<T> wr = c + (bz.b + bz.a*t)*t;
	T dr = dot2(wr);

	T prevT = t;

	{
		T t2 = t*t;
		T t3 = t2*t;

		T vy = bc + p*t + q*t2 + r*t3;
		T dy = p + q*t*2.0f + r*3.0f*t2;
		t = t - vy / dy;
	}

#if 1
	{
		T t2 = t*t;
		T t3 = t2*t;

		T vy = bc + p*t + q*t2 + r*t3;
		T dy = p + q*t*2.0f + r*3.0f*t2;
		t = t - vy / dy;
	}
#endif

	// printf("%f -> %f\n", prevT, t);

	if (abs(t - prevT) > 0.5f) {
		return distanceToBezier(pos, bz.A, bz.B, bz.C);
	}
	t = saturate(t);

	// t = saturate(t - kx);
	Vec2<T> w = c + (bz.b + bz.a*t)*t;
	return min(dr, dot2(w));
}

template <typename T>
void calculateDistanceNR(Rasterizer &r)
{
	constexpr uint32_t step = vectorWidth<T>();
	float hits[step * 2];
	float ts[step * 2];

	int32_t padding = 2;

	uint32_t pixelCount = r.width * r.height;
	for (uint32_t i = 0; i < pixelCount; i++) {
		r.distances[i] = Inf;
	}

	float errorSum = 0.0f;
	float errorMax = 0.0f;
	float errorCount = 0.0f;

	uint32_t bezierIndex = 0;
	for (const Bezier &bezier : r.font.beziers) {
		Vec2f boundsMin = min(min(bezier.a, bezier.b), bezier.c);
		Vec2f boundsMax = max(max(bezier.a, bezier.b), bezier.c);

		float boundX = (float)(r.width - 1);
		float boundY = (float)(r.height - 1);

		uint32_t px0 = (uint32_t)clamp(inverseTransformX(r, boundsMin.x), 0.0f, boundX);
		uint32_t px1 = (uint32_t)clamp(inverseTransformX(r, boundsMax.x) + 1, 0.0f, boundX);
		uint32_t py0 = (uint32_t)clamp(inverseTransformY(r, boundsMax.y), 0.0f, boundY);
		uint32_t py1 = (uint32_t)clamp(inverseTransformY(r, boundsMin.y) + 1, 0.0f, boundY);

		BezierEx<float> bex = setupBezierEx<float>(bezier.a, bezier.b, bezier.c);
		BezierNR<float> bnr = setupBezierNR<float>(bezier.a, bezier.b, bezier.c);
		for (uint32_t py = py0; py < py1; py++) {
			float y = transformY(r, (float)py);
			if (!intersectBezierEx(y, bex, hits, ts)) continue;

			for (uint32_t i = 0; i < 2; i++) {
				if (hits[i] < Inf) {
					float x = hits[i];
					int32_t pxb = (int32_t)clamp(inverseTransformX(r, x), 0.0f, boundX);

					for (int32_t dx = -padding; dx <= padding; dx++) {
						int32_t pxi = pxb + dx;
						if (pxi < 0 || pxi >= (int32_t)r.width) continue;
						uint32_t px = (uint32_t)pxi;

						float x = transformX(r, (float)px);
						Vec2<T> pos = vec2(x, y);
						T d2 = distanceToBezierNR(pos, bnr, ts[i]);
						// T d2 = distanceToBezierNR(pos, bezier.a, bezier.b, bezier.c, ts[i]);
						bezierCalls--;
#if STATS
						T dref = distanceToBezier(pos, bezier.a, bezier.b, bezier.c);
						float error = abs(sqrt(dref) - sqrt(d2)) / r.options.scale.x;
						errorSum += error;
						errorMax = max(errorMax, error);
						errorCount += 1.0f;
#endif
						float &dr = r.distances[py * r.width + px];
						dr = min(dr, d2);
					}
				}
			}
		}

		bex.inPlaceSwapXY();
		for (uint32_t px = px0; px < px1; px++) {
			float x = transformX(r, (float)px);
			if (!intersectBezierEx(x, bex, hits, ts)) continue;

			for (uint32_t i = 0; i < 2; i++) {
				if (hits[i] < Inf) {
					float y = hits[i];
					int32_t pyb = (int32_t)clamp(inverseTransformY(r, y), 0.0f, boundY);

					for (int32_t dy = -padding; dy <= padding; dy++) {
						int32_t pyi = pyb + dy;
						if (pyi < 0 || pyi >= (int32_t)r.height) continue;
						uint32_t py = (uint32_t)pyi;

						float y = transformY(r, (float)py);
						Vec2<T> pos = vec2(x, y);
						T d2 = distanceToBezierNR(pos, bnr, ts[i]);
						// T d2 = distanceToBezierNR(pos, bezier.a, bezier.b, bezier.c, ts[i]);
#if STATS
						bezierCalls--;
						T dref = distanceToBezier(pos, bezier.a, bezier.b, bezier.c);
						float error = abs(sqrt(dref) - sqrt(d2)) / r.options.scale.x;
						errorSum += error;
						errorMax = max(errorMax, error);
						errorCount += 1.0f;
#endif
						float &dr = r.distances[py * r.width + px];
						dr = min(dr, d2);
					}
				}
			}
		}

		bezierIndex++;
	}

#if STATS
	printf("error: %f (max %f)\n", errorSum / errorCount, errorMax);
#endif

	for (uint32_t i = 0; i < pixelCount; i += step) {
		T t = vectorLoad<T>(r.distances + i);
		t = sqrt(t);
		vectorStore(r.distances + i, t);
	}
}

#define VERBOSE 0

#if VERBOSE
	#define verbosef(...) printf(__VA_ARGS__)
#else
	#define verbosef(...)
#endif

template <typename T>
__declspec(noinline) void linePassX(Rasterizer &r, const Line &line, uint32_t py, T y, T x)
{
	constexpr uint32_t step = vectorWidth<T>();

	uint32_t width = r.width;

	int32_t windingDelta = line.b.y - line.a.y > 0.0f ? +1 : -1;

	float pxs[step];
	vectorStore(pxs, clamp(inverseTransformX(r, x), T(0.0f), T((float)r.width)));

	// Accumulate winding
	for (uint32_t i = 0; i < step; i++) {
		uint32_t px = (uint32_t)pxs[i];
		if (px < width) {
			verbosef("line %u,%u : %+d\n", px, py+i, windingDelta);
			r.winding[(py + i) * r.safeWidth + px] += windingDelta;
		}
	}
}

static int serial = 0;

template <typename T>
__declspec(noinline) void bezierPassX(Rasterizer &r, uint32_t py, T y, T x, T t, T winding)
{
	constexpr uint32_t step = vectorWidth<T>();

	uint32_t width = r.width;

	uint32_t windingBits = vectorBitmask(winding >= 0.0f);

	float pxs[step];
	vectorStore(pxs, clamp(inverseTransformX(r, x), T(0.0f), T((float)r.width)));

	// Accumulate winding
	for (uint32_t i = 0; i < step; i++) {
		uint32_t px = (uint32_t)pxs[i];
		if (px < width) {
			int32_t windingDelta = (windingBits & (1 << i)) != 0 ? +1 : -1;
			verbosef("bezier %u,%u: %+d\n", px, py+i, windingDelta);
			r.winding[(py + i) * r.safeWidth + px] += windingDelta;
		}
	}
}

void expandWinding(Rasterizer &r)
{
	uint32_t width = r.width;

	for (uint32_t py = 0; py < r.height; py++) {
		int32_t *windings = r.winding.data() + py * r.safeWidth;
		float *dists = r.distances + py * r.width;

		int32_t winding = 0;
		for (uint32_t px = 0; px < width; px++) {
			winding += windings[px];
			windings[px] = winding;

			dists[px] = winding != 0 ? 1.0f : 0.0f;
		}
	}
}

template <typename T>
void renderSdfNew(Rasterizer &r)
{
	constexpr uint32_t step = vectorWidth<T>();

	r.safeHeight = r.height + step * 2;
	r.safeWidth = r.width + step * 2;
	uint32_t safeSize = r.safeWidth * r.safeHeight;

	r.distSq.resize(safeSize, Inf);
	r.offsetX.resize(safeSize);
	r.offsetY.resize(safeSize);
	r.winding.resize(safeSize);
	r.rasterYs.resize(r.safeHeight);

	float maximumX = (float)(r.width - 1);
	float maximumY = (float)(r.height - 1);

	uint32_t height = r.height;

	float *rasterYs = r.rasterYs.data();

	// Precompute Y values
	for (uint32_t py = 0; py < height; py += step) {
		T y = T(float(py)) + vectorSequence<T>();
		y = transformY(r, y);
		y = nudgeOffInteger(y);
		vectorStore(rasterYs + py, y);
	}

	// -- Lines
	int lineIndex = 0;
	for (const Line &line : r.font.lines) {
		verbosef("== line %d\n", lineIndex++);
		Vec2f boundsMin = min(line.a, line.b);
		Vec2f boundsMax = max(line.a, line.b);

		uint32_t px0 = (uint32_t)clamp(inverseTransformX(r, boundsMin.x), 0.0f, maximumX);
		uint32_t px1 = (uint32_t)clamp(inverseTransformX(r, boundsMax.x) + 1, 0.0f, maximumX);
		uint32_t py0 = (uint32_t)clamp(inverseTransformY(r, boundsMax.y), 0.0f, maximumY);
		uint32_t py1 = (uint32_t)clamp(inverseTransformY(r, boundsMin.y) + 1, 0.0f, maximumY);

		for (uint32_t py = py0; py < py1; py += step) {
			T y = vectorLoad<T>(rasterYs + py);

			T xs[1];
			if (!intersectLineEx2(y, line.a, line.b, xs)) continue;

			linePassX<T>(r, line, py, y, xs[0]);
		}

	}

	// -- Beziers
	int bezierIndex = 0;
	for (const Bezier &bezier : r.font.beziers) {
		verbosef("== bezier %d\n", bezierIndex++);

		Vec2f boundsMin = min(min(bezier.a, bezier.b), bezier.c);
		Vec2f boundsMax = max(max(bezier.a, bezier.b), bezier.c);

		uint32_t px0 = (uint32_t)clamp(inverseTransformX(r, boundsMin.x), 0.0f, maximumX);
		uint32_t px1 = (uint32_t)clamp(inverseTransformX(r, boundsMax.x) + 1, 0.0f, maximumX);
		uint32_t py0 = (uint32_t)clamp(inverseTransformY(r, boundsMax.y), 0.0f, maximumY);
		uint32_t py1 = (uint32_t)clamp(inverseTransformY(r, boundsMin.y) + 1, 0.0f, maximumY);

		BezierEx<T> bex = setupBezierEx<T>(bezier.a, bezier.b, bezier.c);
		for (uint32_t py = py0; py < py1; py += step) {
			T y = vectorLoad<T>(rasterYs + py);

			T xs[2], ts[2], windings[2];
			if (!intersectBezierEx2<T>(y, bex, xs, ts, windings)) continue;

			bezierPassX<T>(r, py, y, xs[0], ts[0], windings[0]);
			if (anyTrue(xs[1] < Inf)) {
				bezierPassX<T>(r, py, y, xs[1], ts[1], windings[1]);
			}
		}
	}

	expandWinding(r);
}

template <typename T>
void renderSdf(Rasterizer &r)
{
	bezierCalls = 0;
	bezierCallsNR = 0;

	// buildBvh(r.bvh, r.font);
	// setupWinding<T>(r);
	// Precomputed<T> pre;
	// precompute(pre, r.font);
	// calculateDistance<T>(r);
	// calculateDistanceBvh<T>(r, pre);

	// setupJobs<T>(r);
	// calculateDistanceJobs<T>(r);
	// calculateDistanceNR<T>(r);
	renderSdfNew<T>(r);

	//calculateDistancePre(r, pre);

#if STATS
	printf(">>> Bezier %llu (%.2f/px)\n", bezierCalls, (double)bezierCalls / double(r.width*r.height));
	printf(">>> Bezier NR %llu (%.2f/px)\n", bezierCallsNR, (double)bezierCallsNR / double(r.width*r.height));
#endif
}

void rasterizeFont(const Font &font, const RasterizeOptions &options, uint32_t width, uint32_t height, float *distances)
{
	Rasterizer r = { font, options, width, height, distances };
	renderSdf<float>(r);
	// renderSdf<SseFloat4>(r);
}
