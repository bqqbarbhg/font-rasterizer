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

template <typename T>
struct LineEx
{
	Vec2<T> a;
	Vec2<T> ba;
	Vec2<T> scaledBa;
	float rcpBay;
	bool vertical;
	uint32_t winding;
};

template <typename T>
LineEx<T> setupLineEx(const Line &line)
{
	Vec2f ba = line.b - line.a;
	LineEx<T> lex;
	lex.a = line.a.cast<T>();
	lex.ba = ba.cast<T>();
	lex.scaledBa = (ba / dot2(ba)).cast<T>();
	lex.rcpBay = 1.0f / ba.y;
	lex.vertical = abs(ba.y) >= abs(ba.x);
	lex.winding = ba.y >= 0.0f ? ~0u : 0;
	return lex;
}

template <typename T>
static bool intersectLineEx2(T y, const LineEx<T> &lex, T xs[1])
{
	T a = lex.ba.y;
	T b = T(lex.a.y) - y;

	T ax = lex.ba.x;
	T bx = lex.a.x;

	T t = -b * lex.rcpBay;
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

template <typename T>
struct CombinedArray {
	T *data;
	size_t size;

	const T &operator[](size_t ix) const {
		assert(ix < size);
		return data[ix];
	}
	T &operator[](size_t ix) {
		assert(ix < size);
		return data[ix];
	}
};

template <int MaxArrays>
struct CombinedAllocation {
	void *data = nullptr;

	struct ArrayReference
	{
		void **pointer;
		size_t offset;
	};

	ArrayReference arrays[MaxArrays];
	size_t arrayCount = 0;
	size_t totalSize = 0;
	bool doFree = false;

	template <typename T>
	void add(CombinedArray<T> &base, size_t size) {
		size_t align = alignof(T);
		totalSize += (totalSize + (align - 1)) & ~align;

		ArrayReference &ref = arrays[arrayCount++];
		ref.pointer = (void**)&base;
		ref.offset = totalSize;
		base.size = size;

		totalSize += sizeof(T) * size;
	}

	void impUpdateOffsets() {
		char *p = (char*)data;
		for (size_t i = 0; i < arrayCount; i++) {
			*arrays[i].pointer = p + arrays[i].offset;
		}
	}

	void allocate() {
		data = malloc(totalSize);
		doFree = true;
		impUpdateOffsets();
	}

	void setPointer(void *pointer, size_t size) {
		assert(size >= totalSize);
		data = pointer;
		impUpdateOffsets();
	}

	void setPointerOrAllocate(void *pointer, size_t size) {
		if (size >= totalSize) {
			setPointer(pointer, size);
		} else {
			allocate();
		}
	}

	~CombinedAllocation() {
		if (doFree) {
			free(data);
		}
	}
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

	CombinedArray<float> distSqData[2];
	CombinedArray<float> offsetXData[2];
	CombinedArray<float> offsetYData[2];
	CombinedArray<float> rasterXs;
	CombinedArray<float> rasterYs;
	CombinedArray<int32_t> winding;

	float *distSq;
	float *offsetY;
	float *offsetX;

	float *distSqDst;
	float *offsetYDst;
	float *offsetXDst;

	std::vector<uint32_t> scratch;
	std::vector<BezierJob> bezierJobs;
	int32_t padding;
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
__forceinline T distanceToBezierExNR(const Vec2<T> &pos, const BezierNR<T> &bz, T t, Vec2<T> &closest)
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

#if 0
	if (abs(t - prevT) > 0.5f) {
		return distanceToBezier(pos, bz.A, bz.B, bz.C);
	}
#endif

	t = saturate(t);

	// t = saturate(t - kx);
	Vec2<T> w = c + (bz.b + bz.a*t)*t;
	closest = w + pos;
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

struct PixelBounds
{
	uint32_t px0, px1;
	uint32_t py0, py1;
};

PixelBounds setupBounds(Rasterizer &r, const Vec2f &boundsMin, const Vec2f &boundsMax)
{
	float maximumX = (float)(r.width - 1);
	float maximumY = (float)(r.height - 1);
	uint32_t px0 = (uint32_t)clamp(inverseTransformX(r, boundsMin.x), 0.0f, maximumX);
	uint32_t px1 = (uint32_t)clamp(inverseTransformX(r, boundsMax.x) + 1.0f, 0.0f, maximumX);
	uint32_t py0 = (uint32_t)clamp(inverseTransformY(r, boundsMax.y), 0.0f, maximumY);
	uint32_t py1 = (uint32_t)clamp(inverseTransformY(r, boundsMin.y) + 1.0f, 0.0f, maximumY);
	return { px0, px1, py0, py1 };
}

PixelBounds setupBounds(Rasterizer &r, const Line &line)
{
	Vec2f boundsMin = min(line.a, line.b);
	Vec2f boundsMax = max(line.a, line.b);
	return setupBounds(r, boundsMin, boundsMax);
}

PixelBounds setupBounds(Rasterizer &r, const Bezier &bezier)
{
	Vec2f boundsMin = min(min(bezier.a, bezier.b), bezier.c);
	Vec2f boundsMax = max(max(bezier.a, bezier.b), bezier.c);
	return setupBounds(r, boundsMin, boundsMax);
}

template <typename T>
__forceinline T distanceToLineEx(const Vec2<T> &pos, const LineEx<T> &lex, Vec2<T> &closest)
{
    Vec2<T> pa = pos - lex.a;
    T t = saturate(dot(lex.scaledBa, pa));
	Vec2<T> p = lex.ba*t - pa;
	closest = pos + p;
    return dot2(p);
}

#define VALIDATE 0

__forceinline void storeClosestDistance(Rasterizer &r, uint32_t px, uint32_t py, float distSq, float x, float y, bool windingBit)
{
#if VALIDATE
	{
		Vec2f pixel = vec2(transformX(r, (float)px), transformY(r, (float)py));
		Vec2f closest = vec2(x, y);
		float distance = sqrt(dot2(closest - pixel));
		float refDistance = sqrt(distSq);
		assert(abs(distance - refDistance) <= 1.0f);
	}
#endif

	uint32_t ix = py * r.safeWidth + px;
	if (distSq > r.distSq[ix]) return;

	int32_t winding = r.winding[ix];
	if (winding < -1 || winding > 1) return;
	winding = windingBit ? winding : -winding;
	if (winding >= 1) return;

	r.distSq[ix] = distSq;
	r.offsetX[ix] = x;
	r.offsetY[ix] = y;
}

template <typename T>
__forceinline void storeClosestDistance(Rasterizer &r, uint32_t px0[4], uint32_t py0, uint32_t dx, const T& distSq, const Vec2<T> &closest, uint32_t windingBits)
{
	constexpr uint32_t step = vectorWidth<T>();

	float distSqs[step];
	float posX[step];
	float posY[step];
	vectorStore(distSqs, distSq);
	vectorStore(posX, closest.x);
	vectorStore(posY, closest.y);

	for (uint32_t i = 0; i < step; i++) {
		uint32_t py = py0 + i;
		uint32_t px = px0[i] + dx;
		if (px < r.width) {
			storeClosestDistance(r, px, py, distSqs[i], posX[i], posY[i], (windingBits & (1 << i)) != 0);
		}
	}
}

template <typename T>
__forceinline void storeClosestDistance(Rasterizer &r, uint32_t px0, uint32_t py0[4], uint32_t dy, const T& distSq, const Vec2<T> &closest, uint32_t windingBits)
{
	constexpr uint32_t step = vectorWidth<T>();

	float distSqs[step];
	float posX[step];
	float posY[step];
	vectorStore(distSqs, distSq);
	vectorStore(posX, closest.x);
	vectorStore(posY, closest.y);

	for (uint32_t i = 0; i < step; i++) {
		uint32_t px = px0 + i;
		uint32_t py = py0[i] + dy;
		if (py < r.height) {
			storeClosestDistance(r, px, py, distSqs[i], posX[i], posY[i], (windingBits & (1 << i)) != 0);
		}
	}
}

constexpr float MaxPixelY = 1048576.0f;

template <typename T>
__declspec(noinline) void linePassX(Rasterizer &r, const LineEx<T> &lex, uint32_t py0, T y, T x, uint32_t windingBits)
{
	constexpr uint32_t step = vectorWidth<T>();

	T pxRaw = inverseTransformX(r, x);
	T pxFloat = rint(clamp(pxRaw, T(0.0f), T(MaxPixelY)));

	uint32_t px0[step];
	vectorStoreUint(px0, pxFloat);

	int32_t padding = r.padding;
	for (int32_t dx = -padding; dx <= padding; dx++) {
		T pxShifted = pxFloat + float(dx);

		uint32_t winding = vectorBitmask(pxShifted > pxRaw) ^ windingBits;

		T x = transformX(r, pxShifted);
		Vec2<T> closest, pos = vec2(x, y);
		T distSq = distanceToLineEx(pos, lex, closest);
		storeClosestDistance(r, px0, py0, dx, distSq, closest, winding);
	}
}

template <typename T>
__declspec(noinline) void linePassY(Rasterizer &r, const LineEx<T> &lex, uint32_t px0, T x, T y, uint32_t windingBits)
{
	constexpr uint32_t step = vectorWidth<T>();

	T pyRaw = maskSelect(y < Inf, inverseTransformY(r, y), T(Inf));
	T pyFloat = rint(clamp(pyRaw, T(0.0f), T(MaxPixelY)));

	uint32_t py0[step];
	vectorStoreUint(py0, pyFloat);

	int32_t padding = r.padding;
	for (int32_t dy = -padding; dy <= padding; dy++) {
		T pyShifted = pyFloat + float(dy);

		uint32_t winding = vectorBitmask(pyShifted > pyRaw) ^ windingBits;

		T y = transformY(r, pyShifted);
		Vec2<T> closest, pos = vec2(y, x);;
		T distSq = distanceToLineEx(pos, lex, closest);
		storeClosestDistance(r, px0, py0, dy, distSq, closest.yx(), winding);
	}
}

template <typename T>
__declspec(noinline) void bezierPassX(Rasterizer &r, BezierNR<T> &bnr, uint32_t py0, T y, T x, T t, T winding)
{
	constexpr uint32_t step = vectorWidth<T>();

	T pxRaw = inverseTransformX(r, x);
	T pxFloat = rint(clamp(pxRaw, T(0.0f), T(MaxPixelY)));

	uint32_t px0[step];
	vectorStoreUint(px0, pxFloat);

	int32_t padding = r.padding;
	for (int32_t dx = -padding; dx <= padding; dx++) {
		T pxShifted = pxFloat + float(dx);

		uint32_t windingBits = vectorBitmask((pxShifted > pxRaw) ^ (winding >= 0.0f));

		T x = transformX(r, pxShifted);
		Vec2<T> closest, pos = vec2(x, y);
		T distSq = distanceToBezierExNR(pos, bnr, t, closest);
		storeClosestDistance(r, px0, py0, dx, distSq, closest, windingBits);
	}
}

template <typename T>
__declspec(noinline) void bezierPassY(Rasterizer &r, BezierNR<T> &bnr, uint32_t px0, T x, T y, T t, T winding)
{
	constexpr uint32_t step = vectorWidth<T>();

	T pyRaw = maskSelect(y < Inf, inverseTransformY(r, y), T(Inf));
	T pyFloat = floor(clamp(pyRaw, T(0.0f), T(MaxPixelY)));

	uint32_t py0[step];
	vectorStoreUint(py0, pyFloat);

	int32_t padding = r.padding;
	for (int32_t dy = -padding; dy <= padding; dy++) {
		T pyShifted = pyFloat + float(dy);

		uint32_t windingBits = vectorBitmask((pyShifted > pyRaw) ^ (winding >= 0.0f));

		T y = transformY(r, pyShifted);
		Vec2<T> closest, pos = vec2(x, y);
		T distSq = distanceToBezierExNR(pos, bnr, t, closest);
		storeClosestDistance(r, px0, py0, dy, distSq, closest, windingBits);
		storeClosestDistance(r, px0, py0, dy, distSq, closest, windingBits);
	}
}

void vertexPass(Rasterizer &r, const Vec2f &point)
{
	float pxf = max(rint(inverseTransformX(r, point.x)), 0.0f);
	float pyf = max(rint(inverseTransformY(r, point.y)), 0.0f);
	uint32_t px0 = (uint32_t)pxf;
	uint32_t py0 = (uint32_t)pyf;
	int32_t padding = r.padding;

	// TODO: SIMD this!

	// Check if the vertex is in a controversial area
	bool controversial = false;
	for (int32_t dy = -padding; dy <= padding; dy++) {
		uint32_t py = py0 + dy;
		if (py >= r.height) continue;

		for (int32_t dx = -padding; dx <= padding; dx++) {
			uint32_t px = px0 + dx;
			if (px >= r.width) continue;

			uint32_t ix = py * r.safeWidth + px;
			uint32_t winding = r.winding[ix];
			if (winding < 1 || winding > 1) {
				controversial = true;
			}
		}
	}

	for (int32_t dy = -padding; dy <= padding; dy++) {
		uint32_t py = py0 + dy;
		if (py >= r.height) continue;
		float refY = r.rasterYs[py];

		for (int32_t dx = -padding; dx <= padding; dx++) {
			uint32_t px = px0 + dx;
			if (px >= r.width) continue;
			float refX = r.rasterXs[px];

			uint32_t ix = py * r.safeWidth + px;

			uint32_t winding = r.winding[ix];
			if (controversial && winding != 0) continue;

			Vec2f pf = vec2(refX, refY);
			float d = dot2(pf - point);

			if (d < r.distSq[ix]) {
				r.distSq[ix] = d;
				r.offsetX[ix] = point.x;
				r.offsetY[ix] = point.y;
			}
		}
	}
}

__declspec(noinline) void expandWinding(Rasterizer &r)
{
	uint32_t width = r.width;

#if 0
	for (uint32_t py = 0; py < r.height; py += 4) {
		int32_t *p0 = r.winding.data + (py + 0) * r.safeWidth;
		int32_t *p1 = r.winding.data + (py + 1) * r.safeWidth;
		int32_t *p2 = r.winding.data + (py + 2) * r.safeWidth;
		int32_t *p3 = r.winding.data + (py + 3) * r.safeWidth;

		int32_t w0 = 0;
		int32_t w1 = 0;
		int32_t w2 = 0;
		int32_t w3 = 0;

		for (uint32_t px = 0; px < width; px++) {
			w0 += p0[px]; p0[px] = w0;
			w1 += p1[px]; p1[px] = w1;
			w2 += p2[px]; p2[px] = w2;
			w3 += p3[px]; p3[px] = w3;
		}
	}
#else
	constexpr int ma = (int)0x03020100;
	constexpr int mb = (int)0x07060504;
	constexpr int mc = (int)0x0b0a0908;

	const __m128i mask_0aba = _mm_setr_epi32(-1, ma, mb, ma);
	const __m128i mask_00ac = _mm_setr_epi32(-1,-1, ma, mc);

	for (uint32_t py = 0; py < r.height; py += 2) {
		int32_t *ptr0 = r.winding.data + (py + 0) * r.safeWidth;
		int32_t *ptr1 = r.winding.data + (py + 1) * r.safeWidth;

		__m128i prefix0 = _mm_setzero_si128();
		__m128i prefix1 = _mm_setzero_si128();
		for (uint32_t px = 0; px < width; px += 4) {
			__m128i v0 = _mm_loadu_si128((const __m128i*)ptr0);
			__m128i v1 = _mm_loadu_si128((const __m128i*)ptr1);

			v0 = _mm_add_epi32(v0, _mm_shuffle_epi8(v0, mask_0aba));
			v0 = _mm_add_epi32(v0, _mm_shuffle_epi8(v0, mask_00ac));
			v0 = _mm_add_epi32(v0, prefix0);
			prefix0 = _mm_shuffle_epi32(v0, _MM_SHUFFLE(3,3,3,3));

			v1 = _mm_add_epi32(v1, _mm_shuffle_epi8(v1, mask_0aba));
			v1 = _mm_add_epi32(v1, _mm_shuffle_epi8(v1, mask_00ac));
			v1 = _mm_add_epi32(v1, prefix1);
			prefix1 = _mm_shuffle_epi32(v1, _MM_SHUFFLE(3,3,3,3));

			_mm_storeu_si128((__m128i*)ptr0, v0);
			_mm_storeu_si128((__m128i*)ptr1, v1);

			ptr0 += 4;
			ptr1 += 4;
		}
	}
#endif
}

void scanFill(Rasterizer &r)
{
	int32_t height = r.height;
	int32_t width = r.width;
	int32_t safeWidth = int32_t(r.safeWidth);

	float *srcD = r.distSq;
	float *srcX = r.offsetX;
	float *srcY = r.offsetY;

	for (int32_t py = 0; py < height; py++) {
		for (int32_t px = 0; px < width; px++) {
			int32_t ix = py * safeWidth + px;

			float selfX = r.rasterXs[px];
			float selfY = r.rasterYs[py];

			float offX = srcX[ix];
			float offY = srcY[ix];
			float selfD = srcD[ix];
			bool changed = false;

			int32_t nb = ix - 1;
			if (srcD[nb] < selfD) {
				float nbX = srcX[nb];
				float nbY = srcY[nb];
				float dx = nbX - selfX;
				float dy = nbY - selfY;
				float nbD = max(srcD[nb], dx*dx + dy*dy);
				if (nbD < selfD) {
					selfD = nbD;
					offX = nbX;
					offY = nbY;
					changed = true;
				}
			}

			nb = ix - safeWidth;
			if (srcD[nb] < selfD) {
				float nbX = srcX[nb];
				float nbY = srcY[nb];
				float dx = nbX - selfX;
				float dy = nbY - selfY;
				float nbD = max(srcD[nb], dx*dx + dy*dy);
				if (nbD < selfD) {
					selfD = nbD;
					offX = nbX;
					offY = nbY;
					changed = true;
				}
			}

			if (changed) {
				srcD[ix] = selfD;
				srcX[ix] = offX;
				srcY[ix] = offY;
			}
		}
	}

	for (int32_t py = height - 1; py >= 0; py--) {
		for (int32_t px = width - 1; px >= 0; px--) {
			int32_t ix = py * safeWidth + px;

			float selfX = r.rasterXs[px];
			float selfY = r.rasterYs[py];

			float offX = srcX[ix];
			float offY = srcY[ix];
			float selfD = srcD[ix];
			bool changed = false;

			int32_t nb = ix + 1;
			if (srcD[nb] < selfD) {
				float nbX = srcX[nb];
				float nbY = srcY[nb];
				float dx = nbX - selfX;
				float dy = nbY - selfY;
				float nbD = max(srcD[nb], dx*dx + dy*dy);
				if (nbD < selfD) {
					selfD = nbD;
					offX = nbX;
					offY = nbY;
					changed = true;
				}
			}

			nb = ix + safeWidth;
			if (srcD[nb] < selfD) {
				float nbX = srcX[nb];
				float nbY = srcY[nb];
				float dx = nbX - selfX;
				float dy = nbY - selfY;
				float nbD = max(srcD[nb], dx*dx + dy*dy);
				if (nbD < selfD) {
					selfD = nbD;
					offX = nbX;
					offY = nbY;
					changed = true;
				}
			}

			if (changed) {
				srcD[ix] = selfD;
				srcX[ix] = offX;
				srcY[ix] = offY;
			}
		}
	}
}

void scanFillConservative(Rasterizer &r)
{
	int32_t height = r.height;
	int32_t width = r.width;
	int32_t safeWidth = int32_t(r.safeWidth);

	float *srcD = r.distSq;
	float *srcX = r.offsetX;
	float *srcY = r.offsetY;

	float pixelScale = r.options.scale.x;

	for (int32_t py = 0; py < height; py++) {
		for (int32_t px = 0; px < width; px++) {
			int32_t ix = py * safeWidth + px;

			float selfX = r.rasterXs[px];
			float selfY = r.rasterYs[py];

			float selfD = srcD[ix];
			bool changed = false;

			int32_t nb = ix - 1;
			if (srcD[nb] < selfD) {
				float nbD = sqrt(srcD[nb]) + pixelScale;
				nbD *= nbD;
				if (nbD < selfD) {
					selfD = nbD;
					changed = true;
				}
			}

			nb = ix - safeWidth;
			if (srcD[nb] < selfD) {
				float nbD = sqrt(srcD[nb]) + pixelScale * sqrt(2.0f);
				nbD *= nbD;
				if (nbD < selfD) {
					selfD = nbD;
					changed = true;
				}
			}

			if (changed) {
				srcD[ix] = selfD;
			}
		}
	}

	for (int32_t py = height - 1; py >= 0; py--) {
		for (int32_t px = width - 1; px >= 0; px--) {
			int32_t ix = py * safeWidth + px;

			float selfX = r.rasterXs[px];
			float selfY = r.rasterYs[py];

			float offX = srcX[ix];
			float offY = srcY[ix];
			float selfD = srcD[ix];
			bool changed = false;

			int32_t nb = ix + 1;
			if (srcD[nb] < selfD) {
				float nbD = sqrt(srcD[nb]) + pixelScale;
				nbD *= nbD;
				if (nbD < selfD) {
					selfD = nbD;
					changed = true;
				}
			}

			nb = ix + safeWidth;
			if (srcD[nb] < selfD) {
				float nbD = sqrt(srcD[nb]) + pixelScale * sqrt(2.0f);
				nbD *= nbD;
				if (nbD < selfD) {
					selfD = nbD;
					changed = true;
				}
			}

			if (changed) {
				srcD[ix] = selfD;
			}
		}
	}
}

template <typename T>
void finalizeDistance(Rasterizer &r)
{
	constexpr uint32_t step = vectorWidth<T>();

	uint32_t width = r.width;

	for (uint32_t py = 0; py < r.height; py++) {
		float *dists = r.distances + py * r.width;
		float *distSq = r.distSq + py * r.safeWidth;
		int32_t *windings = r.winding.data + py * r.safeWidth;

		int32_t winding = 0;
		uint32_t px = 0;
		for (; px + step <= width; px += step) {
			T d = vectorLoad<T>(dists + px);
			T dsq = vectorLoad<T>(distSq + px);
			T result = flipSignIfNonZero(d * sqrt(dsq), windings + px);

			vectorStore(dists + px, result);
		}
		for (; px < width; px++) {
			float d = dists[px];
			float dsq = distSq[px];
			d = d * sqrt(dsq);
			d = windings[px] != 0 ? -d : d;
			dists[px] = d;
		}
	}
}

template <typename T>
__forceinline void addWinding(Rasterizer &r, uint32_t py0, T x, uint32_t windingBits)
{
	constexpr uint32_t step = vectorWidth<T>();

	float pxs[step];
	T hitPx = clamp(inverseTransformX(r, x) + 1.0f, T(0.0f), T((float)r.width));
	vectorStore(pxs, hitPx);

	// Accumulate winding
	for (uint32_t i = 0; i < step; i++) {
		uint32_t py = py0 + i;
		uint32_t px = (uint32_t)pxs[i];
		if (px < r.width) {
			int32_t windingDelta = (windingBits & (1 << i)) != 0 ? +1 : -1;
			r.winding[py * r.safeWidth + px] += windingDelta;
		}
	}
}

template <typename T>
__declspec(noinline) void computeContourWinding(Rasterizer &r)
{
	constexpr uint32_t step = vectorWidth<T>();

	float *rasterXs = r.rasterXs.data;
	float *rasterYs = r.rasterYs.data;

	for (const Line &line : r.font.lines) {
		PixelBounds bounds = setupBounds(r, line);
		LineEx<T> lex = setupLineEx<T>(line);
		for (uint32_t py = bounds.py0; py < bounds.py1; py += step) {
			T y = vectorLoad<T>(rasterYs + py);

			T xs[1];
			if (!intersectLineEx2(y, lex, xs)) continue;

			addWinding<T>(r, py, xs[0], lex.winding);
		}
	}

	for (const Bezier &bezier : r.font.beziers) {
		PixelBounds bounds = setupBounds(r, bezier);
		BezierEx<T> bex = setupBezierEx<T>(bezier.a, bezier.b, bezier.c);
		for (uint32_t py = bounds.py0; py < bounds.py1; py += step) {
			T y = vectorLoad<T>(rasterYs + py);

			T xs[2], ts[2], windings[2];
			if (!intersectBezierEx2<T>(y, bex, xs, ts, windings)) continue;

			addWinding<T>(r, py, xs[0], vectorBitmask(windings[0] >= 0.0f));
			if (anyTrue(xs[1] < Inf)) {
				addWinding<T>(r, py, xs[1], vectorBitmask(windings[1] >= 0.0f));
			}
		}
	}
}

template <typename T>
__declspec(noinline) void computeContourDistance(Rasterizer &r)
{
	constexpr uint32_t step = vectorWidth<T>();

	for (const Line &line : r.font.lines) {
		PixelBounds bounds = setupBounds(r, line);

		Vec2f delta = line.b - line.a;
		bool vertical = abs(delta.y) >= abs(delta.x);

		if (vertical) {
			LineEx<T> lex = setupLineEx<T>(line);
			for (uint32_t py = bounds.py0; py < bounds.py1; py += step) {
				T y = vectorLoad<T>(r.rasterYs.data + py);

				T xs[1];
				if (!intersectLineEx2(y, lex, xs)) continue;

				linePassX<T>(r, lex, py, y, xs[0], lex.winding);
			}
		} else {
			Line lineYX = { line.a.yx(), line.b.yx() };
			LineEx<T> lex = setupLineEx<T>(lineYX);
			for (uint32_t px = bounds.px0; px < bounds.px1; px += step) {
				T x = vectorLoad<T>(r.rasterXs.data + px);

				T ys[1];
				if (!intersectLineEx2(x, lex, ys)) continue;

				linePassY<T>(r, lex, px, x, ys[0], lex.winding);
			}
		}
	}

	for (const Bezier &bezier : r.font.beziers) {
		PixelBounds bounds = setupBounds(r, bezier);

		BezierNR<T> bnr = setupBezierNR<T>(bezier.a.cast<T>(), bezier.b.cast<T>(), bezier.c.cast<T>());
		BezierEx<T> bex = setupBezierEx<T>(bezier.a, bezier.b, bezier.c);

		for (uint32_t py = bounds.py0; py < bounds.py1; py += step) {
			T y = vectorLoad<T>(r.rasterYs.data + py);

			T xs[2], ts[2], windings[2];
			if (!intersectBezierEx2<T>(y, bex, xs, ts, windings)) continue;

			bezierPassX<T>(r, bnr, py, y, xs[0], ts[0], windings[0]);
			if (anyTrue(xs[1] < Inf)) {
				bezierPassX<T>(r, bnr, py, y, xs[1], ts[1], windings[1]);
			}
		}

		bex.inPlaceSwapXY();

		for (uint32_t px = bounds.px0; px < bounds.px1; px += step) {
			T x = vectorLoad<T>(r.rasterXs.data + px);

			T ys[2], ts[2], windings[2];
			if (!intersectBezierEx2<T>(x, bex, ys, ts, windings)) continue;

			bezierPassY<T>(r, bnr, px, x, ys[0], ts[0], windings[0]);
			if (anyTrue(ys[1] < Inf)) {
				bezierPassY<T>(r, bnr, px, x, ys[1], ts[1], windings[1]);
			}
		}
	}
}

template <typename T>
__declspec(noinline) void clearInf(float *dst, size_t size)
{
	constexpr uint32_t step = vectorWidth<T>();
	constexpr uint32_t align = step * 4;
	assert(size % align == 0);

	float *end = dst + size;
	T inf = T(Inf);
	for (; dst != end; dst += align) {
		vectorStore(dst + step*0, inf);
		vectorStore(dst + step*1, inf);
		vectorStore(dst + step*2, inf);
		vectorStore(dst + step*3, inf);
	}
}

template <typename T>
void renderSdfNew(Rasterizer &r)
{
	constexpr uint32_t step = vectorWidth<T>();
	uint32_t safeStep = max(step, 4u);

	r.safeHeight = r.height + safeStep * 2;
	r.safeWidth = r.width + safeStep * 2;
	uint32_t safeSize = r.safeWidth * r.safeHeight;

	uint32_t safeAlign = 4 * safeStep;
	safeSize = safeSize + (safeAlign - safeSize % safeAlign) % safeAlign;

	CombinedAllocation<16> alloc;

	alloc.add(r.distSqData[0], safeSize);
	alloc.add(r.offsetXData[0], safeSize);
	alloc.add(r.offsetYData[0], safeSize);

	alloc.add(r.distSqData[1], safeSize);
	alloc.add(r.offsetXData[1], safeSize);
	alloc.add(r.offsetYData[1], safeSize);

	alloc.add(r.winding, safeSize);
	alloc.add(r.rasterXs, r.safeWidth);
	alloc.add(r.rasterYs, r.safeHeight);

	alloc.setPointerOrAllocate(r.options.scratchMemory, r.options.scratchMemorySize);

	clearInf<T>(r.distSqData[0].data, safeSize);
	// clearInf<T>(r.offsetXData[0].data, safeSize);
	// clearInf<T>(r.offsetYData[0].data, safeSize);
	memset(r.winding.data, 0, r.winding.size * sizeof(int32_t));

	r.distSq = r.distSqData[0].data + (r.safeWidth + 1) * step;
	r.offsetX = r.offsetXData[0].data + (r.safeWidth + 1) * step;
	r.offsetY = r.offsetYData[0].data + (r.safeWidth + 1) * step;
	r.distSqDst = r.distSqData[1].data + (r.safeWidth + 1) * step;
	r.offsetXDst = r.offsetXData[1].data + (r.safeWidth + 1) * step;
	r.offsetYDst = r.offsetYData[1].data + (r.safeWidth + 1) * step;
	r.padding = 2;

	float maximumX = (float)(r.width - 1);
	float maximumY = (float)(r.height - 1);

	uint32_t width = r.width;
	uint32_t height = r.height;

	float *rasterXs = r.rasterXs.data;
	float *rasterYs = r.rasterYs.data;

	// Precompute X values
	for (uint32_t px = 0; px < width; px += step) {
		T x = T(float(px)) + vectorSequence<T>();
		x = transformX(r, x);
		x = nudgeOffInteger(x);
		vectorStore(rasterXs + px, x);
	}

	// Precompute Y values
	for (uint32_t py = 0; py < height; py += step) {
		T y = T(float(py)) + vectorSequence<T>();
		y = transformY(r, y);
		y = nudgeOffInteger(y);
		vectorStore(rasterYs + py, y);
	}

	computeContourWinding<T>(r);
	expandWinding(r);
	computeContourDistance<T>(r);

	// -- Vertices
	for (const Line &line : r.font.lines) {
		vertexPass(r, line.a);
		vertexPass(r, line.b);
	}
	for (const Bezier &bezier : r.font.beziers) {
		vertexPass(r, bezier.a);
		vertexPass(r, bezier.c);
	}

	// scanFillConservative(r);
	// scanFill(r);

	finalizeDistance<T>(r);
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
	// renderSdf<float>(r);
	renderSdf<SseFloat4>(r);
}
