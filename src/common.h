#pragma once

#include <xmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>

#include <cmath>
#include <stdint.h>
#include <string.h>
#include <assert.h>

// -- Scalar

template <typename T>
struct VectorTraits { };

template <>
struct VectorTraits<float> {
	static constexpr uint32_t width = 1;

	typedef float FloatType;
	typedef int32_t IntType;
};

template <typename T>
using VectorInt = typename VectorTraits<T>::IntType;

inline bool anyTrue(bool mask) { return mask; }
inline bool anyFalse(bool mask) { return !mask; }
inline float maskSelect(bool mask, float a, float b) { return mask ? a : b; }

inline float min(float a, float b) { return a < b ? a : b; }
inline float max(float a, float b) { return a < b ? b : a; }

template <typename T>
inline T min(T a, T b) { return a < b ? a : b; }
template <typename T>
inline T max(T a, T b) { return a < b ? b : a; }

template <typename T>
inline T vectorLoad(const float *memory) { return *memory; }

template <typename T>
inline VectorInt<T> vectorLoad(const int32_t *memory) { return *memory; }

inline void vectorStore(float *memory, float v) { *memory = v; }

inline int32_t vectorToBits(float v) {
	int32_t r;
	memcpy(&r, &v, sizeof(float));
	return r;
}

inline float vectorFromBits(int32_t v) {
	float r;
	memcpy(&r, &v, sizeof(int32_t));
	return r;
}

inline float vectorSum(float v) { return v; }

template <typename T>
inline T vectorSequence() { return 0.0f; }

inline uint32_t vectorBitmask(bool v) { return v ? 1 : 0; }

inline void vectorStoreUint(uint32_t *dst, float t) {
	dst[0] = (uint32_t)(int32_t)t;
}

template <typename T>
constexpr uint32_t vectorWidth() { return VectorTraits<T>::width; }

template <typename T>
inline T flipSignIfNonZero(T value, const int32_t *ints) {
	auto i = vectorLoad<T>(ints);
	return maskSelect(i == 0, value, -value);
}

// -- SSE 2

struct SseFloat4
{
	struct Mask
	{
		__m128 v;

		Mask(__m128 v) : v(v) { }

		Mask operator!() const { return _mm_cmpeq_ps(v, _mm_setzero_ps()); }
		Mask operator&(const Mask &rhs) const { return _mm_and_ps(v, rhs.v); }
		Mask operator|(const Mask &rhs) const { return _mm_or_ps(v, rhs.v); }
		Mask operator^(const Mask &rhs) const { return _mm_xor_ps(v, rhs.v); }
	};

	__m128 v;

	SseFloat4() { }
	SseFloat4(float v) : v(_mm_set1_ps(v)) { }
	SseFloat4(__m128 v) : v(v) { }
	SseFloat4(float a, float b, float c, float d) : v(_mm_set_ps(d, c, b, a)) { }

	inline float x() { return _mm_cvtss_f32(v); }
	inline float y() { return _mm_cvtss_f32(_mm_shuffle_ps(v, v, _MM_SHUFFLE(3,2,1,1))); }
	inline float z() { return _mm_cvtss_f32(_mm_shuffle_ps(v, v, _MM_SHUFFLE(3,2,1,2))); }
	inline float w() { return _mm_cvtss_f32(_mm_shuffle_ps(v, v, _MM_SHUFFLE(3,2,1,3))); }

	SseFloat4 operator-() const { return _mm_sub_ps(_mm_setzero_ps(), v); }

	SseFloat4 operator+(SseFloat4 rhs) const { return _mm_add_ps(v, rhs.v); }
	SseFloat4 operator-(SseFloat4 rhs) const { return _mm_sub_ps(v, rhs.v); }
	SseFloat4 operator*(SseFloat4 rhs) const { return _mm_mul_ps(v, rhs.v); }
	SseFloat4 operator/(SseFloat4 rhs) const { return _mm_div_ps(v, rhs.v); }

	Mask operator==(SseFloat4 rhs) const { return _mm_cmpeq_ps(v, rhs.v); }
	Mask operator!=(SseFloat4 rhs) const { return _mm_cmpneq_ps(v, rhs.v); }
	Mask operator<(SseFloat4 rhs) const { return _mm_cmplt_ps(v, rhs.v); }
	Mask operator>(SseFloat4 rhs) const { return _mm_cmpgt_ps(v, rhs.v); }
	Mask operator<=(SseFloat4 rhs) const { return _mm_cmple_ps(v, rhs.v); }
	Mask operator>=(SseFloat4 rhs) const { return _mm_cmpge_ps(v, rhs.v); }
};

struct SseInt4
{
	struct Mask
	{
		__m128i v;

		Mask(__m128i v) : v(v) { }

		operator SseFloat4::Mask() { return _mm_castsi128_ps(v); }
	};

	__m128i v;

	SseInt4() { }
	SseInt4(int32_t v) : v(_mm_set1_epi32(v)) { }
	SseInt4(__m128i v) : v(v) { }

	SseInt4 operator+(SseInt4 rhs) const { return _mm_add_epi32(v, rhs.v); }
	SseInt4 operator-(SseInt4 rhs) const { return _mm_sub_epi32(v, rhs.v); }
	SseInt4 operator*(SseInt4 rhs) const { return _mm_mullo_epi32(v, rhs.v); }

	SseInt4 operator&(SseInt4 rhs) const { return _mm_and_si128(v, rhs.v); }
	SseInt4 operator|(SseInt4 rhs) const { return _mm_or_si128(v, rhs.v); }
	SseInt4 operator^(SseInt4 rhs) const { return _mm_xor_si128(v, rhs.v); }

	SseInt4 operator>>(int32_t shift) const { return _mm_srai_epi32(v, shift); }

	Mask operator==(SseInt4 rhs) const { return _mm_cmpeq_epi32(v, rhs.v); }
};

template <>
struct VectorTraits<SseFloat4> {
	static constexpr uint32_t width = 4;
	typedef SseFloat4 FloatType;
	typedef SseInt4 IntType;
};

inline SseFloat4 sqrt(SseFloat4 v) { return _mm_sqrt_ps(v.v); }
inline SseFloat4 abs(SseFloat4 v) { return _mm_andnot_ps(_mm_set1_ps(-0.0f), v.v); }
inline SseFloat4 floor(SseFloat4 v) { return _mm_floor_ps(v.v); }
inline SseFloat4 rint(SseFloat4 v) { return _mm_round_ps(v.v, _MM_FROUND_RINT); }
inline SseFloat4 min(SseFloat4 a, SseFloat4 b) { return _mm_min_ps(a.v, b.v); }
inline SseFloat4 max(SseFloat4 a, SseFloat4 b) { return _mm_max_ps(a.v, b.v); }
inline SseFloat4 copysign(SseFloat4 v, SseFloat4 sign) {
	__m128 signBit = _mm_set1_ps(-0.0f);
	return _mm_or_ps(_mm_andnot_ps(signBit, v.v), _mm_and_ps(signBit, sign.v));
}
inline SseFloat4 pow(SseFloat4 a, float b) {
#if 1
	float x = powf(a.x(), b);
	float y = powf(a.y(), b);
	float z = powf(a.z(), b);
	float w = powf(a.w(), b);
	return { x, y, z, w };
#else
	return pow_ps(a.v, _mm_set1_ps(b));
#endif
}
inline bool anyTrue(SseFloat4::Mask mask) {
	__m128i t = _mm_castps_si128(mask.v);
	return !_mm_test_all_zeros(t, t);
}
inline bool anyFalse(SseFloat4::Mask mask) {
	__m128i t = _mm_castps_si128(mask.v);
	return !_mm_test_all_ones(t);
}
inline SseFloat4 maskSelect(SseFloat4::Mask mask, SseFloat4 a, SseFloat4 b) {
	return _mm_blendv_ps(b.v, a.v, mask.v);
}

template<>
inline SseFloat4 vectorLoad<SseFloat4>(const float *memory) {
	return _mm_loadu_ps(memory);
}
template<>
inline SseInt4 vectorLoad<SseFloat4>(const int32_t *memory) {
	return _mm_loadu_si128((const __m128i*)memory);
}
inline void vectorStore(float *memory, SseFloat4 v) {
	_mm_storeu_ps(memory, v.v);
}

inline SseInt4 vectorToBits(SseFloat4 v) { return _mm_castps_si128(v.v); }
inline SseFloat4 vectorFromBits(SseInt4 v) { return _mm_castsi128_ps(v.v); }

inline float vectorSum(SseFloat4 v) {
	// TODO: Optimize
	return v.x() + v.y() + v.z() + v.w();
}

template <>
inline SseFloat4 vectorSequence<SseFloat4>() {
	return SseFloat4(0.0f, 1.0f, 2.0f, 3.0f);
}

inline uint32_t vectorBitmask(SseFloat4::Mask v) {
	return (uint32_t)_mm_movemask_ps(v.v);
}

inline void vectorStoreUint(uint32_t *dst, SseFloat4 v) {
	_mm_storeu_si128((__m128i*)dst, _mm_cvtps_epi32(v.v));
}

// -- Math

// https://www.musicdsp.org/en/latest/Other/206-fast-cube-root-square-root-and-reciprocal-for-x86-sse-cpus.html
template <typename T>
inline T approxCbrt(T t) {
	auto bits = vectorToBits(t);
	auto sign = bits & 0x80000000;
	bits = (bits & 0x7fffffff) - 0x3f800000; // unbias
	bits = (bits >> 10) * 341; // logarithmic multiply by 1/3
	bits = (bits + 0x3f800000) & 0x7fffffff; // rebias
	bits = bits | sign;
	T z = vectorFromBits(bits);

	// Newton–Raphson steps
	z = z - (z*z*z - t) / (z*z*3.0f);
	// z = z - (z*z*z - t) / (z*z*3.0f);

	return z;
}

// -- Vector

template <typename T>
inline T clamp(T v, T minV, T maxV) {
	return min(max(v, minV), maxV);
}

template <typename T>
inline T saturate(T v) {
	return clamp(v, T(0.0f), T(1.0f));
}

template <typename T>
struct Vec2
{
	T x, y;

	Vec2(): x(), y() { }
	Vec2(T x, T y) : x(x), y(y) { }

	Vec2<T> operator+(Vec2<T> rhs) const { return { x + rhs.x, y + rhs.y }; }
	Vec2<T> operator-(Vec2<T> rhs) const { return { x - rhs.x, y - rhs.y }; }
	Vec2<T> operator*(Vec2<T> rhs) const { return { x * rhs.x, y * rhs.y }; }
	Vec2<T> operator/(Vec2<T> rhs) const { return { x / rhs.x, y / rhs.y }; }
	Vec2<T> operator*(T rhs) const { return { x * rhs, y * rhs }; }
	Vec2<T> operator/(T rhs) const { return { x / rhs, y / rhs }; }

	template <typename U>
	Vec2<U> cast() const { return Vec2<U> { U(x), U(y) }; }

	Vec2 yx() const { return { y, x }; }
};

template <typename T>
inline Vec2<T> vec2(T x, T y) { return { x, y }; }

template <typename T>
inline T dot(Vec2<T> a, Vec2<T> b) {
	return a.x*b.x + a.y*b.y;
}

template <typename T>
inline T dot2(Vec2<T> a) {
	return a.x*a.x + a.y*a.y;
}

template <typename T>
inline Vec2<T> min(const Vec2<T> &a, const Vec2<T> &b) { return { min(a.x, b.x), min(a.y, b.y) }; }
template <typename T>
inline Vec2<T> max(const Vec2<T> &a, const Vec2<T> &b) { return { max(a.x, b.x), max(a.y, b.y) }; }
template <typename T>
inline Vec2<T> clamp(const Vec2<T> &a, const Vec2<T> &minV, const Vec2<T> &maxV) { return { clamp(a.x, minV.x, maxV.x), clamp(a.y, minV.y, maxV.y) }; }

static const float Inf = (float)INFINITY;

using Vec2f = Vec2<float>;
