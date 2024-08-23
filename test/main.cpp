#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdint.h>

#include "../src/rasterizer.h"

#define STB_TRUETYPE_IMPLEMENTATION
#include "../external/stb_truetype.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../external/stb_image_write.h"

#define CPUTIME_IMPLEMENTATION
#include "../external/cputime.h"

unsigned char buffer[24<<20];

unsigned char scratchMemory[24<<22];

#if 0
#error PLAN: \
	To avoid internal distances we must do the winding pass first to identify filled areas. \
	With that we can ignore any shapes which reduce the winding when the pixel in questions is still wound that way. \
	ie. pixel in winding number 1 ignores line giving it -1 winding in that direction. \
	Similarly, if the magnitude of the winding number is >=2 then we must ignore any distances. \
	This still seems a bit sketchy but hopefully the filling algorithm can cover for that, and even then the worst that \
	can happen is some bad antialiasing.
#endif

#error TODO: Figure out cheap spreading, can be janky as it should only apply to controversial areas
#error Also SIMD vertex pass or something

void test()
{
	// fread(buffer, 1, sizeof(buffer), fopen("c:/windows/fonts/arialbd.ttf", "rb"));
	// fread(buffer, 1, sizeof(buffer), fopen("C:\\Unity\\Kabe\\Assets\\UI\\Fonts\\Darumadrop_One\\DarumadropOne-Regular.ttf", "rb"));
	fread(buffer, 1, sizeof(buffer), fopen("C:\\Unity\\Kabe\\Assets\\UI\\Fonts\\Noto_Sans\\NotoSansJP-Medium.ttf", "rb"));

	stbtt_fontinfo font_info;
	stbtt_InitFont(&font_info, buffer, 0);

	int codepoint = 0x6f22;
	// int codepoint = 'P';

	stbtt_vertex *vertices;
	int num_vertices = stbtt_GetCodepointShape(&font_info, codepoint,  &vertices);

	Vec2<int16_t> contourStart, prevPos;

	Font font;
	bool verbose = false;

	Vec2f boundsMin = vec2(Inf, Inf);
	Vec2f boundsMax = vec2(-Inf, -Inf);

	for (int i = 0; i < num_vertices; i++) {
		stbtt_vertex v = vertices[i];
		Vec2<int16_t> pos = vec2(v.x, v.y);
		Vec2<int16_t> control = vec2(v.cx, v.cy);

		boundsMin = min(boundsMin, pos.cast<float>());
		boundsMax = max(boundsMax, pos.cast<float>());

		switch (v.type) {
		case STBTT_vmove:
			if (verbose) printf("==\nmove %d %d\n", v.x, v.y);
			if (prevPos.x != contourStart.x || prevPos.y != contourStart.y) {
				font.lines.push_back(Line{
					prevPos.cast<float>(),
					contourStart.cast<float>(),
				});
			}

			contourStart = pos;
			break;
		case STBTT_vline:
			if (verbose) printf("line %d %d\n", v.x, v.y);
			font.lines.push_back(Line{
				prevPos.cast<float>(),
				pos.cast<float>(),
			});
			break;
		case STBTT_vcurve:
			if (verbose) printf("curve %d %d %d %d\n", v.x, v.y, v.cx, v.cy);
			font.beziers.push_back(Bezier{
				prevPos.cast<float>(),
				control.cast<float>(),
				pos.cast<float>(),
			});
			break;
		case STBTT_vcubic:
			assert(0);
			break;
		}

		prevPos = pos;
	}

	uint32_t width = 128;
	uint32_t height = 0;
	if (height == 0) height = width;
	RasterizeOptions opts = { };

	opts.scratchMemory = scratchMemory;
	opts.scratchMemorySize = sizeof(scratchMemory);

	Vec2f boundsExtent = boundsMax - boundsMin;
	float boundsSize = max(boundsExtent.x, boundsExtent.y) * 1.2f;
	Vec2f boundsMid = (boundsMin + boundsMax) * 0.5f;

	opts.offset = vec2(boundsMid.x - boundsSize * 0.5f, boundsMid.y + boundsSize * 0.5f);
	opts.scale = vec2(boundsSize / float(width), -boundsSize / float(height));

	std::vector<float> distances;

	cputime_begin_init();

	uint64_t minTime = UINT64_MAX;

	uint32_t runs = 10000;

	for (uint32_t i = 0; i < runs; i++) {
		distances.clear();
		distances.resize(width * height, 1.0f);

		uint64_t rasterizeBegin = cputime_cpu_tick();
		rasterizeFont(font, opts, width, height, distances.data());
		uint64_t rasterizeEnd = cputime_cpu_tick();
		uint64_t duration = rasterizeEnd - rasterizeBegin;
		if (minTime > duration) {
			minTime = duration;
		}
	}

	cputime_end_init();
	printf("Took %.3fms (%ux%u)\n", cputime_cpu_delta_to_sec(NULL, minTime) * 1e3, width, height);

#if 1
	for (float &d : distances) {
		d = d / 100.0f + 0.5f;
	}
#endif

	stbtt_FreeShape(&font_info, vertices);

	std::vector<uint8_t> pixels;
	pixels.resize(width * height * 4);

	for (uint32_t i = 0; i < width * height; i++)
	{
		uint8_t v = (uint8_t)clamp(distances[i] * 255.5f, 0.0f, 255.0f);
		pixels[i * 4 + 0] = v;
		pixels[i * 4 + 1] = v;
		pixels[i * 4 + 2] = v;
		pixels[i * 4 + 3] = 0xff;
	}

	stbi_write_png("result.png", (int)width, (int)height, 4, pixels.data(), 0);

	return;

#if 1
	{
		unsigned char *sdf = nullptr;
		int w, h, x, y;
		uint64_t stbMinTime = UINT64_MAX;

		for (uint32_t i = 0; i < runs; i++) {
			if (sdf) stbtt_FreeSDF(sdf, NULL);

			uint64_t stbBegin = cputime_cpu_tick();
			sdf = stbtt_GetCodepointSDF(&font_info, 0.25f, codepoint, 20, 128, 1.0f, &w, &h, &x, &y);
			uint64_t stbEnd = cputime_cpu_tick();
			uint64_t duration = stbEnd - stbBegin;
			if (stbMinTime > duration) {
				stbMinTime = duration;
			}
		}

		printf("stbtt took %.2fms (%dx%d)\n", cputime_cpu_delta_to_sec(NULL, stbMinTime) * 1e3, w, h);
		stbi_write_png("sdf.png", w, h, 1, sdf, 0);
	}
#endif

	int refWidth, refHeight, refX, refY;
	unsigned char *ref = stbtt_GetCodepointBitmap(&font_info, 1.0f, 1.0f, codepoint, &refWidth, &refHeight, &refX, &refY);
	stbi_write_png("ref.png", (int)refWidth, (int)refHeight, 1, ref, 0);

}

int main(int argc, char **argv)
{
	test();
	return 0;
}
