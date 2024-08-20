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

void test()
{
	fread(buffer, 1, sizeof(buffer), fopen("c:/windows/fonts/arialbd.ttf", "rb"));
	// fread(buffer, 1, sizeof(buffer), fopen("C:\\Unity\\Kabe\\Assets\\UI\\Fonts\\Darumadrop_One\\DarumadropOne-Regular.ttf", "rb"));

	stbtt_fontinfo font_info;
	stbtt_InitFont(&font_info, buffer, 0);

	stbtt_vertex *vertices;
	int num_vertices = stbtt_GetCodepointShape(&font_info, 'P',  &vertices);

	Vec2<int16_t> contourStart, prevPos;

	Font font;
	bool verbose = false;

	for (int i = 0; i < num_vertices; i++) {
		stbtt_vertex v = vertices[i];
		Vec2<int16_t> pos = vec2(v.x, v.y);
		Vec2<int16_t> control = vec2(v.cx, v.cy);

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

	uint32_t width = 256;
	uint32_t height = 256;
	RasterizeOptions opts = { };
	opts.offset = vec2(-200.0f, 700.0f);
	opts.scale = vec2(1000.0f / float(width), -1000.0f / float(height));

	std::vector<float> distances;

	cputime_begin_init();

	uint64_t minTime = UINT64_MAX;

	uint32_t runs = 5;

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
	printf("Took %.2fms\n", cputime_cpu_delta_to_sec(NULL, minTime) * 1e3);

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

#if 1
	{
		unsigned char *sdf = nullptr;
		int w, h, x, y;
		uint64_t stbMinTime = UINT64_MAX;

		for (uint32_t i = 0; i < runs; i++) {
			if (sdf) stbtt_FreeSDF(sdf, NULL);

			uint64_t stbBegin = cputime_cpu_tick();
			sdf = stbtt_GetCodepointSDF(&font_info, 0.25f, 'P', 20, 128, 1.0f, &w, &h, &x, &y);
			uint64_t stbEnd = cputime_cpu_tick();
			uint64_t duration = stbEnd - stbBegin;
			if (stbMinTime > duration) {
				stbMinTime = duration;
			}
		}

		printf("stbtt took %.2fms (%dx%d)\n", cputime_cpu_delta_to_sec(NULL, stbMinTime) * 1e3, w, h);
		// stbi_write_png("sdf.png", w, h, 1, sdf, 0);
	}
#endif

	int refWidth, refHeight, refX, refY;
	unsigned char *ref = stbtt_GetCodepointBitmap(&font_info, 1.0f, 1.0f, 'P', &refWidth, &refHeight, &refX, &refY);
	stbi_write_png("ref.png", (int)refWidth, (int)refHeight, 1, ref, 0);

}

int main(int argc, char **argv)
{
	test();
	return 0;
}
