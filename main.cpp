#include <algorithm>
#include <iostream>

#include "model.h"
#include "tgaimage.h"

using namespace std;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
const TGAColor blue = TGAColor(0, 0, 255, 255);

void line(vec2 v0, vec2 v1, TGAImage &image, const TGAColor &color) {
    bool steep = false;

    // if the line is steep, take it's transpose
    if (abs(v1.x - v0.x) < abs(v1.y - v0.y)) {
        steep = true;
        swap(v0.x, v0.y);
        swap(v1.x, v1.y);
    }

    // make the line left-to-right
    if (v1.x < v0.x) {
        swap(v0, v1);
    }

    const int dx = v1.x - v0.x;
    const int dy = v1.y - v0.y;
    const int derror2 = abs(dy) * 2;

    int error2 = 0;
    int y = v0.y;
    for (int x = v0.x; x <= v1.x; x += 1) {
        // if the line is steep, de-transpose
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }

        // slope error calculations
        error2 += derror2;
        if (error2 > dx) {
            y += v1.y > v0.y ? 1 : -1;
            error2 -= dx * 2;
        }
    }
}

void triangle(vec2 t0, vec2 t1, vec2 t2, TGAImage &image, const TGAColor &color) {
    // ignore degenerate triangles
    if (t0.y == t1.y && t1.y == t2.y) return;

    // sort vertices by y in ascending order
    if (t0.y > t1.y) swap(t0, t1);
    if (t0.y > t2.y) swap(t0, t2);
    if (t1.y > t2.y) swap(t1, t2);

    int totalHeight = t2.y - t0.y;
    for (int y = t0.y; y <= t0.y + totalHeight; y += 1) {
        bool secondHalf = y >= t1.y;
        int segmentHeight = secondHalf ? t2.y - t1.y : t1.y - t0.y;

        float alpha = (y - t0.y) / float(totalHeight);
        float beta = 0;
        if (segmentHeight > 0) {
            beta = (secondHalf ? y - t1.y : y - t0.y) / float(segmentHeight);
        }

        vec2 A = t0 + (t2 - t0) * alpha;
        vec2 B = secondHalf ? t1 + (t2 - t1) * beta : t0 + (t1 - t0) * beta;

        if (B.x < A.x) swap(A, B);

        for (int x = A.x; x <= B.x; x += 1) {
            image.set(x, y, color);
        }
    }
}

int main(int argc __attribute__((unused)), char* argv[] __attribute__((unused))) {
    TGAImage image(200, 200, TGAImage::RGB);

    vec2 t0[3] = {vec2(10, 70), vec2(50, 160), vec2(70, 80)};
    vec2 t1[3] = {vec2(180, 50), vec2(150, 1), vec2(70, 180)};
    vec2 t2[3] = {vec2(180, 150), vec2(120, 160), vec2(130, 180)};
    vec2 t3[3] = {vec2(25, 25), vec2(50, 50), vec2(75, 25)};

    triangle(t0[0], t0[1], t0[2], image, red);
    triangle(t1[0], t1[1], t1[2], image, white);
    triangle(t2[0], t2[1], t2[2], image, green);
    triangle(t3[0], t3[1], t3[2], image, blue);

    image.write_tga_file("output.tga");
}
