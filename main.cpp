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
    if (t0.x == t1.x && t1.x == t2.x) return;

    // bubble sort vertices by y in ascending order
    vec2 &bottom = t0, &middle = t1, &top = t2;
    if (bottom.y > middle.y) swap(top, middle);
    if (bottom.y > top.y) swap(bottom, top);
    if (middle.y > top.y) swap(middle, top);

    // the triangle is made up of three line segments:
    //  - the "long" segment, which runs directly from bottom to top
    //    it's the longest segment on the y axis, not necessarily the longest overall
    //  - the "bottom" segment, which runs from bottom to middle
    //  - the "top" segment, which runs from middle to top
    vec2 longSegment = top - bottom;
    assert(longSegment.y > 1e-3); // height of triangle is non-zero
    vec2 bottomSegment = middle - bottom;
    vec2 topSegment = top - middle;

    // sweep lines from bottom to top
    for (int y = bottom.y; y <= top.y; y += 1) {
        // find the point where the long segment intersects with y
        float longSegmentRatio = (y - bottom.y) / longSegment.y;
        vec2 longSegmentIntercept = bottom + longSegment * longSegmentRatio;

        // find the point where the short segment intersects with y
        vec2 shortSegmentIntercept;
        if (y <= middle.y && bottom.y != middle.y) {
            // bottom segment
            float bottomSegmentRatio = (y - bottom.y) / bottomSegment.y;
            shortSegmentIntercept = bottom + bottomSegment * bottomSegmentRatio;
        } else {
            // top segment
            float topSegmentRatio = (y - middle.y) / topSegment.y;
            shortSegmentIntercept = middle + topSegment * topSegmentRatio;
        }

        // find the left and right x-intercepts for this y value
        int xInterceptLeft = round(longSegmentIntercept.x);
        int xInterceptRight = round(shortSegmentIntercept.x);
        if (xInterceptLeft > xInterceptRight) swap(xInterceptLeft, xInterceptRight);

        // sweep pixels from left to right
        for (int x = xInterceptLeft; x <= xInterceptRight; x += 1) {
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
    vec2 t4[3] = {vec2(55, 50), vec2(80, 25), vec2(105, 50)};

    triangle(t0[0], t0[1], t0[2], image, red);
    triangle(t1[0], t1[1], t1[2], image, white);
    triangle(t2[0], t2[1], t2[2], image, green);
    triangle(t3[0], t3[1], t3[2], image, blue);
    triangle(t4[0], t4[1], t4[2], image, blue);

    image.write_tga_file("output.tga");
}
