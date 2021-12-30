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

struct BarycentricCoords {
    double u, v, w;
};

// find the barycentric coordinates of a point P given a triangle ABC
BarycentricCoords barycentric(vec2 A, vec2 B, vec2 C, vec2 P) {
    // assume valid triangles
    assert(A.x != B.x || B.x != C.x);
    assert(A.y != B.y || B.y != C.y);

    // calculate the edge vectors
    vec2 AB = B - A;
    vec2 AC = C - A;
    vec2 AP = P - A;

    double denom = AB.x * AC.y - AC.x * AB.y;
    double u = (AP.x * AC.y - AC.x * AP.y) / denom;
    double v = (AB.x * AP.y - AP.x * AB.y) / denom;
    double w = 1 - u - v;

    return BarycentricCoords{u, v, w};
}

// draw a triangle by testing if individual pixels are within the triangle
void triangle(vec2 t0, vec2 t1, vec2 t2, TGAImage &image, const TGAColor &color) {
    // ignore degenerate triangles
    if (t0.x == t1.x && t1.x == t2.x) return;
    if (t0.y == t1.y && t1.y == t2.y) return;

    // find lower left and upper right corners of minimal bounding box
    vec2 bboxmin(
        min({t0.x, t1.x, t2.x, image.get_width() - 1.0}),
        min({t0.y, t1.y, t2.y, image.get_height() - 1.0})
    );
    vec2 bboxmax(
        max({t0.x, t1.x, t2.x, 0.0}),
        max({t0.y, t1.y, t2.y, 0.0})
    );

    // clamp bounding box to within image dimensions
    int xStart = round(max(bboxmin.x, 0.0));
    int xEnd = round(min(bboxmax.x, image.get_width() - 1.0));
    int yStart = round(max(bboxmin.y, 0.0));
    int yEnd = round(min(bboxmax.y, image.get_height() - 1.0));

    // loop through all pixels in bounding box
#pragma omp parallel for
    for (int x = xStart; x <= xEnd; x += 1) {
        for (int y = yStart; y <= yEnd; y += 1) {
            // check that pixel is within the triangle (simplex)
            BarycentricCoords coords = barycentric(t0, t1, t2, vec2(x, y));
            if (coords.u < 0 || coords.v < 0 || coords.w < 0) continue;

            image.set(x, y, color);
        }
    }
}

int main(int argc __attribute__((unused)), char* argv[] __attribute__((unused))) {
    TGAImage image(200, 200, TGAImage::RGB);

    vec2 t0[3] = {vec2(10, 70), vec2(50, 160), vec2(70, 80)};
    vec2 t1[3] = {vec2(180, 50), vec2(150, 1), vec2(70, 180)};
    vec2 t2[3] = {vec2(180, 150), vec2(120, 160), vec2(130, 180)};
    vec2 t3[3] = {vec2(25, 25), vec2(30, 50), vec2(75, 25)};
    vec2 t4[3] = {vec2(65, 60), vec2(105, 15), vec2(110, 60)};
    vec2 t5[3] = {vec2(160, 120), vec2(250, 140), vec2(250, 50)};
    vec2 t6[3] = {vec2(-5, 150), vec2(-5, 230), vec2(150, 260)};

    // degenerate triangles
    vec2 d0[3] = {vec2(5, 5), vec2(5, 5), vec2(5, 5)};
    vec2 d1[3] = {vec2(5, 5), vec2(100, 5), vec2(150, 5)};
    vec2 d2[3] = {vec2(5, 5), vec2(5, 100), vec2(5, 150)};

    triangle(t0[0], t0[1], t0[2], image, red);
    triangle(t1[0], t1[1], t1[2], image, white);
    triangle(t2[0], t2[1], t2[2], image, green);
    triangle(t3[0], t3[1], t3[2], image, blue);
    triangle(t4[0], t4[1], t4[2], image, green);
    triangle(t5[0], t5[1], t5[2], image, red);
    triangle(t6[0], t6[1], t6[2], image, blue);

    triangle(d0[0], d0[1], d0[2], image, red);
    triangle(d1[0], d1[1], d1[2], image, red);
    triangle(d2[0], d2[1], d2[2], image, red);

    image.write_tga_file("output.tga");
}
