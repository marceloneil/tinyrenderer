#include <algorithm>
#include <iostream>
#include <limits>

#include "model.h"
#include "tgaimage.h"

using namespace std;

const TGAColor white = TGAColor(255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0);
const TGAColor green = TGAColor(0, 255, 0);
const TGAColor blue = TGAColor(0, 0, 255);

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

void rasterize(vec2 v0, vec2 v1, TGAImage &image, const TGAColor &color, int ybuffer[], int renderHeight) {
    // make the line left-to-right
    if (v1.x < v0.x) {
        swap(v0, v1);
    }

    for (int x = v0.x; x <= v1.x; x += 1) {
        // calculate the y-value of the line for this x-value
        float t = (x - v0.x) / (float) (v1.x - v0.x);
        int y = v0.y + (v1.y - v0.y) * t;

        // only draw a pixel if it's the current closest pixel to the camera for that x-value
        if (ybuffer[x] < y) {
            ybuffer[x] = y;

            for (int renderY = 0; renderY < renderHeight; renderY += 1) {
                image.set(x, renderY, color);
            }
        }
    }
}

int main(int argc __attribute__((unused)), char* argv[] __attribute__((unused))) {
    const int width = 800;
    const int height = 500;

    // three line segments from triangles intersecting with the plane
    vec2 l1[2] = {vec2(20, 34), vec2(744, 400)};
    vec2 l2[2] = {vec2(120, 434), vec2(444, 400)};
    vec2 l3[2] = {vec2(330, 463), vec2(594, 200)};

    // screen line segment (below the triangles) interesecting with the plane
    vec2 sl[2] = {vec2(10, 10), vec2(790, 10)};

    // draw 2D scene of intersection with plane
    {
        TGAImage scene(width, height, TGAImage::RGB);

        line(l1[0], l1[1], scene, red);
        line(l2[0], l2[1], scene, green);
        line(l3[0], l3[1], scene, blue);
        line(sl[0], sl[1], scene, white);

        scene.write_tga_file("scene.tga");
    }

    // draw the render, a top down "1D" view of the slice of triangles
    {
        const int renderHeight = 16;
        TGAImage render(width, renderHeight, TGAImage::RGB);

        // initialize a ybuffer
        int ybuffer[width];
        for (int x = 0; x < width; x += 1) {
            ybuffer[x] = numeric_limits<int>::min();
        }

        rasterize(l1[0], l1[1], render, red, ybuffer, renderHeight);
        rasterize(l2[0], l2[1], render, green, ybuffer, renderHeight);
        rasterize(l3[0], l3[1], render, blue, ybuffer, renderHeight);
        rasterize(sl[0], sl[1], render, white, ybuffer, renderHeight);

        render.write_tga_file("render.tga");
    }
}
