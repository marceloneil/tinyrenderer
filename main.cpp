#include <algorithm>
#include <iostream>

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

int main(int argc, char* argv[]) {
    string filename = "obj/african_head/african_head.obj";
    int width = 800;
    int height = 800;

    try {
        switch (argc) {
            case 4:
                height = stoi(argv[3]);
                if (height < 1) throw 1;
                // fall through
            case 3:
                width = stoi(argv[2]);
                if (width < 1) throw 1;
                // fall through
            case 2:
                filename = argv[1];
                // fall through
            case 1:
                break;
            default:
                throw 1;
        }
    } catch (...) {
        cerr << "Usage: " << argv[0]
             << " [ filename [ width (> 0) [ height (> 0) ] ] ]" << endl;
    }

    Model model(filename);
    TGAImage image(width, height, TGAImage::RGB);

    // light vector faces directly into the image
    vec3 lightVector(0, 0, -1);

    // draw each face
    for (int face = 0; face < model.nfaces(); face += 1) {
        // populate 3D "world" vertices and 2D "screen" vertices
        vec3 worldVertices[3];
        vec2 screenVertices[3];
        for (int vIdx = 0; vIdx < 3; vIdx += 1) {
            worldVertices[vIdx] = model.vert(face, vIdx);

            // adjust world vertices which are between -1.0 and 1.0 to be
            // within the dimensions of the image
            screenVertices[vIdx] = vec2(
                (worldVertices[vIdx].x / 2.0 + 0.5) * width,
                (worldVertices[vIdx].y / 2.0 + 0.5) * height
            );
        }

        // calculate the surface normal of the face (counterclockwise)
        // the surface normal points out of the visible side of the face
        vec3 edge01 = worldVertices[1] - worldVertices[0];
        vec3 edge02 = worldVertices[2] - worldVertices[0];
        vec3 surfaceNormal = cross(edge01, edge02).normalize();

        // calculate light intensity as the negative dot product of the surface normal
        // and the light vector
        // the light intensity will be positive if and only if the light vector
        // is pointing in the opposite z-direction of the normal. in other words,
        // it will be positive if the light vector points into the visible part of
        // the face, and it's intensity is determined by what angle it hits the face
        double lightIntensity = -(surfaceNormal * lightVector);

        // back-face culling
        // skip face if the surface normal faces away from the direction of light
        if (lightIntensity <= 0) continue;

        // draw
        int shade = lightIntensity * 255;
        triangle(
            screenVertices[0], screenVertices[1], screenVertices[2],
            image, TGAColor(shade, shade, shade)
        );
    }

    image.write_tga_file("output.tga");
}
