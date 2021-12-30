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

// draw a triangle by filling in lines of pixels between the line segments of the triangle
void triangle_lineSweep(vec2 t0, vec2 t1, vec2 t2, TGAImage &image, const TGAColor &color) {
    // ignore degenerate triangles
    if (t0.x == t1.x && t1.x == t2.x) return;
    if (t0.y == t1.y && t1.y == t2.y) return;

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

    // clamp bottom and top of triangle to height of image
    int yBottomClamped = max((int) round(bottom.y), 0);
    int yTopClamped = min((int) round(top.y), image.get_height() - 1);

    // sweep lines from bottom to top
    for (int y = yBottomClamped; y <= yTopClamped; y += 1) {
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

        // clamp left and right intercepts to width of image
        xInterceptLeft = max(xInterceptLeft, 0);
        xInterceptRight = min(xInterceptRight, image.get_width() - 1);

        // sweep pixels from left to right
        for (int x = xInterceptLeft; x <= xInterceptRight; x += 1) {
            image.set(x, y, color);
        }
    }
}

struct BarycentricCoords {
    double u, v, w;
};

// find the barycentric coordinates using a cross product
BarycentricCoords barycentric_crossProduct(vec2 A, vec2 B, vec2 C, vec2 P) {
    // solve the equation:
    //      (1 - u - v) * A + u * B + v * C = P
    // <=>  u * (B - A) + v * (C - A) + (A - P) = 0
    // <=>  u * (B_x - A_x) + v * (C_x - A_x) + (A_x - P_x) = 0
    //      u * (B_y - A_y) + v * (C_y - A_y) + (A_y - P_y) = 0
    // <=>  [u, v, 1] * [B_x - A_x, C_x - A_x, A_x - P_x] = 0
    //      [u, v, 1] * [B_y - A_y, C_y - A_y, A_y - P_y] = 0
    // <=>  [B_x - A_x, C_x - A_x, A_x - P_x] ⨯ [B_y - A_y, C_y - A_y, A_y - P_y] = [u, v, 1]
    vec3 solution = cross(
        vec3(B.x - A.x, C.x - A.x, A.x - P.x),
        vec3(B.y - A.y, C.y - A.y, A.y - P.y)
    );

    // divide solution by some scaling factor 1 * solution[2] 
    double u = solution[0] / solution[2];
    double v = solution[1] / solution[2];
    double w = 1 - u - v;

    return BarycentricCoords{u, v, w};
}

// find the barycentric coordinates by solving for the edges using a 2D inverse matrix
BarycentricCoords barycentric_edgeMatrix(vec2 A, vec2 B, vec2 C, vec2 P) {
    // solve the equation:
    //      (1 - u - v) * A + u * B + v * C = P
    // <=>  u * (B - A) + v * (C - A) = P - A
    // <=>  u * (B_x - A_x) + v * (C_x - A_x) = P_x - A_x
    //      u * (B_y - A_y) + v * (C_y - A_y) = P_y - A_y
    // we can create a matrix with this system of equations, giving us:
    //      system * [u, v] = P - A
    // if we multiply from the left by the inverse of this matrix, we have:
    //      [u, v] = system^-1 * (P - A)
    // we also know that w = 1 - u - v
    mat<2,2> system = {
        vec2(B.x - A.x, C.x - A.x),
        vec2(B.y - A.y, C.y - A.y)
    };
    vec2 solution = system.invert() * (P - A);

    double u = solution[0];
    double v = solution[1];
    double w = 1 - u - v;

    return BarycentricCoords{u, v, w};
}

// find the barycentric coordinates by solving for the vertices using a 3D inverse matrix
BarycentricCoords barycentric_vertexMatrix(vec2 A, vec2 B, vec2 C, vec2 P) {
    // solve the system of equations:
    //      u * B_x + v * C_x + w * A_x = P_x
    //      u * B_y + v * C_y + w * A_y = P_y
    //      u + v + w = 1
    // we can create a matrix with this system, giving us the equation:
    //      system * [u, v, w] = [P_x, P_y, 1]
    // if we multiply from the left by the inverse of this matrix, we have:
    //      [u, v, w] = system^-1 * [P_x, P_y, 1]
    mat<3,3> system = {
        vec3(B.x, C.x, A.x),
        vec3(B.y, C.y, A.y),
        vec3(1, 1, 1)
    };
    vec3 result = vec3(P.x, P.y, 1);
    vec3 solution = system.invert() * result;

    double u = solution[0];
    double v = solution[1];
    double w = solution[2];

    return BarycentricCoords{u, v, w};
}

// find the barycentric coordinates using signed ratios of areas
BarycentricCoords barycentric_area(vec2 A, vec2 B, vec2 C, vec2 P) {
    // the signed areas are calculated in the basis (AB, AC, k)
    vec2 AB = B - A;
    vec2 AC = C - A;
    vec3 k = vec3(0, 0, 1); // unit vector, as in (i, j, k)

    // we can derive an equation for AP in this basis
    vec2 AP = P - A;

    // the equation is given by:
    //      AP = 1 / (AB, AC, k) * ((AP, AC, k) * AB + (AB, AP, k) * AC + (AB, AC, AP) * k)
    // where (e, f, g) = (e ⨯ f) * g is the mixed product
    // note that:
    //  - (AB, AC, k) = (AB ⨯ AC) * k is twice the signed area of the triangle ABC
    //  - (AP, AC, k) = (AP ⨯ AC) * k is twice the signed area of the triangle APC
    //  - (AB, AP, k) = (AB ⨯ AP) * k is twice the signed area of the triangle ABP
    //  - (AB, AC, AP) = (AB ⨯ AC) * AP is zero, as AP is orthogonal AB ⨯ AC
    // we can then rearrange our equation as
    //      AP = 1 / (2 * ABC) * ((2 * APC) * AB + (2 * ABP) * AC)
    //      AP = (APC / ABC) * AB + (ABP / ABC) * AC
    //      P = A + (APC / ABC) * AB + (ABP / ABC) * AC
    //      P = (1 - APC / ABC - ABP / ABC) * A + (APC / ABC) * B + (ABP / ABC) * C
    // finally, we have derived equations for the barycentric coordinates in terms of
    // the signed areas of ABC, APC and ABP:
    //      u = APC / ABC
    //      v = ABP / ABC
    //      w = 1 - u - v

    // the signed area of the triangle ABC is given by (AB ⨯ AC) * k / 2
    // the cross product is calculated counterclockwise within our basis (A,B,C)
    double ABC = cross(vec3(AB.x, AB.y, 0), vec3(AC.x, AC.y, 0)) * k / 2;

    // the signed area of the triangle APC is given by (AP ⨯ AC) * k / 2
    // the cross product is calculated counterclockwise within our basis (A,P,C)
    // if P is outside of the triangle ABC, this area will be negative since
    // the cross product AP ⨯ AC will be facing in the direction -k
    double APC = cross(vec3(AP.x, AP.y, 0), vec3(AC.x, AC.y, 0)) * k / 2;

    // the signed area of the triangle ABP is given by (AB ⨯ AP) * k / 2
    // the cross product is calculated counterclockwise within our basis (A,B,P)
    // if P is outside of the triangle ABC, this area will be negative since
    // the cross product AB ⨯ AP will be facing in the direction -k
    double ABP = cross(vec3(AB.x, AB.y, 0), vec3(AP.x, AP.y, 0)) * k / 2;

    // the coordinates are given by the signed area ratio of the opposite subtriangle
    // for u * B this is ACP, for v * C this is ABP, and for w * A this is BCP
    double u = APC / ABC;
    double v = ABP / ABC;
    double w = 1 - u - v;

    return BarycentricCoords{u, v, w};
}

// find the barycentric coordinates directly
// note: all other methods can be simplified down to this method
BarycentricCoords barycentric_simple(vec2 A, vec2 B, vec2 C, vec2 P) {
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

inline BarycentricCoords barycentric(vec2 A, vec2 B, vec2 C, vec2 P) {
    // assume valid triangles
    assert(A.x != B.x || B.x != C.x);
    assert(A.y != B.y || B.y != C.y);

#if defined(BARYCENTRIC_CROSSPRODUCT)
    return barycentric_crossProduct(A, B, C, P);
#elif defined(BARYCENTRIC_EDGEMATRIX)
    return barycentric_edgeMatrix(A, B, C, P);
#elif defined(BARYCENTRIC_VERTEXMATRIX)
    return barycentric_vertexMatrix(A, B, C, P);
#elif defined(BARYCENTRIC_AREA)
    return barycentric_area(A, B, C, P);
#elif defined(BARYCENTRIC_SIMPLE)
    return barycentric_simple(A, B, C, P);
#else
    #error no barycentric method specified
#endif
}

// draw a triangle by testing if individual pixels are within the triangle
void triangle_barycentric(vec2 t0, vec2 t1, vec2 t2, TGAImage &image, const TGAColor &color) {
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

inline void triangle(vec2 t0, vec2 t1, vec2 t2, TGAImage &image, const TGAColor &color) {
#if defined(TRIANGLE_LINESWEEP)
    triangle_lineSweep(t0, t1, t2, image, color);
#elif defined(TRIANGLE_BARYCENTRIC)
    triangle_barycentric(t0, t1, t2, image, color);
#else
    #error no triangle implementation specified
#endif
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
