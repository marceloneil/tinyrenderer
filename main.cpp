#include <iostream>

#include "model.h"
#include "tgaimage.h"

using namespace std;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);

void line(int x0, int y0, int x1, int y1, TGAImage &image, const TGAColor &color) {
    bool steep = false;

    // if the line is steep, take it's transpose
    if (abs(x1 - x0) < abs(y1 - y0)) {
        steep = true;
        swap(x0, y0);
        swap(x1, y1);
    }

    // make the line left-to-right
    if (x1 < x0) {
        swap(x0, x1);
        swap(y0, y1);
    }

    const int dx = x1 - x0;
    const int dy = y1 - y0;
    const int derror2 = abs(dy) * 2;

    int error2 = 0;
    int y = y0;
    for (int x = x0; x <= x1; x += 1) {
        // if the line is steep, de-transpose
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }

        // slope error calculations
        error2 += derror2;
        if (error2 > dx) {
            y += y1 > y0 ? 1 : -1;
            error2 -= dx * 2;
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

    for (int face = 0; face < model.nfaces(); face += 1) {
        // draw a line from each vertex in a face to the next one
        // there are three vertices in each face
        for (int vertex = 0; vertex < 3; vertex += 1) {
            vec3 v0 = model.vert(face, vertex);
            vec3 v1 = model.vert(face, (vertex + 1) % 3);

            // normalize vertex coordinates which are between -1.0 and 1.0
            int x0 = (v0.x / 2.0 + 0.5) * width;
            int y0 = (v0.y / 2.0 + 0.5) * height;
            int x1 = (v1.x / 2.0 + 0.5) * width;
            int y1 = (v1.y / 2.0 + 0.5) * height;

            line(x0, y0, x1, y1, image, white);
        }
    }

    image.write_tga_file("output.tga");
}
