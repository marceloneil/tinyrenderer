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

int main(int argc __attribute__((unused)), char* argv[] __attribute__((unused))) {
    TGAImage image(100, 100, TGAImage::RGB);
    for (int i = 0; i < 1e6; i += 1) {
        line(13, 20, 80, 40, image, white);
        line(20, 13, 40, 80, image, red);
        line(80, 40, 13, 20, image, red);
    }
    image.write_tga_file("output.tga");
}
