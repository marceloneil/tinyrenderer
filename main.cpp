#include "tgaimage.h"

using namespace std;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
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
    // absolute slope. derror <= 1 pixel as we know that dy <= dx
    const float derror = abs(dy / float(dx));

    float error = 0;
    int y = y0;
    for (int x = x0; x <= x1; x += 1) {
        if (steep) {
            image.set(y, x, color); // de-transpose
        } else {
            image.set(x, y, color);
        }

        error += derror;
        if (error > 0.5) { // error exceeds a half-pixel
            y += y1 > y0 ? 1 : -1;
            error -= 1.0;
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
