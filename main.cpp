#include "tgaimage.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
    for (int x = x0; x <= x1; x += 1) {
        float t = (x - x0) / (float) (x1 - x0);
        int y = y0 + (y1 - y0) * t;
        image.set(x, y, color);
    }
}

int main(int argc __attribute__((unused)), char* argv[] __attribute__((unused))) {
    TGAImage image(100, 100, TGAImage::RGB);
    line(13, 20, 80, 40, image, white); 
    line(20, 13, 40, 80, image, red); 
    line(80, 40, 13, 20, image, red);
    image.write_tga_file("output.tga");
}
