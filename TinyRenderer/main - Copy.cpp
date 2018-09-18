#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
const TGAColor blue = TGAColor(0, 0, 255, 255);
Model *model = NULL;
const int width  = 800;
const int height = 800;
const int depth = 255;
int *zbuffer = NULL;

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
    bool steep = false;
    if (std::abs(x0-x1)<std::abs(y0-y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0>x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    for (int x=x0; x<=x1; x++) {
        float t = (x-x0)/(float)(x1-x0);
        int y = int(y0*(1.-t) + y1*t);
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
    }
}

void triangle(Vec3i t0, Vec3i t1, Vec3i t2, TGAImage &image, TGAColor color, int* zbuffer) {
    if (t0.y == t1.y && t0.y == t2.y) return;

    if (t0.y > t1.y) std::swap(t0, t1);
    if (t0.y > t2.y) std::swap(t0, t2);
    if (t1.y > t2.y) std::swap(t1, t2);
    int total_height = t2.y - t0.y;

    for (int i = 0; i < total_height; i++)
    {
        bool second_half = i > t1.y - t0.y || t1.y == t0.y;
        int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
        float alpha = (float)i / total_height;
        float beta = (float)(i - (second_half ? t1.y-t0.y : 0)) / segment_height;
        Vec3i A = CastToVec3i(CastToVec3f(t0) + CastToVec3f(t2 - t0)*alpha);
        Vec3i B = second_half ? CastToVec3i(CastToVec3f(t1) + CastToVec3f(t2 - t1)*beta) : CastToVec3i(CastToVec3f(t0) + CastToVec3f(t1 - t0)*beta);

        if (A.x > B.x) std::swap(A, B);
        
        for (int j = A.x; j <= B.x; j++)
        {
            float phi = B.x == A.x ? 1.f : (float)(j - A.x) / (float)(B.x - A.x);
            Vec3f p = CastToVec3f(A) + CastToVec3f(B - A)*phi;
            int idx = p.x + p.y * width;

            if (zbuffer[idx] < p.z)
            {
                image.set(p.x, p.y, color);
            }
        }
    }
}

int main(int argc, char** argv) {
    TGAImage image(width, height, TGAImage::RGB);
    Vec3f light_dir(0, 0, -1);
    zbuffer = new int[width*height];
    for (int i = 0; i < width*height; i++) {
        zbuffer[i] = std::numeric_limits<int>::min();
    }

    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    for (int i=0; i<model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        Vec3i screen_coords[3];
        Vec3f world_coords[3];
        for (int j=0; j<3; j++) {
            Vec3f v = model->vert(face[j]);
            screen_coords[j] = Vec3i((v.x + 1.f)*width / 2.f, (v.y + 1.f)*height / 2.f, (v.z+1.f)*depth/2.f);
            world_coords[j] = v;
        }
        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
        n.normalize();
        float intensity = n * light_dir;
        if (intensity > 0)
        {
            triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255), zbuffer);
        }
    }
    {
        TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                zbimage.set(i, j, TGAColor(zbuffer[i + j * width], 1));
            }
        }
        zbimage.flip_vertically(); // i want to have the origin at the left bottom corner of the image
        zbimage.write_tga_file("zbuffer.tga");
    }

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;
    delete[] zbuffer;

    system("pause");
    return 0;
}

