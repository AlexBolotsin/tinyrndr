#include <vector>
#include <cmath>
#include <limits>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

Model *model = NULL;
const int width = 800;
const int height = 800;
const int depth = 255;
int *zbuffer = NULL;
Vec3f light_dir(0, 0, -1);

void triangle(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, TGAImage &image, float intensity, int* zbuffer) {
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
        float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height;
        Vec3f A = Vec3f(t0) + Vec3f(t2 - t0)*alpha;
        Vec3f B = second_half ? Vec3f(t1) + Vec3f(t2 - t1)*beta : Vec3f(t0) + Vec3f(t1 - t0)*beta;
        Vec2i uvA = uv0 + (uv2 - uv0)*alpha;
        Vec2i uvB = second_half ? uv1 + (uv2 - uv1)*beta : uv0 + (uv1 - uv0)*beta;

        if (A.x > B.x) std::swap(A, B);

        for (int j = A.x; j <= B.x; j++)
        {
            float phi = B.x == A.x ? 1.f : (float)(j - A.x) / (float)(B.x - A.x);
            Vec3f p = A + (B - A)*phi;
            int idx = p.x + p.y * width;
            Vec2i uvP = uvA + (uvB - uvA)*phi;

            if (zbuffer[idx] < p.z)
            {
                zbuffer[idx] = p.z;
                TGAColor color = model->diffuse(uvP);
                image.set(p.x, p.y, TGAColor(color.r*intensity, color.g*intensity, color.b*intensity, 255));
            }
        }
    }
}

int main(int argc, char** argv) {
    if (2 == argc) {
        model = new Model(argv[1]);
    }
    else {
        model = new Model("obj/african_head.obj");
    }
    zbuffer = new int[width*height];
    for (int i = 0; i < width*height; i++) {
        zbuffer[i] = std::numeric_limits<int>::min();
    }

    {
        TGAImage image(width, height, TGAImage::RGB);
        for (int i = 0; i < model->nfaces(); i++) {
            std::vector<int> face = model->face(i);
            Vec3i screen_coords[3];
            Vec3f world_coords[3];
            for (int j = 0; j < 3; j++) {
                Vec3f v = model->vert(face[j]);
                screen_coords[j] = Vec3i((v.x + 1.f)*width / 2.f, (v.y + 1.f)*height / 2.f, (v.z + 1.f)*depth / 2.f);
                world_coords[j] = v;
            }
            Vec3f n = cross<float>(world_coords[2] - world_coords[0], world_coords[1] - world_coords[0]);
            n.normalize();
            float intensity = n * light_dir;
            if (intensity > 0)
            {
                Vec2i uv[3];
                for (int k = 0; k < 3; k++) {
                    uv[k] = model->uv(i, k);
                }
                triangle(screen_coords[0], screen_coords[1], screen_coords[2], uv[0], uv[1], uv[2], image, intensity, zbuffer);
            }
        }
        image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
        image.write_tga_file("output.tga");
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

    delete model;
    delete[] zbuffer;

    system("pause");
    return 0;
}

