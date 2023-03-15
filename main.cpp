#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <eigen3/Eigen/eigen>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
using namespace Eigen;
const float radius = 1.5;
Vector3f center = { 0,0,0 };
const float noise_amplitude = 1.0;

template <typename T> inline T lerp(const T& v0, const T& v1, float t) {
    return v0 + (v1 - v0) * std::max(0.f, std::min(1.f, t));
}

float hash(const float n) {// 伪随机数  
    float x = sin(n) * 43758.5453f;
    return x - floor(x);
}

float noise(const Vector3f& x) {// 噪声函数
    Vector3f p(floor(x.x()), floor(x.y()), floor(x.z()));
    Vector3f f(x.x() - p.x(), x.y() - p.y(), x.z() - p.z()); // t = x - x0
    f = f * (f.dot(Vector3f(3.f, 3.f, 3.f) - f * 2.f)); // 插值函数 t = 3t^2 - 2t^3
    float n = p.dot(Vector3f(1.f, 57.f, 113.f));
    return lerp(lerp(
        lerp(hash(n + 0.f), hash(n + 1.f), f.x()),
        lerp(hash(n + 57.f), hash(n + 58.f), f.x()), f.y()),
        lerp(
            lerp(hash(n + 113.f), hash(n + 114.f), f.x()),
            lerp(hash(n + 170.f), hash(n + 171.f), f.x()), f.y()), f.z());
}

Vector3f rotate(const Vector3f& v) {
    return Vector3f(Vector3f(0.00, 0.80, 0.60).dot(v), Vector3f(-0.80, 0.36, -0.48).dot(v), Vector3f(-0.60, -0.48, 0.64).dot(v));
}

float fractal_brownian_motion(const Vector3f& x) {// 分型布朗运动this is a bad noise function with lots of artifacts. TODO: find a better one
    Vector3f p = rotate(x);                              // 把噪声叠加
    float f = 0;
    f += 0.5000 * noise(p); p = p * 2.32;
    f += 0.2500 * noise(p); p = p * 3.03;
    f += 0.1250 * noise(p); p = p * 2.61;
    f += 0.0625 * noise(p);
    return f / 0.9375;
}
float signed_distance(const Vector3f& p) { // 距离物体表面距离的函数 f(x, y, z) = x^2+y^2+z^2 - r^2 - sin(x)*sin(y)*sin(z)
    //Vector3f s = Vector3f(p).normalized() * radius;/*
    //float displacement = sin(16 * s.x()) * sin(16 * s.y()) * sin(16 * s.z()) * noise_a*/mplitude;
    float displacement = -fractal_brownian_motion(p * 3.4) * noise_amplitude;
    return p.norm() - (radius + displacement);
}
//bool sphere_trace(const Vector3f& orig, const Vector3f& dir, Vector3f& pos) {
//    float ra = radius;
//    if (pos != Vector3f(0, 0, 0)) {
//        ra = ra + sin(16 * pos.x()) * sin(16 * pos.y()) * sin(16 * pos.z()) * noise_amplitude;
//    }
//    Vector3f L = center - orig;
//    float L2 = L.dot(L);
//    float tca = L.dot(dir);
//    float d2 = L.dot(L) - tca * tca;
//    if (d2 > ra * ra) return false;
//    if (L2 > ra * ra && tca < 0)return false;
//    float thc = sqrtf(ra * ra - d2);
//    float t0 = tca - thc;// 球内钝角，球外
//    float t1 = tca + thc; // 球内锐角
//    if (t0 < 0) t0 = t1;
//    pos = orig + dir * t0;
//    return true;
//}
bool sphere_trace(const Vector3f& orig, const Vector3f& dir, Vector3f& pos) {
    pos = orig;
    if (orig.dot(orig) - pow(orig.dot(dir), 2) > pow(signed_distance(pos), 2)) return false;
    for (size_t i = 0; i < 128; i++) {
        float d = signed_distance(pos);
        if (d < 0) return true;
        pos = pos + dir * std::max(d * 0.1f, .01f);
    }
    return false;
}
Vector3f distance_field_normal(const Vector3f& pos) {// 通过点的微扰来算法向量
    const float eps = 0.1;
    float d = signed_distance(pos);
    float nx = signed_distance(pos + Vector3f(eps, 0, 0)) - d;
    float ny = signed_distance(pos + Vector3f(0, eps, 0)) - d;
    float nz = signed_distance(pos + Vector3f(0, 0, eps)) - d;
    return Vector3f(nx, ny, nz).normalized();
}

Vector3f palette_fire(const float d) {
    const Vector3f   yellow(1.7, 1.3, 1.0); // note that the color is "hot", i.e. has components >1
    const Vector3f   orange(1.0, 0.6, 0.0);
    const Vector3f      red(1.0, 0.0, 0.0);
    const Vector3f darkgray(0.2, 0.2, 0.2);
    const Vector3f     gray(0.4, 0.4, 0.4);

    float x = std::max(0.f, std::min(1.f, d));
    if (x < .25f)
        return lerp(gray, darkgray, x * 4.f);
    else if (x < .5f)
        return lerp(darkgray, red, x * 4.f - 1.f);
    else if (x < .75f)
        return lerp(red, orange, x * 4.f - 2.f);
    return lerp(orange, yellow, x * 4.f - 3.f);
}

int main() {
    int cnt = 0;
    const int   width = 640;
    const int   height = 480;
    const float fov = M_PI / 3.;
    std::vector<Vector3f> framebuffer(width * height);
    const Vector3f eye_pos = { 0,0,3 };
    float scale = std::tan(fov * 0.5f*M_PI / 180.0);
    float imageAspectRatio = width / (float)height;
    #pragma omp parallel for
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            float dir_x = (i + 0.5) - width / 2.;
            float dir_y = -(j + 0.5) + height / 2.;
            float dir_z = -height / (2. * tan(fov / 2.));
            Vector3f hit = { 0,0,0 };
            if (sphere_trace(Vector3f(0, 0, 3), Vector3f(dir_x, dir_y, dir_z).normalized(), hit)) {
                Vector3f light_dir = (Vector3f(10, 10, 10) - hit).normalized();  //光源(10,10,10)
                float noise_level = (radius - hit.norm()) / noise_amplitude;
                hit = distance_field_normal(hit);
                float light_intensity = std::max(0.4f, light_dir.dot(hit));
                framebuffer[i + j * width] = palette_fire((-.2 + noise_level) * 2) * light_intensity;
            }
            else {
                framebuffer[i + j * width] = Vector3f(0.2, 0.7, 0.8); // 背景颜色
            }
        }
        std::cout << cnt++<<std::endl;
    }

    std::vector<unsigned char> pixmap(width * height * 3);
    for (size_t i = 0; i < height * width; ++i) {
        Vector3f& c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1) c = c * (1. / max);
        for (size_t j = 0; j < 3; j++) {
            pixmap[i * 3 + j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    std::string fileName = "./output/out_.jpg";
    stbi_write_jpg(fileName.c_str(), width, height, 3, pixmap.data(), 100);
}