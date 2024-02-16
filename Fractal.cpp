#include <fstream>
#include <complex.h>
#include <vector>
#include <omp.h>

namespace
{
  constexpr float radius_sqr = 4.f; // squared of escape radius
  constexpr int bitdepth = 16;
  constexpr int rgb_pixel_size = 3;
  constexpr int max_thead_num = 8;

  // implementation using complex numbers
  int mandelbrot_comp(const float real, const float imag, const int max_count)
  {
    const std::complex<float> c(real, imag);
    std::complex<float> z(real, imag);
    int count = 0;
    while (std::norm(z) < radius_sqr and ++count < max_count)
    {
      z = z * z + c;
    }
    return count;
  };

  int julia_comp(const float real, const float imag, const int max_count)
  {
    const std::complex<float> c(-0.8, 0.156);
    std::complex<float> z(real, imag);
    int count = 0;
    while (std::norm(z) < radius_sqr and ++count < max_count)
    {
      z = z * z + c;
    }
    return count;
  }
  
  // implementation using c-like approach  
  int mandelbrot(const float real, const float imag, const int max_count)
  {
    float z_real = real;
    float z_imag = imag;
    int count = 0;
    while (z_real * z_real + z_imag * z_imag < radius_sqr and ++count < max_count)
    {
      const float temp = z_real * z_real - z_imag * z_imag + real;
      z_imag = 2.0f * z_real * z_imag + imag;
      z_real = temp;
    }
    return count;
  };

  int julia(const float real, const float imag, const int max_count)
  {
    const float c_real = -0.8;  //-0.4;
    const float c_imag = 0.156; // 0.6;

    float z_real = real;
    float z_imag = imag;
    int count = 0;
    while (z_real * z_real + z_imag * z_imag < radius_sqr and ++count < max_count)
    {
      const float temp = z_real * z_real - z_imag * z_imag;
      z_imag = c_imag + 2.f * z_real * z_imag;
      z_real = c_real + temp;
    }
    return count;
  }
};

// generate fractal set using function pointer style
using FractalRule = int (*)(const float, const float, const int);
void generate_fractal_set(const FractalRule rule,
                          const float x_start, const float y_start,
                          const float x_end, const float y_end,
                          const int nx /*width*/, const int ny /*height*/,
                          const int max_count, std::vector<int> &fractal_data)
{
  // 1. determine increment steps
  const float dx = (x_end - x_start) / nx; // cols
  const float dy = (y_end - y_start) / ny; // rows
  fractal_data.resize(0);
  fractal_data.resize(nx * ny, 0);

// 2. loop over coord x, y
#pragma omp parallel for shared(fractal_data)
  for (int j = 0; j < ny; ++j)
  {
    int index = j * nx;
    float x = x_start;
    const float y = y_start + j * dy;
    for (int i = 0; i < nx; ++i)
    {
      fractal_data[index++] = rule(x, y, max_count);
      x += dx;
    }
  }
};

void int_to_rgb(std::vector<char> &out, const std::vector<int> &in)
{
  out.resize(0);
  out.reserve(in.size() * rgb_pixel_size);
  char r, g, b;

  for (const auto var : in)
  {
    // 1. convert int to r, g, b
    const auto temp = var / bitdepth;
    r = temp / bitdepth;
    g = temp % bitdepth;
    b = var % bitdepth;

    // 2. store
    out.push_back(r);
    out.push_back(g);
    out.push_back(b);
  }
}

void write_ppm(const std::string &name, const std::vector<char> &shades, const int nx /*width columns*/, const int ny /*height rows*/)
{
  // 1. validate input data
  std::ofstream stream(name.c_str(), std::ios::out | std::ios::binary);

  if (!stream.is_open())
  {
    throw std::runtime_error("could not open the file");
  }
  
  if (shades.size() != nx * ny * rgb_pixel_size)
  {
    throw std::runtime_error("image buffer size and image size mismatch");
  }

  // 2. write the magic number, image size and bit depth
  stream << "P6\n"
         << nx << ' ' << ny << '\n'
         << bitdepth << '\n';

  // 3. write image data
  stream.write(shades.data(), shades.size() * sizeof(char));
  stream.close();
}

int main()
{
  // 1. data buffer
  std::vector<int> fractal_data;
  std::vector<char> image_data;

  constexpr int nx = 1000;
  constexpr int max_iter = 255;

  // 2. mandel set
  {
    constexpr float minx = -2.5;
    constexpr float miny = -2.0;
    constexpr float maxx = 1.5;
    constexpr float maxy = 2.0;
    constexpr int ny = (maxy - miny) / (maxx - minx) * nx;
    generate_fractal_set(mandelbrot, minx, miny, maxx, maxy, nx, ny, max_iter, fractal_data);
    int_to_rgb(image_data, fractal_data);
    write_ppm("mandelbrot.ppm", image_data, nx, ny);
  }

  // 3. julia set
  {
    constexpr float minx = -2.0;
    constexpr float miny = -2.0;
    constexpr float maxx = 2.0;
    constexpr float maxy = 2.0;
    constexpr int ny = (maxy - miny) / (maxx - minx) * nx;
    generate_fractal_set(julia, minx, miny, maxx, maxy, nx, ny, max_iter, fractal_data);
    int_to_rgb(image_data, fractal_data);
    write_ppm("julia.ppm", image_data, nx, ny);
  }

  return 0;
}
