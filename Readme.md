# Introduction
This is a 2D fractal generator coming with two predefined rules: Mandelbrot and Julia. A '*.ppm' image file visualises the computed result.
You can experiment with other fractal rules. A simple OpenMP directive is also included to enable parallel runs.

# To compile \& run
For Linux & MacOs users:

g++ --std=c++17 -O3 Fractal.cpp -o Fractal

To enable OpenMP:

g++ --std=c++17 -fopenmp -O3 Fractal.cpp -o Fractal

To run:

./Fractal

# License
The example code in this directory is licensed under the terms of the MIT license. 
See [LICENSE-MIT](LICENSE-MIT) for details.
