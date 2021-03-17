#include "bezier.h"

#include <cassert>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

int main() {
  // Koniczyna
  auto fn = bezier::Cap();
  fn = bezier::MovePoint(fn, 0, 0.75, 1.25);
  fn = bezier::MovePoint(fn, 3, -0.75, 1.25);
  fn = bezier::Concatenate(fn,
                           bezier::Rotate(fn, 270),
                           bezier::Rotate(fn, 180),
                           bezier::Rotate(fn, 90));
  fn = bezier::Scale(fn, 1.0, 0.5);
  const bezier::P3CurvePlotter plot1(fn, 4);
  plot1.Print(std::cout, '*', '.');
  std::ofstream file{"clover.out"};
  plot1.Print(file, 'o');

  // Gwiazda
  fn = bezier::ConcaveArc();
  fn = bezier::Concatenate(fn,
                           bezier::Rotate(fn, 270),
                           bezier::Rotate(fn, 180),
                           bezier::Rotate(fn, 90));
  const bezier::P3CurvePlotter plot2(bezier::Scale(bezier::Rotate(fn, 30), 1.0, 0.5), 4, 50);
  plot2.Print(std::cout, 'X', ' ');
  plot2.Print(std::cout, ' ', '0');

  // Buźka
  fn = bezier::ConvexArc();
  fn = bezier::Concatenate(fn,
                           bezier::Rotate(fn, 270),
                           bezier::Rotate(fn, 180),
                           bezier::Rotate(fn, 90),
                           bezier::Translate(bezier::Scale(bezier::Cup(), 0.625, 0.25), 0.0, -0.375),
                           bezier::Scale(bezier::Cap(), 0.125, 0.125),
                           bezier::LineSegment(bezier::types::point_2d{-0.625, 0.375},
                                               bezier::types::point_2d{-0.25, 0.5}),
                           bezier::LineSegment(bezier::types::point_2d{0.625, 0.375},
                                               bezier::types::point_2d{0.25, 0.5}));
  const bezier::P3CurvePlotter plot3(bezier::Scale(fn, 1.0, 0.5), 8, 60);
  plot3.Print();

  // Węzły (punkty kontrolne) krzywej – sprawdzenie poprawności węzłów
  constexpr bezier::types::real_t precision = std::pow(2, -16);
  bezier::types::real_t a, b;
  fn = bezier::Concatenate(bezier::Rotate(bezier::Cup(), 90),
                           bezier::MovePoint(bezier::MovePoint(bezier::Cap(), 0, 1.0, 0.0), 3, -1.0, 0.0));
  std::vector<std::pair<bezier::types::real_t, bezier::types::real_t>>
    v{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {0, -1}, {-1, 1}, {1, 1}, {0, -1}};
  for_each(v.begin(), v.end(),
    [&](auto & p) {
      static bezier::types::node_index_t k = 0;
      a = std::abs(p.first - fn(k).X);
      b = std::abs(p.second - fn(k).Y);
      assert(a <= precision && b <= precision);
      k++;
    }
  );
  std::cout << "Curve node correctness: PASS\n";

  // Sprawdzenie wyjątku
  fn = bezier::Concatenate(bezier::ConvexArc(),
                           bezier::ConcaveArc(),
                           bezier::Cap(),
                           bezier::Cup());
  try {
    a = fn(4 * bezier::constants::NUM_OF_CUBIC_BEZIER_NODES).X;
    std::cout << "Exception: FAIL\n";
  }
  catch (const std::out_of_range &) {
    std::cout << "Exception: PASS\n";
  }
  catch (...) {
    std::cout << "Exception: FAIL\n";
  }
}
