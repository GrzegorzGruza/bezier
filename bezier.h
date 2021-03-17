#ifndef BEZIER_H
#define BEZIER_H

#include <cstddef>
#include <ostream>
#include <vector>
#include <functional>
#include <cmath>
#include <iostream>

namespace bezier {

    namespace constants {
        const int NUM_OF_CUBIC_BEZIER_NODES = 4;
    } // bezier::constants

    class Point2d;

    namespace types {
        using point_2d = Point2d;
        using real_t = double;
        using node_index_t = unsigned int;
    } // bezier::types

    namespace details {
        const types::real_t ARC = 4 * (sqrt(2) - 1) / 3;

        types::real_t degToRad(types::real_t deg) {
            return deg * M_PI / 180;
        }

        using curve_t = std::function<types::point_2d(types::node_index_t)>;
        using oneDimArray_t = std::vector<bool>;
        using twoDimArray_t = std::vector<std::vector<bool>>;
    } // bezier::details

    class Point2d {
    public:
        Point2d(types::real_t X, types::real_t Y) : X(X), Y(Y) {};

        // A linear transformation to a given point.
        Point2d(const Point2d &point,
                types::real_t a11, types::real_t a12,
                types::real_t a21, types::real_t a22) :
                X(point.X * a11 + point.Y * a12),
                Y(point.X * a21 + point.Y * a22) {};

        Point2d operator+(const Point2d &other) const {
            return Point2d{this->X + other.X, this->Y + other.Y};
        }

        Point2d operator*(const types::real_t scalar) const {
            return Point2d{this->X * scalar, this->Y * scalar};
        }

        bool operator==(const Point2d &other) const {
            return (this->X == other.X && this->Y == other.Y);
        }

        std::ostream &operator<<(std::ostream &ostr) const {
            ostr << "(" << this->X << ", " << this->Y << ")";
            return ostr;
        }

        const types::real_t X;
        const types::real_t Y;
    };

    Point2d operator*(const types::real_t scalar, const Point2d &point) {
        return point * scalar;
    }

    types::point_2d nodes(const types::point_2d p0,
                          const types::point_2d p1,
                          const types::point_2d p2,
                          const types::point_2d p3,
                          types::node_index_t index) {
        return (index == 0 ? p0 :
                (index == 1 ? p1 :
                 (index == 2 ? p2 :
                  (index == 3 ? p3 :
                   throw std::out_of_range(
                           "a curve node index is out of range")))));
    }

    details::curve_t Cup() {
        static auto cup =
                [](types::node_index_t index) {
                    return nodes(types::point_2d(-1, 1),
                                 types::point_2d(-1, -1),
                                 types::point_2d(1, -1),
                                 types::point_2d(1, 1),
                                 index);
                };
        return cup;
    }

    details::curve_t Cap() {
        static auto cap =
                [](types::node_index_t index) {
                    return nodes(types::point_2d(-1, -1),
                                 types::point_2d(-1, 1),
                                 types::point_2d(1, 1),
                                 types::point_2d(1, -1),
                                 index);
                };
        return cap;
    }

    details::curve_t ConvexArc() {
        static auto convexArc =
                [](types::node_index_t index) {
                    return nodes(types::point_2d(0, 1),
                                 types::point_2d(details::ARC, 1),
                                 types::point_2d(1, details::ARC),
                                 types::point_2d(1, 0),
                                 index);
                };
        return convexArc;
    }

    details::curve_t ConcaveArc() {
        static auto concaveArc =
                [](types::node_index_t index) {
                    return nodes(types::point_2d(0, 1),
                                 types::point_2d(0, 1 - details::ARC),
                                 types::point_2d(1 - details::ARC, 0),
                                 types::point_2d(1, 0),
                                 index);
                };
        return concaveArc;
    }

    details::curve_t LineSegment(
            types::point_2d p, types::point_2d q) {
        auto lineSegment =
                [p, q](types::node_index_t index) {
                    return nodes(p, p, q, q, index);
                };
        return lineSegment;
    }


    details::curve_t MovePoint(
            const details::curve_t &f,
            types::node_index_t index_no,
            types::real_t x,
            types::real_t y) {
        auto movePoint =
                [f, index_no, x, y](types::node_index_t index) {
                    return (index == index_no
                            ? f(index) + types::point_2d(x, y)
                            : f(index));
                };
        return movePoint;
    }

    details::curve_t Rotate(
            const details::curve_t &f,
            types::real_t a) {
        auto rotate =
                [f, a](types::node_index_t index) {
                    return types::point_2d(f(index),
                                           cos(details::degToRad(a)),
                                           -sin(details::degToRad(a)),
                                           sin(details::degToRad(a)),
                                           cos(details::degToRad(a)));
                };
        return rotate;
    }

    details::curve_t Scale(
            const details::curve_t &f,
            types::real_t x,
            types::real_t y) {
        auto scale =
                [f, x, y](types::node_index_t index) {
                    return types::point_2d(f(index),
                                           x, 0,
                                           0, y);
                };
        return scale;
    }

    details::curve_t Translate(
            const details::curve_t &f,
            types::real_t x,
            types::real_t y) {
        auto translate =
                [f, x, y](types::node_index_t index) {
                    return f(index) + types::point_2d(x, y);
                };
        return translate;
    }

    details::curve_t Concatenate(
            const details::curve_t &f1,
            const details::curve_t &f2) {
        auto concatenate =
                [f1, f2](types::node_index_t index) {
                    return (index < 4 ? f1(index) : f2(index - 4));
                };
        return concatenate;
    }

    template<typename... Arg>
    details::curve_t Concatenate(
            details::curve_t f1, Arg... args) {
        return Concatenate(f1, Concatenate(args...));
    }

    class P3CurvePlotter {
    public:
        explicit P3CurvePlotter(
                const details::curve_t &f,
                size_t numberOfSegments = 1,
                size_t resolution = 80) :
                resolution(resolution) {
            belongsToPlot = details::twoDimArray_t
                    (resolution, details::oneDimArray_t(resolution, false));

            // Rounded up.
            size_t pointsPerSegment =
                    (resolution * resolution + numberOfSegments - 1)
                    / numberOfSegments;

            types::real_t t;

            size_t x_index, y_index;
            for (size_t sg_no = 1; sg_no <= numberOfSegments; sg_no++) {
                for (size_t i = 0; i <= pointsPerSegment; i++) {
                    t = 1.0 * i / pointsPerSegment;
                    types::point_2d p = this->operator()(f, t, sg_no);

                    // First shift by vector (1, 1).
                    x_index = (size_t) ((p.X + 1)
                                        * ((types::real_t) resolution / 2));
                    y_index = (size_t) ((p.Y + 1)
                                        * ((types::real_t) resolution / 2));

                    if (x_index < 0 || x_index >= resolution
                        || y_index < 0 || y_index >= resolution) {
                        continue;
                    }

                    // Printing a vector is done with lines from last to first.
                    belongsToPlot[(resolution - 1) - y_index][x_index] = true;
                }
            }
        }

        void Print(std::ostream &os = std::cout,
                   const char fb = '*',
                   const char bg = ' ') const {

            for (const auto &row : belongsToPlot) {
                for (auto element : row) {
                    if (element) {
                        os << fb;
                    } else {
                        os << bg;
                    }
                }
                os << std::endl;
            }
        }

        types::point_2d operator()(const details::curve_t &f,
                                   types::real_t t,
                                   size_t segmentNumber) const {

            const types::point_2d p0 = f(0 + (segmentNumber - 1) * 4);
            const types::point_2d p1 = f(1 + (segmentNumber - 1) * 4);
            const types::point_2d p2 = f(2 + (segmentNumber - 1) * 4);
            const types::point_2d p3 = f(3 + (segmentNumber - 1) * 4);

            types::point_2d b0 = (((types::real_t) 1.0 - t) * p0) + (t * p1);
            types::point_2d b1 = (((types::real_t) 1.0 - t) * p1) + (t * p2);
            types::point_2d b2 = (((types::real_t) 1.0 - t) * p2) + (t * p3);

            types::point_2d b0_2 = (((types::real_t) 1.0 - t) * b0) + (t * b1);
            types::point_2d b1_2 = (((types::real_t) 1.0 - t) * b1) + (t * b2);

            types::point_2d b0_3 = (((types::real_t) 1.0 - t) * b0_2) + (t * b1_2);

            return b0_3;
        }

    private:
        size_t resolution;
        details::twoDimArray_t belongsToPlot;
    };

} // bezier

#endif // BEZIER_H