#ifndef BEZIER_H
#define BEZIER_H

#include <cstddef>
#include <ostream>
#include <cmath>
#include <array>
#include <vector>
#include <cstdarg>
#include <iostream>

namespace bezier {

    namespace constants {
        const int NUM_OF_CUBIC_BEZIER_NODES = 4;
    } // namespace constants

    struct Point2d;

    namespace types {
        using point_2d = Point2d;
        using real_t = double;
        using node_index_t = size_t;
    } // namespace types

    namespace details {
        const types::real_t ARC = 4 * (sqrt(2) - 1) / 3;

        types::real_t degreesToRadians(types::real_t degrees) {
            return degrees * M_PI / 180;
        }

        types::real_t sinDegrees(types::real_t angle) {
            return sin(degreesToRadians(angle));
        }

        types::real_t cosDegrees(types::real_t angle) {
            return cos(degreesToRadians(angle));
        }
    } // namespace details

    struct Point2d {
        Point2d(types::real_t x, types::real_t y) : X(x), Y(y) {}

        Point2d operator+(const Point2d &point) const {
            return Point2d{this->X + point.X, this->Y + point.Y};
        }

        Point2d operator*(types::real_t scalar) const {
            return Point2d{X * scalar, Y * scalar};
        }

        bool operator==(const Point2d &point) const {
            return (this->X == point.X && this->Y == point.Y);
        }

        friend std::ostream &operator<<(std::ostream &os, const Point2d &point);

        types::real_t X;
        types::real_t Y;
    };

    inline std::ostream &operator<<(std::ostream &os, const Point2d &point) {
        os << "(" << point.X << ", " << point.Y << ")";
        return os;
    }

    template<typename T>
    inline Point2d operator*(const T &scalar, const Point2d &point) {
        return point * scalar;
    }

    class Curve {
    public:
        Curve(const Point2d &p0,
              const Point2d &p1,
              const Point2d &p2,
              const Point2d &p3) : nodes{{p0, p1, p2, p3}} {}

        Curve(const Curve &curve1, const Curve &curve2) : nodes(curve1.nodes) {
            for (auto const &value : curve2.nodes) {
                nodes.push_back(value);
            }
        }

        const Point2d &operator()(types::node_index_t index) const {
            try {
                return nodes.at(index / constants::NUM_OF_CUBIC_BEZIER_NODES).at(
                        index % constants::NUM_OF_CUBIC_BEZIER_NODES);
            }
            catch (const std::out_of_range &) {
                throw std::out_of_range("a curve node index is out of range");
            }
        }

        friend Curve MovePoint(const Curve &curve,
                               types::node_index_t index,
                               types::real_t x,
                               types::real_t y);

        friend Curve Rotate(const Curve &curve, types::real_t angle);

        friend Curve Scale(const Curve &curve, types::real_t x, types::real_t y);

        friend Curve Translate(const Curve &curve,
                               types::real_t x,
                               types::real_t y);

    private:
        Point2d &at(types::node_index_t index) {
            try {
                return nodes.at(index / constants::NUM_OF_CUBIC_BEZIER_NODES).at(
                        index % constants::NUM_OF_CUBIC_BEZIER_NODES);
            }
            catch (const std::out_of_range &) {
                throw std::out_of_range("a curve node index is out of range");
            }
        }

    private:
        std::vector<std::array<Point2d, constants::NUM_OF_CUBIC_BEZIER_NODES>>
                nodes;
    };

    inline Curve Cup() {
        return Curve{Point2d{-1, 1}, Point2d{-1, -1},
                     Point2d{1, -1}, Point2d{1, 1}};
    }

    inline Curve Cap() {
        return Curve{Point2d{-1, -1}, Point2d{-1, 1},
                     Point2d{1, 1}, Point2d{1, -1}};
    }

    inline Curve ConvexArc() {
        return Curve{Point2d{0, 1}, Point2d{details::ARC, 1},
                     Point2d{1, details::ARC}, Point2d{1, 0}};
    }

    inline Curve ConcaveArc() {
        return Curve{Point2d{0, 1}, Point2d{0, 1 - details::ARC},
                     Point2d{1 - details::ARC, 0}, Point2d{1, 0}};
    }

    inline Curve LineSegment(const Point2d &p, const Point2d &q) {
        return Curve{p, p, q, q};
    }

    inline Curve MovePoint(const Curve &curve,
                           types::node_index_t index,
                           types::real_t x,
                           types::real_t y) {
        Curve newCurve(curve);
        Point2d &oldPoint = newCurve.at(index);
        oldPoint.X += x;
        oldPoint.Y += y;
        return newCurve;
    }

    inline Curve Rotate(const Curve &curve, types::real_t angle) {
        Curve newCurve(curve);
        for (types::node_index_t i = 0;
             i < newCurve.nodes.size() * constants::NUM_OF_CUBIC_BEZIER_NODES;
             i++) {
            Point2d &oldPoint = newCurve.at(i);
            Point2d newPoint
                    (oldPoint.X * details::cosDegrees(angle) - oldPoint.Y * details::sinDegrees(angle),
                     oldPoint.X * details::sinDegrees(angle) + oldPoint.Y * details::cosDegrees(angle));
            oldPoint = newPoint;
        }
        return newCurve;
    }

    inline Curve Scale(const Curve &curve, types::real_t x, types::real_t y) {
        Curve newCurve(curve);
        for (types::node_index_t i = 0;
             i < newCurve.nodes.size() * constants::NUM_OF_CUBIC_BEZIER_NODES;
             i++) {
            Point2d &oldPoint = newCurve.at(i);
            oldPoint.X *= x;
            oldPoint.Y *= y;
        }
        return newCurve;
    }

    inline Curve Translate(const Curve &curve, types::real_t x, types::real_t y) {
        Curve newCurve(curve);
        for (types::node_index_t i = 0;
             i < newCurve.nodes.size() * constants::NUM_OF_CUBIC_BEZIER_NODES;
             i++) {
            Point2d &oldPoint = newCurve.at(i);
            oldPoint.X += x;
            oldPoint.Y += y;
        }
        return newCurve;
    }

    inline Curve Concatenate(const Curve &curve1, const Curve &curve2) {
        return Curve(curve1, curve2);
    }

    template<typename T, typename... Rest>
    inline Curve Concatenate(const T &curve, Rest... rest) {
        return Concatenate(curve, Concatenate(rest...));
    }

    class P3CurvePlotter {
    public:
        P3CurvePlotter(const Curve &curve,
                       size_t numberOfSegments,
                       size_t resolution) : resolution(resolution) {
            size_t numberOfPoints = resolution * resolution;
            size_t numberOfPointsPerSegment = numberOfPoints / numberOfSegments;
            size_t numberOfPointsPerLastSegment =
                    numberOfPointsPerSegment + numberOfPoints % numberOfSegments;
            types::real_t t;

            for (size_t segmentNumber = 1; segmentNumber < numberOfSegments;
                 segmentNumber++) {
                for (size_t i = 0; i < numberOfPointsPerSegment; i++) {
                    t = 1.0 * i / numberOfPointsPerSegment;
                    points.push_back(this->operator()(curve, t, segmentNumber));
                }
            }

            for (size_t i = 0; i < numberOfPointsPerLastSegment; i++) {
                t = 1.0 * i / numberOfPointsPerLastSegment;
                points.push_back(this->operator()(curve, t, numberOfSegments));
            }
        }

        P3CurvePlotter(const Curve &curve, size_t numberOfSegments)
                : P3CurvePlotter(curve, numberOfSegments, 80) {}

        explicit P3CurvePlotter(const Curve &curve) : P3CurvePlotter(curve,
                                                                     1,
                                                                     80) {}

        void Print(std::ostream &os, char curveChar, char backgroundChar) const {
            std::vector<std::vector<char>>
                    plot(resolution, std::vector<char>(resolution, backgroundChar));

            for (const auto &point : points) {
                int xIndex;
                int yIndex;
                if (resolution % 2 == 0) {
                    xIndex = static_cast<int>(floor(point.X * resolution / 2))
                             + static_cast<int>(resolution / 2);
                    yIndex = static_cast<int>(floor(point.Y * resolution / 2))
                             + static_cast<int>(resolution / 2);
                } else {
                    xIndex = static_cast<int>(round(point.X * resolution / 2))
                             + static_cast<int>((resolution - 1) / 2);
                    yIndex = static_cast<int>(round(point.Y * resolution / 2))
                             + static_cast<int>((resolution - 1) / 2);
                }
                if (xIndex >= 0 && static_cast<size_t>(xIndex) < resolution
                    && yIndex >= 0 && static_cast<size_t>(yIndex) < resolution) {
                    plot.at(static_cast<size_t>(yIndex)).at(static_cast<size_t>(xIndex)) =
                            curveChar;
                }
            }

            for (auto plotIterator = plot.rbegin(); plotIterator != plot.rend(); ++plotIterator) {
                std::vector<char> &line = *plotIterator;
                for (char character : line) {
                    os << character;
                }
                os << std::endl;
            }
        }

        void Print(std::ostream &os, char curveChar) const {
            this->Print(os, curveChar, ' ');
        }

        void Print(std::ostream &os) const {
            this->Print(os, '*', ' ');
        }

        void Print() const {
            this->Print(std::cout, '*', ' ');
        }

        Point2d operator()(const Curve &curve,
                           types::real_t t,
                           size_t segmentNumber) const {
            const Point2d &p0 =
                    curve(constants::NUM_OF_CUBIC_BEZIER_NODES * (segmentNumber - 1));
            const Point2d &p1 = curve(
                    constants::NUM_OF_CUBIC_BEZIER_NODES * (segmentNumber - 1) + 1);
            const Point2d &p2 = curve(
                    constants::NUM_OF_CUBIC_BEZIER_NODES * (segmentNumber - 1) + 2);
            const Point2d &p3 = curve(
                    constants::NUM_OF_CUBIC_BEZIER_NODES * (segmentNumber - 1) + 3);

            Point2d b0 = (1 - t) * p0 + t * p1;
            Point2d b1 = (1 - t) * p1 + t * p2;
            Point2d b2 = (1 - t) * p2 + t * p3;

            b0 = (1 - t) * b0 + t * b1;
            b1 = (1 - t) * b1 + t * b2;

            b0 = (1 - t) * b0 + t * b1;

            return b0;
        }
    private:
        std::vector<Point2d> points;
        size_t resolution;
    };

} // namespace bezier

#endif // BEZIER_H