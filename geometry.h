#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>

const double eps = 0.000001;

class Line;

struct Point {
	double x = 0.0;
	double y = 0.0;
	Point() = default;
	Point(const double& x, const double& y) : x(x), y(y) {}

	friend std::ostream& operator<<(std::ostream& out, const Point& point);
	bool operator==(const Point& right) const {
		return (std::abs(x - right.x) < eps && std::abs(y - right.y) < eps);
	}

	bool operator!=(const Point& right) const {
		return !(*this == right);
	}

	Point& operator-=(const Point& right) {
		x -= right.x;
		y -= right.y;
		return *this;
	}

	Point& operator+=(const Point& right) {
		x += right.x;
		y += right.y;
		return *this;
	}


	void fromLines(const Line& first, const Line& second);

	double modul() const {
		return x * x + y * y;
	}
};

std::ostream& operator<<(std::ostream& out, const Point& point) {
	out << "(" << point.x << "," << point.y << ")";
	return out;
}

Point operator-(Point left, const Point& right) {
	return left -= right;
}

Point operator+(Point left, const Point& right) {
	return left += right;
}

double operator*(const Point& left, const Point& right) {
	return left.x * right.x + left.y * right.y;
}

struct Vector_;

class Line {
private:
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
public:
	friend std::ostream& operator<<(std::ostream&, const Line&);
	friend Point;
	Line() = default;
	Line(const double& k, const double& b) : a(k), b(-1), c(b) {}
	Line(const Point& point, const double& k) : a(k), b(-1), c(point.y - point.x * k) {}
	Line(const Point& point1, const Point& point2) : a(point1.y - point2.y), b(point2.x - point1.x), c(point1.x* point2.y - point1.y * point2.x) {}
	Line(double a, double b, double c) : a(a), b(b), c(c) {}
	double yVal(const double& x) const { 
		return (a * x + c) / (-b);
	}

	bool operator==(const Line& right) const {
		double ta = 0.0;
		if (abs(right.a) < eps) {
			if (abs(a) > eps) return false;
			else ta = 1;
		}
		else ta = a / right.a;
		double tb = 0.0;
		if (abs(right.b) < eps) {
			if (abs(b) > eps) return false;
			else tb = 1;
		}
		else tb = b / right.b;
		if (abs(ta - tb) > eps) return false;
		double tc = 0.0;
		if (abs(right.c) < eps) {
			if (abs(c) > eps) return false;
			else tc = 1;
		}
		else tc = c / right.c;
		if (abs(ta - tc) > eps) return false;
		return true;
	}

	Vector_ normal() const;

	bool operator!=(const Line& right) const {
		return !(*this == right);
	}
};

std::ostream& operator<<(std::ostream& out, const Line& line) {
	out << "(" << line.a << "," << line.b << "," << line.c << ")";
	return out;
}

void Point::fromLines(const Line& first, const Line& second) {
	x = (first.b * second.c - second.b * first.c) / (first.a * second.b - first.b * second.a);
	y = (first.c * second.a - second.c * first.a) / (first.a * second.b - first.b * second.a);
}

struct Vector_ {
	Point end;
	Vector_() = default;
	Vector_(const Point& begin, const Point& end) : end((end.x - begin.x), (end.y - begin.y)) {}

	double x() const {
		return end.x;
	}

	double y() const {
		return end.y;
	}

	int signOfRotation(const Vector_& right) const { //������ ������� +
		int sign = 0;
		if (x() == 0) {
			if (end.y > 0) sign = (right.x() > 0 ? -1 : 1);
			else sign = (right.x() < 0 ? -1 : 1);
		}
		else {
			Line prev = Line(Point(), end);
			if (x() > 0) {
				if (prev.yVal(right.x()) < right.y()) sign = 1;
				else sign = -1;
			}
			else {
				if (prev.yVal(right.x()) > right.y()) sign = 1;
				else sign = -1;
			}
		}
		return sign;
	}

	double length() const {
		double len = std::sqrt(end.x * end.x + end.y * end.y);
		return len;
	}

	Vector_& operator*=(double k) {
		end.x *= k;
		end.y *= k;
		return *this;
	}
};

Vector_ operator*(Vector_ v, double k) {
	return v *= k;
}

Vector_ Line::normal() const {
	Vector_ to_return = Vector_(Point(0, 0), Point(a, b));
	return to_return;
}

class Shape {
protected:
	mutable double perimetr = 0.0;
	mutable double areaa = 0.0;

public:
	virtual bool operator==(const Shape&) const { return false; }
	virtual bool operator!=(const Shape&) const {  return true; }
	virtual bool isCongruentTo(const Shape&) const { return false; }
	virtual bool isSimilarTo(const Shape&) const {  return false;  }
	virtual double area(bool flag = false) const = 0;
	virtual double perimeter(bool flag = false) const = 0;
	virtual bool containsPoint(const Point& point) const = 0;
	virtual void rotate(const Point& center, const double& angle) = 0;
	virtual void reflex(const Point&) = 0;
	virtual void reflex(const Line&) = 0;
	virtual void scale(const Point& center, const double& coefficient) = 0;

	virtual ~Shape() = 0;
};

Shape::~Shape() = default;

class Polygon : public Shape {
protected:
	mutable std::vector<Vector_> sides; //������ �������� ������ 0->1, 1->2, 2->3 ... n->0
	mutable std::vector<double> sides_length; //������ ���� ������ 0->1, 1->2, 2->3 ... n->0
	std::vector<Point> vertices;
public:
	Polygon() = default;
	Polygon(std::vector<Point> points) : vertices(points) {}
	Polygon(Point first, std::vector<Vector_> sides) : sides(sides) {
		vertices.resize(sides.size());
		vertices[0] = first;
		for (size_t i = 1; i < sides.size(); ++i) {
			vertices[i] = sides[i - 1].end + vertices[i - 1];
		}
	}
	Polygon(std::initializer_list<Point> points) {
		vertices.resize(points.size());
		std::copy(points.begin(), points.end(), vertices.begin());
	}

	const std::vector<Point> getVertices() const { return vertices; }
	
	virtual long long verticesCount() const {
		return vertices.size();
	}

	void setSidesVectors(bool flag = false) const { 
		if (sides.size() == 0 || (sides.size() != 0 && flag == true)) {
			sides.resize(vertices.size());
			for (size_t i = 0; i + 1 < sides.size(); ++i) {
				sides[i] = Vector_(vertices[i], vertices[i + 1]);
			}
			sides[sides.size() - 1] = Vector_(vertices[sides.size() - 1], vertices[0]);
		}
	}

	void setSidesLength(bool flag = false) const {
		if (sides_length.size() == 0 || (sides_length.size() != 0 && flag == true)) {
			setSidesVectors(flag);
			sides_length.resize(vertices.size());
			for (size_t i = 0; i < sides.size(); ++i) {
				sides_length[i] = sides[i].length();
			}
		}
	}

	double perimeter(bool flag = false) const override {
		if (perimetr == 0 || (perimetr != 0 && flag == true)) {
			if (areaa == 0) setSidesLength(flag);
			for (size_t i = 0; i < sides.size(); ++i) {
				perimetr += sides_length[i];
			}
		}
		return perimetr;
	}

	virtual bool isConvex() const {
		bool to_return = true;
		setSidesVectors();
		int sign = sides[0].signOfRotation(sides[1]);
		size_t inc = 1;
		while (inc < sides.size() && to_return) {
			int compare = sides[inc].signOfRotation(sides[(inc + 1) % sides.size()]);
			if (compare != sign) to_return = false;
			++inc;
		}
		return to_return;
	}

	std::vector<double> setAngles() const {
		std::vector<double> angles(vertices.size()); 
		setSidesLength();
		for (size_t i = 0; i < angles.size(); ++i) {
			angles[i] = (sides[(i + 1) % angles.size()].x() * sides[i].x() + sides[(i + 1) % angles.size()].y() * sides[i].y()) / (sides_length[(i + 1) % sides_length.size()] * sides_length[i]);
		}
		return angles;
	}

	bool isCongruentTo(const Polygon& another) const {
		if (vertices.size() != another.vertices.size()) return false;
		if (this == &another) return true;
		setSidesLength();
		another.setSidesLength();
		size_t inc = 0;
		while (inc < sides_length.size()) {
			while (inc < sides_length.size() && abs(sides_length[inc] - another.sides_length[0]) > eps) {
				++inc;
			}
			if (inc == sides_length.size()) return false;
			bool flag = true;
			size_t tmp = 1;                        
			if (abs(sides_length[(inc + 1) % sides_length.size()] - another.sides_length[tmp]) < eps) {
				while (flag && tmp < sides_length.size()) {
					if (abs(sides_length[(inc + tmp) % sides_length.size()] - another.sides_length[tmp]) > eps) flag = false;
					++tmp;
				}
			}
			else if (abs(sides_length[(sides_length.size() + inc - 1) % sides_length.size()] - another.sides_length[tmp]) < eps) {
				while (flag && tmp < sides_length.size()) {
					if (abs(sides_length[(sides_length.size() + inc - tmp) % sides_length.size()] - another.sides_length[tmp]) > eps) flag = false;
					++tmp;
				}
			}
			else {
				++inc;
				continue;
			}
			if (flag == false) {
				++inc;
				continue;
			}
			tmp = 1;
			std::vector<double> angles = setAngles();
			std::vector<double> another_angles = another.setAngles();
			if (abs(sides_length[(inc + 1) % sides_length.size()] - another.sides_length[tmp]) < eps) {
				while (flag && tmp < sides_length.size()) {
					if (abs(angles[(inc + tmp) % sides_length.size()] - another_angles[tmp]) > eps) flag = false;
					++tmp;
				}
			}
			else if (abs(sides_length[(sides_length.size() + inc - 1) % sides_length.size()] - another.sides_length[tmp]) < eps) {
				while (flag && tmp < sides_length.size()) {
					if (abs(angles[(sides_length.size() + inc - tmp - 1) % sides_length.size()] - another_angles[tmp]) > eps) flag = false;
					++tmp;
				}
			}
			if (flag == true) {
				return flag;
			}
			++inc;
		}
		return false;
	}

	bool isSimilarTo(const Polygon& another) const {
		std::vector<double> angles = setAngles();
		std::vector<double> another_angles = another.setAngles();
		size_t inc = 0;
		while (inc < sides_length.size()) {
			while (inc < sides_length.size() && abs(angles[inc] - another_angles[0]) > eps) {
				++inc;
			}
			if (inc == sides_length.size()) return false;
			bool flag = true;
			size_t tmp = 1;                  
			if (abs(angles[(inc + 1) % sides_length.size()] - another_angles[tmp]) < eps) {
				while (flag && tmp < sides_length.size()) {
					if (abs(angles[(inc + tmp) % sides_length.size()] - another_angles[tmp]) > eps) flag = false;
					++tmp;
				}
			}
			else if (abs(angles[(sides_length.size() + inc - 1) % sides_length.size()] - another_angles[tmp]) < eps) {
				while (flag && tmp < sides_length.size()) {
					if (abs(angles[(sides_length.size() + inc - tmp) % sides_length.size()] - another_angles[tmp]) > eps) flag = false;
					++tmp;
				}
			}
			else {
				++inc;
				continue;
			}
			if (flag == false) {
				++inc;
				continue;
			}
			tmp = 1;
			double k = sides_length[inc] / another.sides_length[0];
			if (abs(k - sides_length[(inc + 1) % sides_length.size()] / another.sides_length[tmp]) < eps) {
				while (flag && tmp < sides_length.size()) {
					if (abs(k - sides_length[(inc + tmp) % sides_length.size()] / another.sides_length[tmp]) > eps) flag = false;
					++tmp;
				}
			}
			else if (abs(sides_length[(sides_length.size() + inc - 1) % sides_length.size()] - another.sides_length[tmp]) < eps) {
				while (flag && tmp < sides_length.size()) {
					if (abs(k - sides_length[(sides_length.size() + inc - tmp) % sides_length.size()] / another.sides_length[tmp]) > eps) flag = false;
					++tmp;
				}
			}
			if (flag == true) {
				return flag;
			}
			++inc;
		}
		return false;
	}

	double area(bool) const override;

	double convexArea(bool) const;

	bool operator==(const Polygon& another) const {
		if (this == &another) return true;
		if (vertices.size() != another.vertices.size()) return false;
		size_t inc = 0;
		while (inc < vertices.size() && vertices[inc] != another.vertices[0]) {
			++inc;
		}
		if (inc == vertices.size()) return false;
		bool flag = true;
		size_t tmp = 1;                       
		if (vertices[(inc + 1) % vertices.size()] == another.vertices[tmp]) {
			while (flag && tmp < vertices.size()) {
				if (vertices[(inc + tmp) % vertices.size()] != another.vertices[tmp]) flag = false;
				++tmp;
			}
		}
		else if (vertices[(vertices.size() + inc - 1) % vertices.size()] == another.vertices[tmp]) {
			while (flag && tmp < vertices.size()) {
				if (vertices[(vertices.size() + inc - tmp) % vertices.size()] != another.vertices[tmp]) flag = false;
				++tmp;
			}
		}
		else return false;
		return flag;
	}

	bool operator!=(const Polygon& another) const { return !(*this == another); }
	
	bool operator==(const Shape& another) const override {
		if (dynamic_cast<Polygon*>(const_cast<Shape*>(&another)) == nullptr) return false;
		return (*this == *dynamic_cast<Polygon*>(const_cast<Shape*>(&another)));
	}

	bool operator!=(const Shape& another) const override {
		return !(*this == another);
	}

	bool isSimilarTo(const Shape& another) const override {
		if (dynamic_cast<Polygon*>(const_cast<Shape*>(&another)) == nullptr) return false;
		return (this->isSimilarTo(*dynamic_cast<Polygon*>(const_cast<Shape*>(&another))));
	}

	bool isCongruentTo(const Shape& another) const override {
		if (dynamic_cast<Polygon*>(const_cast<Shape*>(&another)) == nullptr) return false;
		return (this->isCongruentTo(*dynamic_cast<Polygon*>(const_cast<Shape*>(&another))));
	}

	void rotate(const Point& center, const double& angle) override {
		for (size_t i = 0; i < vertices.size(); ++i) {
			Point new_vertice = Point((vertices[i].x - center.x) * cos(angle * M_PI / 180) - (vertices[i].y - center.y) * sin(angle * M_PI / 180) + center.x, (vertices[i].x - center.x) * sin(angle * M_PI / 180) + (vertices[i].y - center.y) * cos(angle * M_PI / 180) + center.y);
			vertices[i] = new_vertice;
		}
		setSidesVectors(true);
		area(true);
	}

	void scale(const Point& center, const double& coefficient) override {
		for (size_t i = 0; i < vertices.size(); ++i) {
			vertices[i] = (Vector_(center, vertices[i]) * coefficient).end + center;
		}
		area(true);
		perimeter(true);
	}
	
	bool containsInConvex(const Point& point) const {
		std::vector<Vector_> vectors_to_virtices(vertices.size() + 1);
		for (size_t i = 0; i < vertices.size() + 1; ++i) {
			vectors_to_virtices[i] = Vector_(point, vertices[i % vertices.size()]);
		}
		size_t inc = 1;
		int sign = vectors_to_virtices[0].signOfRotation(vectors_to_virtices[1]);
		while (inc  < vertices.size()) {
			int new_sign = vectors_to_virtices[inc].signOfRotation(vectors_to_virtices[inc + 1]);
			if (new_sign * sign < -eps) return false;
			++inc;
		}
		return true;
	}
	bool containsPoint(const Point& point) const override;

	void reflex(const Point& center) override {
		for (size_t i = 0; i < vertices.size(); ++i) {
			vertices[i] = (Vector_(vertices[i], center) * 2).end + vertices[i];
		}
		setSidesVectors(true);
	}

	void reflex(const Line& axis) override {
		for (size_t i = 0; i < vertices.size(); ++i) {
			Point point;
			point.fromLines(axis, Line(vertices[i], vertices[i] + axis.normal().end));
			vertices[i] = (Vector_(vertices[i], point) * 2).end;
		}
		setSidesVectors();
	}

	~Polygon() = default;
};

class Ellipse : public Shape {
protected:
	Point f1;
	Point f2;
	double a = 0.0;
	double c = 0.0;
	double b = 0.0;
public:
	Ellipse() = default;
	Ellipse(const Point& f1, const Point& f2, double delta) : f1(f1), f2(f2), a(delta / 2), c(std::sqrt((f1.x - f2.x)* (f1.x - f2.x) + (f1.y - f2.y) * (f1.y - f2.y)) / 2), b(std::sqrt(a* a - c * c)) {}

	std::pair<Point, Point> focuses() const {
		std::pair<Point, Point> to_return;
		to_return.first = f1;
		to_return.second = f2;
		return to_return;
	}

	virtual Point center() const {
		Point to_return = Point((f1.x + f2.x) / 2, (f1.y + f2.y) / 2);
		return to_return;
	}

	virtual double eccentricity() const {
		double to_return = c / a;
		return to_return;
	}

	std::pair<Line, Line> directrices() const {
		std::pair<Line, Line> to_return;
		Point centr = center();
		to_return.first = Line(Point(centr.x - a * a / c, centr.y), Point(centr.x - a * a / c, centr.y + 1));
		to_return.first = Line(Point(centr.x + a * a / c, centr.y), Point(centr.x + a * a / c, centr.y + 1));
		return to_return;
	}

	double area(bool flag = false) const override {
		if (areaa == 0 || (areaa != 0 && flag == true)) {
			areaa = M_PI * a * b;
		}
		return areaa;
	}

	double perimeter(bool flag = false) const override {
		if (perimetr == 0 || (perimetr != 0 && flag == true)) {
			perimetr = M_PI * (3 * (a + b) - std::sqrt((3 * a + b) * (a + 3 * b)));
		}
		return perimetr;
	}

	bool isCongruentTo(const Ellipse& another) const {
		return std::abs(a - another.a) < eps && std::abs(Vector_(f1, f2).length() - Vector_(another.f1, another.f2).length()) < eps;
	}

	bool isSimilarTo(const Ellipse& another) const {
		return std::abs(a / another.a - c / another.c) < eps;
	}
	
	bool isSimilarTo(const Shape& another) const override {
		if (dynamic_cast<Ellipse*>(const_cast<Shape*>(&another)) == nullptr) return false;
		return (this->isSimilarTo(*dynamic_cast<Ellipse*>(const_cast<Shape*>(&another))));
	}

	bool isCongruentTo(const Shape& another) const override {
		if (dynamic_cast<Ellipse*>(const_cast<Shape*>(&another)) == nullptr) return false;
		return (this->isCongruentTo(*dynamic_cast<Ellipse*>(const_cast<Shape*>(&another))));
	}

	bool containsPoint(const Point& point) const override {
		bool flag = (Vector_(point, f1).length() + Vector_(point, f2).length() - 2 * a < eps);
		return flag;
	}

	void rotate(const Point& center, const double& angle) override {
		Point new_focus = Point((f1.x - center.x) * cos(angle * M_PI / 180) - (f1.y - center.y) * sin(angle * M_PI / 180) + center.x, (f1.x - center.x) * sin(angle * M_PI / 180) + (f1.y - center.y) * cos(angle * M_PI / 180) + center.y);
		f1 = new_focus;
		new_focus = Point((f2.x - center.x) * cos(angle * M_PI / 180) - (f2.y - center.y) * sin(angle * M_PI / 180) + center.x, (f2.x - center.x) * sin(angle * M_PI / 180) + (f2.y - center.y) * cos(angle * M_PI / 180) + center.y);
		f2 = new_focus;
	}
	
	void scale(const Point& center, const double& coefficient) override {
		f1 = (Vector_(center, f1) * coefficient).end + center;
		f2 = (Vector_(center, f2) * coefficient).end + center;
		a *= coefficient;
		c *= coefficient;
		b *= coefficient;
		perimeter(true);
		area(true);
	}

	bool operator==(const Ellipse& another) {
		bool flag = f1 == another.f1 && f2 == another.f2 && std::abs(a - another.a) < eps;
		return flag;
	}

	bool operator!=(const Ellipse& another) {
		return !(*this == another);
	}

	bool operator==(const Shape& another) const override {
		if (dynamic_cast<Ellipse*>(const_cast<Shape*>(&another)) == nullptr) return false;
		return (*this == *dynamic_cast<Ellipse*>(const_cast<Shape*>(&another)));
	}

	bool operator!=(const Shape& another) const override {
		return !(*this == another);
	}

	void reflex(const Point& center) override {
		f1 = (Vector_(f1, center) * 2).end + f1;
		f2 = (Vector_(f2, center) * 2).end + f2;
	}

	void reflex(const Line& axis) override {
		Point point;
		point.fromLines(axis, Line(f1, f1 + axis.normal().end));
		f1 = (Vector_(f1, point) * 2).end;
		point.fromLines(axis, Line(f2, f2 + axis.normal().end));
		f2 = (Vector_(f2, point) * 2).end;
	}

	~Ellipse() = default;
};

class Circle : public Ellipse {
public:
	Circle() = default;
	Circle(const Point& center, const double& radius) : Ellipse(center, center, 2 * radius) {}

	double radius() const {
		return a;
	}

	double perimeter(bool flag = false) const override {
		if (perimetr == 0 || (perimetr != 0 && flag == true)) {
			perimetr = 2 * M_PI * a;
		}
		return perimetr;
	}

	double area(bool flag = false) const override {
		if (areaa == 0 || (areaa != 0 && flag == true)) {
			areaa = M_PI * a * a;
		}
		return areaa;
	}

	Point center() const override { return f1; }
	double eccentricity() const override { return 0; }
	~Circle() = default;
};

class Rectangle : public Polygon {
public:
	using Polygon::Polygon;
	Rectangle() = default;
	Rectangle(Point p1, Point p2, Point p3, Point p4) : Polygon({ p1, p2, p3, p4 }) {}

	Point center() const {
		return Point((vertices[0].x + vertices[2].x) / 2, (vertices[0].y + vertices[2].y) / 2);
	}

	Rectangle(const Point& first, const Point& second, double coefficient) { 
		vertices.resize(4);
		vertices[0] = first;
		vertices[2] = second;
		double cosin = (coefficient * coefficient - 1) / (coefficient * coefficient + 1);
		double sinus = 2 * coefficient / (coefficient * coefficient + 1);
		Point centr = center();
		if (coefficient > 1) {
			cosin *= -1;
		}
		vertices[1] = Point((second.x - centr.x) * cosin - (second.y - centr.y) * sinus + centr.x, (second.x - centr.x) * sinus + (second.y - centr.y) * cosin + centr.y);
		vertices[3] = Point((first.x - centr.x) * cosin - (first.y - centr.y) * sinus + centr.x, (first.x - centr.x) * sinus + (first.y - centr.y) * cosin + centr.y);
	}

	double area(bool flag = false) const override {
		if (areaa == 0 || (areaa != 0 && flag == true)) {
			setSidesLength(flag);
			areaa = sides_length[0] * sides_length[1];
		}
		return areaa;
	}

	bool isCongruentTo(const Rectangle& another) const {
		if (this == &another) return true;
		setSidesLength();
		another.setSidesLength();
		bool flag =(std::abs(sides_length[0] - another.sides_length[0]) < eps && std::abs(sides_length[1] - another.sides_length[1]) < eps) || (std::abs(sides_length[1] - another.sides_length[0]) < eps && std::abs(sides_length[0] - another.sides_length[1]) < eps);
		return flag;
	}

	long long verticesCount() const override { return 4; }
	bool isConvex() const override { return true; }

	~Rectangle() = default;
};

class Square : public Rectangle {
public:
	using Rectangle::Rectangle;

	Square(const Point& first, const Point& second) {
		vertices.resize(4);
		vertices[0] = first;
		vertices[2] = second;
		vertices[1] = Point((vertices[2].x - center().x) * cos(M_PI_2) - (vertices[2].y - center().y) * sin(M_PI_2) + center().x, (vertices[2].x - center().x) * sin(M_PI_2) + (vertices[2].y - center().y) * cos(M_PI_2) + center().y);
		vertices[3] = Point((vertices[0].x - center().x) * cos(M_PI_2) - (vertices[0].y - center().y) * sin(M_PI_2) + center().x, (vertices[0].x - center().x) * sin(M_PI_2) + (vertices[0].y - center().y) * cos(M_PI_2) + center().y);
	}

	Circle circumscribedCircle() const {
		Circle to_return = Circle(center(), std::sqrt((center().x - vertices[0].x) * (center().x - vertices[0].x) + (center().y - vertices[0].y) * (center().y - vertices[0].y)));
		return to_return;
	}

	Circle inscribedCircle() const {
		Circle to_return = Circle(center(), std::sqrt(((center().x - vertices[0].x) * (center().x - vertices[0].x) + (center().y - vertices[0].y) * (center().y - vertices[0].y)) / 2));
		return to_return;
	}

	bool isCongruentTo(const Square& another) const {
		setSidesLength();
		another.setSidesLength();
		bool flag = (std::abs(sides_length[0] - another.sides_length[0]) < eps);
		return flag;
	}

	~Square() = default;
};

class Triangle : public Polygon {
private:
	mutable Circle circumscribed;
	mutable Circle inscribed;

public:
	using Polygon::Polygon;

	Triangle(Point p1, Point p2, Point p3) : Polygon({p1, p2, p3}) {}

	Circle circumscribedCircle(bool flag = false) const {
		if (circumscribed == Circle() || (circumscribed != Circle() && flag == true)) {
			double x12 = vertices[0].x - vertices[1].x;
			double x23 = vertices[1].x - vertices[2].x;
			double x31 = vertices[2].x - vertices[0].x;

			double y12 = vertices[0].y - vertices[1].y;
			double y23 = vertices[1].y - vertices[2].y;
			double y31 = vertices[2].y - vertices[0].y;

			double z1 = vertices[0].modul();
			double z2 = vertices[1].modul();
			double z3 = vertices[2].modul();

			Point to_return_center = Point(-(y12 * z3 + y23 * z1 + y31 * z2) / (2 * (x12 * y31 - y12 * x31)), (x12 * z3 + x23 * z1 + x31 * z2) / (2 * (x12 * y31 - y12 * x31)));
			circumscribed = Circle(to_return_center, Vector_(to_return_center, vertices[0]).length());
		}
		circumscribed.radius();
		return circumscribed;
	}

	double area(bool flag = false) const override {
		if (areaa == 0 || (areaa != 0 && flag == true)) {
			setSidesLength();
			double p = perimeter() / 2;
			areaa = std::sqrt(p * (p - sides_length[0]) * (p - sides_length[1]) * (p - sides_length[2]));
		}	
		return areaa;
	}

	Circle inscribedCircle(bool flag = false) const {
		setSidesLength();
		if (inscribed == Circle() || (inscribed != Circle() && flag == true)) {
			double k = sides_length[2] / sides_length[0];
			double x11 = (vertices[2].x + k * vertices[1].x) / (1 + k);
			double y11 = (vertices[2].y + k * vertices[1].y) / (1 + k);
			k = sides_length[1] / sides_length[0];
			double x12 = (vertices[2].x + k * vertices[0].x) / (1 + k);
			double y12 = (vertices[2].y + k * vertices[0].y) / (1 + k);
			Point to_return_center;
			if (vertices[0].x == x11) {
				to_return_center.x = vertices[0].x;
				double k2 = (y12 - vertices[1].y) / (x12 - vertices[1].x);
				double b2 = (x12 * vertices[1].y - vertices[1].x * y12) / (x12 - vertices[1].x);
				to_return_center.y = k2 * to_return_center.x + b2;
			}
			else if (vertices[1].x == x12) {
				to_return_center.x = vertices[1].x;
				double k1 = (y11 - vertices[0].y) / (x11 - vertices[0].x);
				double b1 = (x11 * vertices[0].y - vertices[0].x * y11) / (x11 - vertices[0].x);
				to_return_center.y = k1 * to_return_center.x + b1;
			}
			else {
				double k1 = (y11 - vertices[0].y) / (x11 - vertices[0].x);
				double b1 = (x11 * vertices[0].y - vertices[0].x * y11) / (x11 - vertices[0].x);
				double k2 = (y12 - vertices[1].y) / (x12 - vertices[1].x);
				double b2 = (x12 * vertices[1].y - vertices[1].x * y12) / (x12 - vertices[1].x);
				to_return_center.x = (b1 - b2) / (k2 - k1);
				to_return_center.y = (k2 * b1 - k1 * b2) / (k2 - k1);
			}
			inscribed = Circle(to_return_center, 2 * area() / perimeter());
		}

		inscribed.radius();
		return inscribed;
	}

	Point centroid() const {
		Point to_return = Point((vertices[0].x + vertices[1].x + vertices[2].x) / 3, (vertices[0].y + vertices[1].y + vertices[2].y) / 3);
		return to_return;
	}

	Point orthocenter() const {
		Line a = Line(vertices[0], vertices[1]);
		Line b = Line(vertices[0], vertices[2]);
		Line h1 = Line(vertices[2], vertices[2] + a.normal().end);
		Line h2 = Line(vertices[1], vertices[1] + b.normal().end);
		Point to_return = Point();
		to_return.fromLines(h1, h2);
		return to_return;
	}

	Line EulerLine() const {
		Line to_return = Line(orthocenter(), centroid());
		return to_return;
	}

	Circle ninePointsCircle() const {
		Point ortho = orthocenter();
		Point to_return_center = Point((ortho.x + circumscribedCircle().center().x) / 2, (ortho.y + circumscribedCircle().center().y) / 2);
		Circle to_return = Circle(to_return_center, circumscribedCircle().radius() / 2);
		return to_return;
	}

	void rotate(const Point& center, const double& angle) override {
		for (size_t i = 0; i < vertices.size(); ++i) {
			Point new_vertice = Point((vertices[i].x - center.x) * cos(angle * M_PI / 180) - (vertices[i].y - center.y) * sin(angle * M_PI / 180) + center.x, (vertices[i].x - center.x) * sin(angle * M_PI / 180) + (vertices[i].y - center.y) * cos(angle * M_PI / 180) + center.y);
			vertices[i] = new_vertice;
		}
		setSidesVectors(true); 
		circumscribedCircle(true);
		inscribedCircle(true);
	}

	void scale(const Point& center, const double& coefficient) override {
		for (size_t i = 0; i < vertices.size(); ++i) {
			vertices[i] = (Vector_(center, vertices[i]) * coefficient).end + center;
		}
		area(true);
		perimeter(true);
		circumscribedCircle(true);
		inscribedCircle(true);
	}

	~Triangle() = default;
};

double Polygon::convexArea(bool flag = false) const {
	if (areaa == 0 || (areaa != 0 && flag == true)) {
		areaa = 0;
		setSidesLength(flag);
		for (size_t i = 2; i < vertices.size(); ++i) {
			areaa += Triangle{vertices[0], vertices[i - 1], vertices[i]}.area();
		}
	}
	return areaa;
}

double Polygon::area(bool flag = false) const {
	if (isConvex()) return convexArea(flag);
	if (areaa == 0 || (areaa != 0 && flag)) {
		areaa = 0;
		setSidesLength(flag);
		std::vector<Vector_> convex_polygon;
		std::vector<Point> convex_points;
		std::vector<int> signs(vertices.size());
		std::vector<Polygon> polygons;
		int plus = 0;
		for (size_t i = 0; i < vertices.size(); ++i) {
			signs[i] = sides[i].signOfRotation(sides[(i + 1) % sides.size()]);
			if (signs[i] > 0) ++plus;
			else --plus;
		}
		if (plus == 0) std::cerr << "asssssssssHOLEE";
		convex_points.push_back(vertices[0]);
		for (size_t i = 0; i < vertices.size(); ++i) {
			convex_polygon.push_back(sides[i]);
			if (i + 1 < vertices.size()) convex_points.push_back(vertices[(i + 1) % vertices.size()]);
			if (signs[i] * plus < 0) {
				polygons.push_back(Triangle{ vertices[i], vertices[(i + 1) % vertices.size()], vertices[(i + 2) % vertices.size()] });
				convex_polygon.pop_back();
				convex_points.pop_back();
				sides[(i + 1) % vertices.size()] = Vector_(convex_points[convex_points.size() - 1], vertices[(i + 2) % vertices.size()]);
				signs[(i + 1) % vertices.size()] = sides[(i + 1) % vertices.size()].signOfRotation(sides[(i + 2) % vertices.size()]);
			}
		}
		Polygon convex = Polygon(convex_points[0], convex_polygon);
		areaa += convex.convexArea();
		for (size_t i = 0; i < polygons.size(); ++i) {
			areaa -= polygons[i].area();
		}
	}
	setSidesVectors(true);
	return areaa;
}

bool Polygon::containsPoint(const Point& point) const {
	if (isConvex()) return containsInConvex(point);
	setSidesVectors();
	std::vector<Vector_> convex_polygon;
	std::vector<Point> convex_points;
	std::vector<int> signs(vertices.size());
	std::vector<Polygon> polygons;
	int plus = 0;
	for (size_t i = 0; i < vertices.size(); ++i) {
		signs[i] = sides[i].signOfRotation(sides[(i + 1) % sides.size()]);
		if (signs[i] > 0) ++plus;
		else --plus;
	}
	if (plus == 0) std::cerr << "asssssssssHOLEE";
	convex_points.push_back(vertices[0]);
	for (size_t i = 0; i < vertices.size(); ++i) {
		convex_polygon.push_back(sides[i]);
		if (i + 1 < vertices.size()) convex_points.push_back(vertices[(i + 1) % vertices.size()]);
		if (signs[i] * plus < 0) {
			polygons.push_back(Triangle(vertices[i], vertices[(i + 1) % vertices.size()], vertices[(i + 2) % vertices.size()]));
			convex_polygon.pop_back();
			convex_points.pop_back();
			sides[(i + 1) % vertices.size()] = Vector_(convex_points[convex_points.size() - 1], vertices[(i + 2) % vertices.size()]);
			signs[(i + 1) % vertices.size()] = sides[(i + 1) % vertices.size()].signOfRotation(sides[(i + 2) % vertices.size()]);
		}
	}
	Polygon convex = Polygon(convex_points[0], convex_polygon);
	if (convex.containsInConvex(point) == false) return false;
	for (size_t i = 0; i < polygons.size(); ++i) {
		if (polygons[i].containsInConvex(point) == true) return false;
	}
	setSidesVectors(true);
	return true;
}
