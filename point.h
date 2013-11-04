#ifndef __POINT_H__
#define __POINT_H__

#include <iostream>

template <class T>
struct Point2d {
	T x,y;
  Point2d(T X, T Y);
  Point2d();
};

template <class T>
Point2d<T>::Point2d() {
  x = y = 0;
}

template <class T>
Point2d<T>::Point2d(T X, T Y) {
  x = X;
  y = Y;
}

template <class T>
Point2d<T> operator+(Point2d<T> a, Point2d<T> b) {
  return Point2d<T>(a.x+b.x, a.y+b.y);
}

template <class T>
Point2d<T> operator-(Point2d<T> a, Point2d<T> b) {
  return Point2d<T>(a.x-b.x, a.y-b.y);
}

template <class T>
Point2d<T> operator*(T sf, Point2d<T> a) {
  return Point2d<T>(sf*a.x, sf*a.y);
}

template <class T>
std::ostream& operator<<(std::ostream& os, Point2d<T>& p) {
  return os << "(" << p.x << "," << p.y << ")";
}


#endif