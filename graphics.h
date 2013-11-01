#ifndef __GRAPHICS_H__
#define __GRAPHICS_H__

extern "C" {
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
}

#include <string>
#include <sstream>
#include <map>

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
Point2d<T> operator+(Point2d<T>& a, Point2d<T>& b) {
  return Point2d<T>(a.x+b.x, a.y+b.y);
}

template <class T>
Point2d<T> operator-(Point2d<T>& a, Point2d<T>& b) {
  return Point2d<T>(a.x-b.x, a.y-b.y);
}

template <class T>
Point2d<T> operator*(T sf, Point2d<T> a) {
  return Point2d<T>(sf*a.x, sf*a.y);
}

class XGraphics {
private:
  Display *display;
  int screen_num;
  unsigned int display_width, display_height;
  XEvent report;
  GC gc;
  Window win;
  int border_width;
  unsigned int width, height;
  XFontStruct * font;
  Colormap screen_colormap;
  std::map<std::string, int> color_list;
  
  float scale;
  Point2d<float> translate;
  
public:
  
  XGraphics();
  XGraphics(int w, int h, float s, Point2d<float>& t);
  ~XGraphics();
  int get_color(std::string c);
  void flush();
  void list_fonts();
  void setup_font();
  void erase_field();
  Point2d<int> mouse_location();
  void draw_point(const Point2d<int>& p, long col);
  void draw_line(const Point2d<int>& p1, const Point2d<int>& p2, long col);
  void draw_line(const Point2d<float>& p1, const Point2d<float>& p2, long col);
  void draw_square(int x, int y, int z, long col);  
  void draw_rectangle(int x, int y, int zx, int zy, long col);
  void draw_filled_rectangle(int x, int y, int zx, int zy, long col);
  void draw_faint_line(const Point2d<int>& p1, const Point2d<int>& p2, long col);
  void erase_circle(const Point2d<int>& p, int r);
  void draw_circle(const Point2d<int>& p, int r, long col);
  void draw_concentric_circles(const Point2d<int>& p, int r, long col);
  void draw_path(const std::vector<Point2d<int> >& L, long col);
  void draw_text(const Point2d<int>& p, std::stringstream &T, long col);
  void draw_text(const Point2d<int>& p, std::string &S, long col);
  void draw_label(const Point2d<int>& p, int i, long col);
  std::string wait_for_key();
};




#endif