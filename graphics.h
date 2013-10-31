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

struct Point2d {
	int x,y;
  Point2d(int X, int Y);
  Point2d();
};

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
  
  
public:
  
  XGraphics();
  ~XGraphics();
  int get_color(std::string c);
  void flush();
  void list_fonts();
  void setup_font();
  void erase_field();
  Point2d mouse_location();
  void draw_point(const Point2d& p, long col);
  void draw_line(const Point2d& p1, const Point2d& p2, long col);
  void draw_square(int x, int y, int z, long col);  
  void draw_rectangle(int x, int y, int zx, int zy, long col);
  void draw_filled_rectangle(int x, int y, int zx, int zy, long col);
  void draw_faint_line(const Point2d& p1, const Point2d& p2, long col);
  void erase_circle(const Point2d& p, int r);
  void draw_circle(const Point2d& p, int r, long col);
  void draw_concentric_circles(const Point2d& p, int r, long col);
  void draw_path(const std::vector<Point2d>& L, long col);
  void draw_text(const Point2d& p, std::stringstream &T, long col);
  void draw_text(const Point2d& p, std::string &S, long col);
  void draw_label(const Point2d& p, int i, long col);
  std::string wait_for_key();
};




#endif