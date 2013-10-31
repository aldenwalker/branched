#include <vector>
#include <iostream>
#include <sstream>

#include "graphics.h"

Point2d::Point2d() {
  x = y = 0;
}

Point2d::Point2d(int X, int Y) {
  x = X;
  y = Y;
}


XGraphics::XGraphics() {
  border_width = 4;
  display=XOpenDisplay(NULL);
  screen_num = DefaultScreen(display);  
  display_width = DisplayWidth(display, screen_num);
  display_height = DisplayHeight(display, screen_num);
  width = 800;
  height = 800;
  win = XCreateSimpleWindow(display, RootWindow(display, screen_num), 0, 0, width, 
                            height, border_width, BlackPixel(display, screen_num), WhitePixel(display, screen_num));
  XSelectInput(display, win, ExposureMask | 
                             KeyPressMask | 
                             ButtonPressMask | 
                             PointerMotionMask |
                             StructureNotifyMask);
  gc = DefaultGC(display, screen_num);
  screen_colormap = DefaultColormap(display, screen_num);
  XSetForeground(display, gc, BlackPixel(display, screen_num));
  XSetBackground(display, gc, WhitePixel(display, screen_num));
  XMapWindow(display, win);
  color_list.clear();
  
  while (true) {  //wait until the window is actually mapped
    XEvent e;
    XNextEvent(display, &e);
    if (e.type == Expose) break;
  }
  setup_font();
}

XGraphics::~XGraphics() {
  XCloseDisplay(display);
}


int XGraphics::get_color(std::string c) {
  std::map<std::string, int>::iterator it;
  it = color_list.find(c);
  if (it == color_list.end()) {
    XColor temp;
    Status rc = XAllocNamedColor(display, screen_colormap, c.c_str(), &temp, &temp);
    if (rc == 0) {
      std::cout << "Color error\n";
      return 0;
    }
    color_list[c] = temp.pixel;
  }
  return color_list[c];
}
    
    


void XGraphics::flush() {
  XFlush(display);
}

void XGraphics::list_fonts() {
  std::cout << "Listing fonts: \n";
  int num_returned;
  char** names = XListFonts(display, "*", 10000, &num_returned);
  for (int i=0; i<num_returned; ++i) {
    std::cout << names[i] << "\n";
  }
}


void XGraphics::setup_font(void){
    const char * fontname = "-*-georgia-*-r-*-*-14-*-*-*-*-*-*-*";
 //   const char * fontname = "-*-times-*-r-*-*-16-*-*-*-*-*-*-*";

    font = XLoadQueryFont (display, fontname);
    /* If the font could not be loaded, revert to the "fixed" font. */
    if (! font) {
      font = XLoadQueryFont (display, "fixed");
      std::cout << "couldn't find font!\n";
    }
    XSetFont (display, gc, font->fid);
}

	
void XGraphics::erase_field(void){
	XClearWindow(display, win);
}

Point2d XGraphics::mouse_location(){
//    Bool result;
  Window window_returned;
  int root_x, root_y;
  int win_x, win_y;
  unsigned int mask_return;
  Point2d p;
    
	XQueryPointer(display, win, &window_returned,
                &window_returned, &root_x, &root_y, &win_x, &win_y,
                &mask_return);
  p.x=win_x;
  p.y=win_y;
  return p;
}

void XGraphics::draw_point(const Point2d& p, long col){
  XSetForeground(display, gc, col);
  XDrawPoint(display, win, gc, p.x, height-p.y);
}

void XGraphics::draw_line(const Point2d& p1, const Point2d& p2, long col) {
  XSetForeground(display, gc, col);
  XSetLineAttributes(display, gc, 1.5, LineSolid, 1, 1);
  XDrawLine(display, win, gc, p1.x, height-p1.y, p2.x, height-p2.y);
}

void XGraphics::draw_square(int x, int y, int z, long col){
  Point2d p1(x,y);
  Point2d p2(x+z,y);
  Point2d p3(x+z,y+z);
  Point2d p4(x,y+z);
  this->draw_line(p1, p2, col);
  this->draw_line(p2, p3, col);
  this->draw_line(p3, p4, col);
  this->draw_line(p4, p1, col);
  //XSetForeground(display, gc, col);
  //XSetLineAttributes(display, gc, 2, LineSolid, 1, 1);
  //XDrawLine(display, win, gc, x, y, x+z, y);
  //XDrawLine(display, win, gc, x+z, y, x+z, y+z);
  //XDrawLine(display, win, gc, x+z, y+z, x, y+z);
  //XDrawLine(display, win, gc, x, y+z, x, y);
}

void XGraphics::draw_rectangle(int x, int y, int zx, int zy, long col) {
  Point2d p1(x,y);
  Point2d p2(x+zx,y);
  Point2d p3(x+zx,y+zy);
  Point2d p4(x,y+zy);
  this->draw_line(p1, p2, col);
  this->draw_line(p2, p3, col);
  this->draw_line(p3, p4, col);
  this->draw_line(p4, p1, col);
  //XSetForeground(display, gc, BlackPixel(display, screen_num));
  //XSetLineAttributes(display, gc, 2, LineSolid, 1, 1);
  //XDrawLine(display, win, gc, x, y, x+zx, y);
  //XDrawLine(display, win, gc, x+zx, y, x+zx, y+zy);
  //XDrawLine(display, win, gc, x+zx, y+zy, x, y+zy);
  //XDrawLine(display, win, gc, x, y+zy, x, y);
}

void XGraphics::draw_filled_rectangle(int x, int y, int zx, int zy, long col) {
  XSetForeground(display, gc, col);
  XFillRectangle(display, win, gc, x, height-(y+zy), zx, zy);
}


void XGraphics::draw_faint_line(const Point2d& p1, 
                                const Point2d& p2, 
                                long col){
  XSetForeground(display, gc, (long) 0xDDDDDD);
	XSetLineAttributes(display, gc, 1, LineOnOffDash, 1, 1);
  XDrawLine(display, win, gc, p1.x, p1.y, p2.x, p2.y);
}

void XGraphics::erase_circle(const Point2d& p, int r){
	XSetForeground(display, gc, 0xFFFFFF);
  XSetLineAttributes(display, gc, 1, LineOnOffDash, 1, 1);
	XSetFillStyle(display, gc, FillSolid);
  XFillArc(display, win, gc, p.x-r, p.y-r, 2*r, 2*r, 0, 23040);
}

void XGraphics::draw_circle(const Point2d& p, int r, long col){
    XSetForeground(display, gc, col);
    XSetLineAttributes(display, gc, 1, LineSolid, 1, 1);
    XSetFillStyle(display, gc, FillSolid); 
    XDrawArc(display, win, gc, (p.x-r), (height-p.y)-r, 2*r, 2*r, 0, 23040);
};

void XGraphics::draw_concentric_circles(const Point2d& p, int r, long col){
	int s;
  for(s=3;s>0;s--){
    XSetForeground(display, gc, col*s);
    XSetLineAttributes(display, gc, 1, LineSolid, 1, 1);
    XSetFillStyle(display, gc, FillSolid);
    XDrawArc(display, win, gc, p.x-r*s+1, p.y-r*s+1, 2*r*s-2, 2*r*s-2, 0, 23040);
	}
}

void XGraphics::draw_path(const std::vector<Point2d>& L, long col){
	int i;
	for(i=1; i<(int)L.size(); i++){
		draw_line(L[i-1],L[i],col);
	}
}

void XGraphics::draw_text(const Point2d& p, std::stringstream &T, long col){
  std::string S;
  XSetForeground(display, gc, col);
  S=T.str();
  XDrawString(display,win,gc,p.x,height-p.y,S.c_str(),strlen(S.c_str()));
}

void XGraphics::draw_text(const Point2d& p, std::string &S, long col){
  XSetForeground(display, gc, col);
  XDrawString(display,win,gc,p.x,height-p.y,S.c_str(),strlen(S.c_str()));
}

void XGraphics::draw_label(const Point2d& p, int i, long col){
  XSetForeground(display, gc, col);
  std::stringstream T;
	T << i;
  Point2d draw_loc;
  draw_loc.x = p.x+5;
  draw_loc.y = p.y-5;
	this->draw_text(draw_loc, T, col);
}

std::string XGraphics::wait_for_key() {
  XFlush(display);
  bool finished=false;
  while(finished==false){ 
    XNextEvent(display, &report);
    if (report.type != KeyPress) continue; //ignore the mouse
    if(XLookupKeysym(&report.xkey, 0) == XK_Left){ // left arrow
      return "LA";
    } else if (XLookupKeysym(&report.xkey, 0) == XK_Right){ // right arrow
      return "RA";
    } else if (XLookupKeysym(&report.xkey, 0) == XK_Up){    // up arrow
      return "UA";
    } else if (XLookupKeysym(&report.xkey, 0) == XK_Down){  // down arrow
      return "DA";
    } else if (XLookupKeysym(&report.xkey, 0) == XK_equal){ 
      return "=";
    } else if (XLookupKeysym(&report.xkey, 0) == XK_minus){ 
      return "-";
    } else if (XLookupKeysym(&report.xkey, 0) == XK_o){ 
      return "o";
    } else {
      return "";
    }
  }
  return "u";
}



