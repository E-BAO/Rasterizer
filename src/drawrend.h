#ifndef CGL_DRAWREND_H
#define CGL_DRAWREND_H

#include "CGL/CGL.h"
#include "CGL/renderer.h"
#include "CGL/color.h"
#include <vector>
#include <cstring>
#include "GLFW/glfw3.h"
#include "svg.h"

namespace CGL {

class DrawRend : public Renderer {
 public:
  DrawRend(std::vector<SVG*> svgs_):
  svgs(svgs_), current_svg(0)
  {}

  ~DrawRend( void );

  // inherited Renderer interface functions
  void init();
  void render();
  void resize( size_t w, size_t h );
  std::string name() { return "Draw"; }
  std::string info();
  void cursor_event( float x, float y );
  void scroll_event( float offset_x, float offset_y );
  void mouse_event( int key, int event, unsigned char mods );
  void keyboard_event( int key, int event, unsigned char mods );

  void set_gl(bool gl_) { gl = gl_; }

  // write current pixel buffer to disk
  void write_screenshot();

  // write only framebuffer to disk
  void write_framebuffer();

  // drawing functions
  void redraw();
  void draw_pixels();
  void draw_zoom();

  // view transform functions
  void view_init();
  void set_view(float x, float y, float span);
  void move_view(float dx, float dy, float scale);

  // rasterize a point
  void rasterize_point( float x, float y, Color color );

  // rasterize a line
  void rasterize_line( float x0, float y0,
                       float x1, float y1,
                       Color color);

  // rasterize a triangle
  void rasterize_triangle( float x0, float y0,
                           float x1, float y1,
                           float x2, float y2,
                           Color color, Triangle *tri = NULL );

    bool isInside(Vector2D& p, vector<Vector2D>& vertexes, vector<Vector2D>& normals);
    bool isInside(double x, double y, vector<Vector2D>& vertexes, vector<Vector2D>& normals);
    Vector3D get_p_bary(float x0, float y0,
                       float x1, float y1,
                       float x2, float y2,
                       float xp, float yp);
//    void rasterize_triangle_half(float x0, float y0, float dy, float Ldx, float Rdx, Color color);

private:
  // Global state variables for SVGs, pixels, and view transforms
  std::vector<SVG*> svgs; size_t current_svg;
  std::vector<Matrix3x3> svg_to_ndc;
  float view_x, view_y, view_span;

  Matrix3x3 ndc_to_screen;

  std::vector<unsigned char> framebuffer;
  size_t width, height;

  // UI state info
  float cursor_x; float cursor_y;
  bool left_clicked;
  int show_zoom;
  int sample_rate;

  PixelSampleMethod psm;
  LevelSampleMethod lsm;

  bool gl;

  typedef std::vector<unsigned char> PixelColorStorage;

  // Intuitively, a sample buffer instance is a pixel,
  // or (samples_per_side x samples_per_side) sub-pixels.
  struct SampleBuffer {
    std::vector<std::vector<PixelColorStorage> > sub_pixels;
    size_t samples_per_side;

    SampleBuffer(size_t sps): samples_per_side(sps) {
      clear();
    }
    
    // Fill the subpixel at i,j with the Color c
    void fill_color(int i, int j, Color c) {
      PixelColorStorage &p = sub_pixels[i][j];
      // Part 1: Overwrite PixelColorStorage p using Color c.
      //         Pay attention to different data types.
        p[0] = max(0.,min(c.r * 255.0, 255.0));
        p[1] = max(0.,min(c.g * 255.0, 255.0));
        p[2] = max(0.,min(c.b * 255.0, 255.0));
      return;
    }

    void fill_pixel(Color c) {
      for (int i = 0; i < samples_per_side; ++i)
        for (int j = 0; j < samples_per_side; ++j)
          fill_color(i, j, c);
    }

    Color get_pixel_color() {
//        return Color(sub_pixels[0][0].data());
      // Part 2: Implement get_pixel_color() for supersampling.
        double r = 0, g = 0, b = 0;
        for(int i = 0; i < samples_per_side; i ++){
            for(int j = 0; j < samples_per_side; j ++){
                r += sub_pixels[i][j][0];
                g += sub_pixels[i][j][1];
                b += sub_pixels[i][j][2];
            }
        }
        double div = (double)samples_per_side * (double)samples_per_side * 255.0;
        r /= div;
        g /= div;
        b /= div;
        return Color(r,g,b);
    }
    
    void clear() {
      if (sub_pixels.size() == samples_per_side) {
        for (int i = 0; i < samples_per_side; ++i)
          for (int j = 0; j < samples_per_side; ++j)
            sub_pixels[i][j].assign(3, (unsigned char)255);
        return;
      }

      sub_pixels.clear();
      PixelColorStorage white = std::vector<unsigned char>(3, 255);
      std::vector<PixelColorStorage> row;
      row.reserve(samples_per_side);
      for (int i = 0; i < samples_per_side; ++i)
        row.push_back(white);
      sub_pixels.reserve(samples_per_side);
      for (int i = 0; i < samples_per_side; ++i)
        sub_pixels.push_back(row);
    }
  };

  std::vector<std::vector<SampleBuffer> > samplebuffer;

  // This function takes the collected sub-pixel samples and
  // combines them together to fill in the framebuffer in preparation
  // for posting pixels to the screen.
  void resolve() {
    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col = samplebuffer[y][x].get_pixel_color();
        for (int k = 0; k < 3; ++k) {
          framebuffer[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }
  }


};

} // namespace CGL

#endif // CGL_DRAWREND_H
