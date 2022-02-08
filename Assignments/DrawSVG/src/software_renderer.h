#ifndef CMU462_SOFTWARE_RENDERER_H
#define CMU462_SOFTWARE_RENDERER_H

#include <stdio.h>
#include <vector>

#include "CMU462.h"
#include "texture.h"
#include "svg_renderer.h"

namespace CMU462 { // CMU462

    //struct Color4i {
    //    float r, g, b, a;
    //    Color4i(){r=0,g=0,b=0,a=0;}
    //    Color4i(float mr, float mg, float mb, float ma) :r(mr), g(mg), b(mb), a(ma) {}
    //    Color4i operator+(Color4i& other){
    //      return Color4i(r+other.r,g+other.g,b+other.b,a+other.a);
    //    }
    //    Color4i operator/(int nums){
    //      return Color4i(r/nums,g/nums,b/nums,a/nums);
    //    }
    //};

class SoftwareRenderer : public SVGRenderer {
 public:

  SoftwareRenderer( ) : sample_rate (1) { }

  // Free used resources
  virtual ~SoftwareRenderer( ) { }

  // Draw an svg input to render target
  virtual void draw_svg( SVG& svg ) = 0;

  // Set sample rate
  virtual void set_sample_rate( size_t sample_rate ) = 0;
  
  // Set render target
  virtual void set_render_target( unsigned char* render_target,
                                  size_t width, size_t height ) = 0;

  // Clear render target
  inline void clear_target() {
    memset(render_target, 255, 4 * target_w * target_h);
  }

  // Set texture sampler
  inline void set_tex_sampler( Sampler2D* sampler ) {
    this->sampler = sampler;
  }

  // Set svg to screen transformation
  inline void set_svg_2_screen( Matrix3x3 svg_2_screen ) {
    this->svg_2_screen = svg_2_screen;
  }

 protected:

  // Sample rate (square root of samples per pixel)
  size_t sample_rate;

  // Render target memory location
  unsigned char* render_target; 

  // Target buffer dimension (in pixels)
  size_t target_w; size_t target_h;

  // Texture sampler being used
  Sampler2D* sampler;

  // SVG coordinates to screen space coordinates
  Matrix3x3 svg_2_screen;

}; // class SoftwareRenderer


class SoftwareRendererImp : public SoftwareRenderer {
 public:

  SoftwareRendererImp( ) : SoftwareRenderer( ) { }

  // draw an svg input to render target
  void draw_svg( SVG& svg );

  // set sample rate
  void set_sample_rate( size_t sample_rate );
  
  // set render target
  void set_render_target( unsigned char* target_buffer,
                          size_t width, size_t height );
  
  inline void clear_target() {
    memset(render_target, 255, 4 * target_w * target_h);
    super_sample_buffer.resize(4*supersample_h * supersample_w);
    memset(&super_sample_buffer[0], (float)255, 4 * supersample_h * supersample_w * sizeof(float));
  }
 private:

  // Primitive Drawing //

  // Draws an SVG element
  void draw_element( SVGElement* element );

  // Draws a point
  void draw_point( Point& p );

  // Draw a line
  void draw_line( Line& line );

  // Draw a polyline
  void draw_polyline( Polyline& polyline );

  // Draw a rectangle
  void draw_rect ( Rect& rect );

  // Draw a polygon
  void draw_polygon( Polygon& polygon );

  // Draw a ellipse
  void draw_ellipse( Ellipse& ellipse );

  // Draws a bitmap image
  void draw_image( Image& image );

  // Draw a group
  void draw_group( Group& group );

  // Rasterization //

  // rasterize a point
  void rasterize_point(double x, double y, const Color& color);

  // rasterize a line
  void rasterize_line( float x0, float y0,
                       float x1, float y1,
                       Color color);
  //Using Bresenham's algorithm to rasterize a line
  void rasterize_line_Bresenham( float x0, float y0,
                       float x1, float y1,
                       Color color);
  void rasterize_line_xiaolinwu(float x0, float y0,
          float x1, float y1,
          Color color);

  //Judge if a point is inside the triangle
  bool inside_triangle(const Vector2D& pointP, const std::vector<Vector2D>& Triangle);

  // rasterize a triangle
  void rasterize_triangle( float x0, float y0,
                           float x1, float y1,
                           float x2, float y2,
                           Color color );

  // rasterize an image
  void rasterize_image( float x0, float y0,
                        float x1, float y1,
                        Texture& tex );

  // resolve samples to render target
  void resolve( void );
  //For Task 4
  //unsigned char* super_sample_buffer;
  //std::vector<unsigned char> super_sample_buffer;
  //std::vector<uint8_t> super_sample_buffer;
  //try use float
  bool SSAA = true;
  std::vector<float> super_sample_buffer;//intermidate buffer
  size_t supersample_w; size_t supersample_h;
  inline void fill_sample(int sx, int sy, const Color& c);
  inline void fill_pixel(int x, int y, const Color& c);
}; // class SoftwareRendererImp


class SoftwareRendererRef : public SoftwareRenderer {
 public:

  SoftwareRendererRef( ) : SoftwareRenderer( ) { }

  // draw an svg input to render target
  void draw_svg( SVG& svg );

  // set sample rate
  void set_sample_rate( size_t sample_rate );
  
  // set render target
  void set_render_target( unsigned char* target_buffer,
                          size_t width, size_t height );

 private:

  // Primitive Drawing //

  // Draws an SVG element
  void draw_element( SVGElement* element );

  // Draws a point
  void draw_point( Point& p );

  // Draw a line
  void draw_line( Line& line );

  // Draw a polyline
  void draw_polyline( Polyline& polyline );

  // Draw a rectangle
  void draw_rect ( Rect& rect );

  // Draw a polygon
  void draw_polygon( Polygon& polygon );

  // Draw a ellipse
  void draw_ellipse( Ellipse& ellipse );

  // Draws a bitmap image
  void draw_image( Image& image );

  // Draw a group
  void draw_group( Group& group );

  // Rasterization //

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
                           Color color );

  // rasterize an image
  void rasterize_image( float x0, float y0,
                        float x1, float y1,
                        Texture& tex );

  // resolve samples to render target
  void resolve( void );

  // Helpers //
  // HINT: you may want to have something similar //
  std::vector<unsigned char> sample_buffer; int w; int h;
  void fill_sample( int sx, int sy, const Color& c );
  void fill_pixel( int x, int y, const Color&c );

}; // class SoftwareRendererRef

} // namespace CMU462

#endif // CMU462_SOFTWARE_RENDERER_H
