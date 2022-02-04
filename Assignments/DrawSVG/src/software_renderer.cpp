#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  transformation = svg_2_screen;

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;

}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;

}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 5 (part 1):
  // Modify this to implement the transformation stack

  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int) floor(x);
  int sy = (int) floor(y);

  // check bounds
  if ( sx < 0 || sx >= target_w ) return;
  if ( sy < 0 || sy >= target_h ) return;

  // fill sample - NOT doing alpha blending!
  render_target[4 * (sx + sy * target_w)    ] = (uint8_t) (color.r * 255);
  render_target[4 * (sx + sy * target_w) + 1] = (uint8_t) (color.g * 255);
  render_target[4 * (sx + sy * target_w) + 2] = (uint8_t) (color.b * 255);
  render_target[4 * (sx + sy * target_w) + 3] = (uint8_t) (color.a * 255);

}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Task 2: 
  // Implement line rasterization
  // Requires 1 deal with non-int 2 any slope 3 based on length of line
  rasterize_line_Bresenham(x0,y0,x1,y1,color);
}

void SoftwareRendererImp::render_one_point_withoutAlpha( int x, int y, Color& color){
        render_target[4 * (x + y * target_w)    ] = (uint8_t) (color.r * 255);
        render_target[4 * (x + y * target_w) + 1] = (uint8_t) (color.g * 255);
        render_target[4 * (x + y * target_w) + 2] = (uint8_t) (color.b * 255);
        render_target[4 * (x + y * target_w) + 3] = (uint8_t) (color.a * 255);
}

void SoftwareRendererImp::rasterize_line_Bresenham( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {
    //Get the screen coordinate
    int sx0 = (int)floor(x0);      
    int sx1 = (int)floor(x1);    
    int sy0 = (int)floor(y0);      
    int sy1 = (int)floor(y1);  
    //Deal with different slope
    int dx = abs(sx1 - sx0), sx = sx0 < sx1 ? 1 : -1;
    int dy = abs(sy1 - sy0), sy = sy0 < sy1 ? 1 : -1;
    int error = (dx > dy ? dx : -dy) / 2;
    int error_tmp;


    while(true) {
        render_one_point_withoutAlpha(sx0,sy0,color);
        if (sx0 == sx1 && sy0 == sy1) break;//draw line ends
        error_tmp = error;
        if (error_tmp > -dx) { error -= dy; sx0 += sx; }
        if (error_tmp < dy) { error += dx; sy0 += sy; }
    }
}


void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  // Task 3: 
  // Implement triangle rasterization
    //Simple solution: a big 2D-bounding box
    int xmin = (int)std::min(x0, std::min(x1, x2));
    int ymin = (int)std::min(y0, std::min(y1, y2));
    int xmax = (int)std::max(x0, std::max(x1, x2)) + 1;
    int ymax = (int)std::max(y0, std::max(y1, y2)) + 1;

    vector<Vector2D> Triangle; 
    
    Triangle[0] = Vector2D(x0, y0);
    Triangle[1] = Vector2D(x1, y1);
    Triangle[2] = Vector2D(x2, y2);
    for (float x = xmin; x <= xmax; x++) {
        for (float y = ymin; y <= ymax; y++) {
            if (inside_triangle(Vector2D(x, y), Triangle)) {
                render_one_point_withoutAlpha(x, y, color);
            }
        }
    }
}


bool SoftwareRendererImp::inside_triangle(const Vector2D& pointP, const vector<Vector2D>& Triangle) {
    //Calculate the 3 edges
    Vector2D t0t1 = Triangle[1] - Triangle[0];
    Vector2D t1t2 = Triangle[2] - Triangle[1];
    Vector2D t2t0 = Triangle[0] - Triangle[2];
    //Calculate vectors form the vertex to the point p 
    Vector2D v0p = pointP - Triangle[0];
    Vector2D v1p = pointP - Triangle[1];
    Vector2D v2p = pointP - Triangle[2];
    //Calculate the cross product
    float cross1 = t0t1[0] * t0p[1] - t0t1[1] * t0p[0];
    float cross2 = t1t2[0] * t1p[1] - t1t2[1] * t1p[0];
    float cross3 = t2t0[0] * t2p[1] - t2t0[1] * t2p[0];
    //! a sample point on a triangle edge is covered by the triangle. In this case, cross product would be 0
    //! the vertices may be in counter-clockwise or clockwise order
    bool counter_cw = (cross1 * cross2 >= 0) && (cross2 * cross3 >= 0) && (cross1 * cross3 >= 0);
    bool cw = (cross1 * cross2 <= 0) && (cross2 * cross3 <= 0) && (cross1 * cross3 <= 0);
    return counter_cw || cw;
}


void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 6: 
  // Implement image rasterization

}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 4: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 4".
  return;

}


} // namespace CMU462
