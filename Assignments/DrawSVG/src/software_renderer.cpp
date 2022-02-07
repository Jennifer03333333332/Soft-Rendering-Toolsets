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
  //clear
  clear_target();
  // set top level transformation
  transformation = svg_2_screen;
  //this->super_sample_buffer = std::vector<uint8_t>(4 * this->supersample_h * this->supersample_w, 255);

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

//Called when changed sample_rate(- and =)
void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 4: 
  // You may want to modify this for supersampling support

    if (this->sample_rate == sample_rate)return;
    this->sample_rate = sample_rate;
    supersample_h = target_h * sample_rate;
    supersample_w = target_w * sample_rate;
    //allocate memory for super_sample_buffer, sizeof rgba * number of samples
    //memset(super_sample_buffer, 0, 4 * sizeof(uint8_t) * this->supersample_h * this->supersample_w);
    //this->super_sample_buffer = std::vector<uint8_t>(4 * this->supersample_h * this->supersample_w, 255);
    
    
    super_sample_buffer.resize(4*supersample_h * supersample_w);

    memset(&super_sample_buffer[0],(float)255,4 * supersample_h * supersample_w*sizeof(float));

}

//Called when resizes
void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  target_w = width;
  target_h = height;

  supersample_h = height * sample_rate;
  supersample_w = width * sample_rate;
  //initialize the buffer with white(255)
  // this->super_sample_buffer.clear();
  super_sample_buffer.resize(4*supersample_h * supersample_w);//if 'this' is needed??
  // //construct could take much room
  //this->super_sample_buffer = std::vector<uint8_t>(4 * this->supersample_h * this->supersample_w, 255);
  //memset(&super_sample_buffer[0],255,4 * supersample_h * supersample_w);
  
  memset(&super_sample_buffer[0], (float)255, 4 * supersample_h * supersample_w * sizeof(float));


  //super_sample_buffer.resize(4 * this->supersample_h * this->supersample_w);
  //memset(super_sample_buffer, 255, 4 * this->supersample_h * this->supersample_w*sizeof(Color4i));
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
// after I changed rasterize_triangle, x,y turn to sample position in a screen
void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int) floor(x);
  int sy = (int) floor(y);

  // check bounds
  if ( sx < 0 || sx >= target_w ) return;
  if ( sy < 0 || sy >= target_h ) return;

  //For super-sample

  // fill sample - NOT doing alpha blending!
  //render_target[4 * (sx + sy * target_w)    ] = (uint8_t) (color.r * 255);
  //render_target[4 * (sx + sy * target_w) + 1] = (uint8_t) (color.g * 255);
  //render_target[4 * (sx + sy * target_w) + 2] = (uint8_t) (color.b * 255);
  //render_target[4 * (sx + sy * target_w) + 3] = (uint8_t) (color.a * 255);

  fill_sample(x*sample_rate, y * sample_rate, color);
}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Task 2: 
  // Implement line rasterization
  // Requires 1 deal with non-int 2 any slope 3 based on length of line
  rasterize_line_Bresenham(x0,y0,x1,y1,color);
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
        rasterize_point(sx0,sy0,color);
        if (sx0 == sx1 && sy0 == sy1) break;//draw line ends
        error_tmp = error;
        if (error_tmp > -dx) { error -= dy; sx0 += sx; }
        if (error_tmp < dy) { error += dx; sy0 += sy; }
    }
}

//already screen coordinates
void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  // Task 3: 
  // Implement triangle rasterization
    //Build an array to store the triangle's vertices
    std::vector<Vector2D> Triangle; 
    Triangle.push_back(Vector2D((double)x0, (double)y0));
    Triangle.push_back(Vector2D((double)x1, (double)y1));
    Triangle.push_back(Vector2D((double)x2, (double)y2));
    //Simple solution: a big 2D-bounding box
    float xmin = floor(std::min(x0, std::min(x1, x2)));
    float ymin = floor(std::min(y0, std::min(y1, y2)));
    float xmax = ceil(std::max(x0, std::max(x1, x2)));
    float ymax = ceil(std::max(y0, std::max(y1, y2)));
    
    //Seperate the normal implementation and the SSAA implementation just to show the difference
    if (sample_rate == 1) {
        //For each pixel in this bounding box
        for (double x = xmin; x <= xmax; x++) {
            for (double y = ymin; y <= ymax; y++) {
                //center of the pixels: x+0.5,y+0.5
                if (inside_triangle(Vector2D(x + 0.5, y + 0.5), Triangle)) {
                    rasterize_point(x, y, color);
                }
            }
        }
    }
    //SSAA
    else {
        for (double x = xmin; x <= xmax; x += 1.0/ sample_rate ) {
            for (double y = ymin; y <= ymax; y += 1.0 / sample_rate) {
                //center of the samples x + 0.5 / sample_rate,y+0.5 / sample_rate
                if (inside_triangle(Vector2D(x + 0.5 / sample_rate, y + 0.5 / sample_rate), Triangle)) {
                    rasterize_point(x, y, color);
                }
            }
        }
    }

}

//Judge if the point is in the triangle
bool SoftwareRendererImp::inside_triangle(const Vector2D& pointP, const std::vector<Vector2D>& Triangle) {
    //Calculate the 3 edges
    Vector2D t0t1 = Triangle[1] - Triangle[0];
    Vector2D t1t2 = Triangle[2] - Triangle[1];
    Vector2D t2t0 = Triangle[0] - Triangle[2];
    //Calculate vectors form the vertex to the point p 
    Vector2D t0p = pointP - Triangle[0];
    Vector2D t1p = pointP - Triangle[1];
    Vector2D t2p = pointP - Triangle[2];
    //Calculate the cross product
    float cross1 = t0t1.x * t0p.y - t0t1.y * t0p.x;
    float cross2 = t1t2.x * t1p.y - t1t2.y * t1p.x;
    float cross3 = t2t0.x * t2p.y - t2t0.y * t2p.x;
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


    //1 Have render_target already. Create a new super_sample_buffer. the length would be 2*2 times than render_target
    //2 Rasterize based on super_sample_buffer
    //Resolve it to render_target

    //For each pixel, project the super_sample_buffer to it. samplerate*samplerate -> 1


    //For each sample
    for (size_t x = 0; x < supersample_w; x += sample_rate) {
        for (size_t y = 0; y < supersample_h; y += sample_rate) {
            //uint8_t r = 0, g = 0, b = 0, a = 0;
            float r = 0, g = 0, b = 0, a = 0;
            //Calculate all samples in the pixel where this sample in
            for (size_t i = 0; i < sample_rate; ++i){
                for (size_t j = 0; j < sample_rate; ++j){
                    size_t samplePos = 4 * (x + i + (y + j) * supersample_w);
                    r += super_sample_buffer[samplePos];//(sample_rate * sample_rate);lost precision
                    g += super_sample_buffer[samplePos + 1];
                    b += super_sample_buffer[samplePos + 2];
                    a += super_sample_buffer[samplePos + 3];
                }
            }
            r/=sample_rate * sample_rate;g/=sample_rate * sample_rate;b/=sample_rate * sample_rate;a/=sample_rate * sample_rate;
            //SSAABuffer position to pixel position. (x,y) is the left-buttom sample of this pixel
            size_t pixelPos = 4 * ((x / sample_rate) + (y / sample_rate) * target_w);
            render_target[pixelPos] = (uint8_t)(r);
            render_target[pixelPos + 1] = (uint8_t)(g);
            render_target[pixelPos + 2] = (uint8_t)(b);
            render_target[pixelPos + 3] = (uint8_t)(a);
        }
    }
    //this->super_sample_buffer = std::vector<uint8_t>(4 * this->supersample_h * this->supersample_w, 255);
    //memset(&super_sample_buffer[0],255,4 * this->supersample_h * this->supersample_w);
    return;

}

//fill the sample buffer
//sx,sy are the screen coordinates for ssaa screen
inline void SoftwareRendererImp::fill_sample(int sx, int sy, const Color& c) {
    //super_sample_buffer[4*(sx + sy * supersample_w)] = (uint8_t)c.r*255;
    //super_sample_buffer[4*(sx + sy * supersample_w)+1] = (uint8_t)c.g*255;
    //super_sample_buffer[4*(sx + sy * supersample_w)+2] = (uint8_t)c.b*255;
    //super_sample_buffer[4*(sx + sy * supersample_w)+3] = (uint8_t)c.a*255;

    //Method n: try store float
    super_sample_buffer[4*(sx + sy * supersample_w)] = c.r*255.0;
    super_sample_buffer[4*(sx + sy * supersample_w)+1] = c.g*255.0;
    super_sample_buffer[4*(sx + sy * supersample_w)+2] = c.b*255.0;
    super_sample_buffer[4*(sx + sy * supersample_w)+3] = c.a*255.0;
} // namespace CMU462

//inline void SoftwareRendererImp::fill_pixel(int x, int y, const CMU462::Color4i& c){
//    render_target[4 * (x + y * target_w)] = (uint8_t)c.r*255;
//    render_target[4 * (x + y * target_w) + 1] = (uint8_t)c.g*255;
//    render_target[4 * (x + y * target_w) + 2] = (uint8_t)c.b*255;
//    render_target[4 * (x + y * target_w) + 3] = (uint8_t)c.a*255;
//}

}