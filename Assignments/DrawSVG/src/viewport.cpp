#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float centerX, float centerY, float vspan ) {

  // Task 5 (part 2): 
  // Set svg coordinate to normalized device coordinate transformation. Your input
  // arguments are defined as normalized SVG canvas coordinates.
  this->centerX = centerX;
  this->centerY = centerY;
  this->vspan = vspan; 

  Matrix3x3 Torigin = Matrix3x3::identity();
  Torigin[2][0] = 0 - centerX;
  Torigin[2][1] = 0 - centerY;
  Matrix3x3 Scale = Matrix3x3::identity();
  Scale[0][0] = 0.5 / vspan;//equals to y scaling
  Scale[1][1] = 0.5 / vspan;//from 2*vspan -> 1
  Matrix3x3 Thalfhalf = Matrix3x3::identity();
  Thalfhalf[2][0] = 0.5;
  Thalfhalf[2][1] = 0.5;
  svg_2_norm = Thalfhalf * Scale * Torigin;
}
//When update_viewbox, it would call DrawSVG::redraw().
void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->centerX -= dx;//float dx = (x - cursor_x) / width  * tabs[current_tab]->width;
  this->centerY -= dy;//    float dy = (y - cursor_y) / height * tabs[current_tab]->height;
  this->vspan *= scale;
  set_viewbox( centerX, centerY, vspan );
}

} // namespace CMU462
