#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float centerX, float centerY, float vspan ) {

  // Task 5 (part 2): 
  // Set svg coordinate to normalized device coordinate transformation. Your input
  // arguments are defined as normalized SVG canvas coordinates.

	//This function should construct a Matrix3x3 object (Homogenous coordinate for 2D)
	//that transforms the viewport area to normalized screen space.

	//when update_viewbox, it would call DrawSVG::redraw().
	//Matrix3x3 m_imp = norm_to_screen * viewport_imp[current_tab]->get_svg_2_norm();
	//norm_to_screen£º normalized coordiantes to screen coordinates
	//get_svg_2_norm(): svg_2_norm: SVG coordinate to normalized display coordinates
  this->centerX = centerX;
  this->centerY = centerY;
  this->vspan = vspan; 
  //SVG coordinate to normalized display coordinates
  //svg_2_norm is row-major
  //first move to origin, then scale, then translate to (0.5,0.5)

  //my code
  //Matrix3x3 Torigin = Matrix3x3::identity();
  //Torigin[2][0] = 0 - centerX;
  //Torigin[2][1] = 0 - centerY;
  //Matrix3x3 Scale = Matrix3x3::identity();
  //Scale[0][0] = 0.5 / centerX;
  //Scale[1][1] = 0.5 / centerY;
  //Matrix3x3 Thalfhalf = Matrix3x3::identity();
  //Thalfhalf[2][0] = 0.5;
  //Thalfhalf[2][1] = 0.5;
  //svg_2_norm = Thalfhalf * Scale * Torigin;

  svg_2_norm = Matrix3x3::identity();
  svg_2_norm[0][0] = 0.5 / vspan;
  svg_2_norm[1][1] = 0.5 / vspan;
  svg_2_norm[2][0] = 0.5 - 0.5 * centerX / vspan;
  svg_2_norm[2][1] = 0.5 - 0.5 * centerY / vspan;

}
//    
void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->centerX -= dx;//float dx = (x - cursor_x) / width  * tabs[current_tab]->width;
  this->centerY -= dy;//    float dy = (y - cursor_y) / height * tabs[current_tab]->height;
  this->vspan *= scale;
  set_viewbox( centerX, centerY, vspan );
}

} // namespace CMU462
