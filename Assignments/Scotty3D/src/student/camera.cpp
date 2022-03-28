
#include "../util/camera.h"
#include "../rays/samplers.h"
#include "debug.h"

//screen_coord is in [0,1]^2
Ray Camera::generate_ray(Vec2 screen_coord) const {

    // TODO (PathTracer): Task 1
    // compute the position of the input sensor sample coordinate on the
    // canonical sensor plane one unit away from the pinhole.

    // Tip: Compute the ray direction in camera space and use
    // the camera transform to transform it back into world space.
    

    /// FOV is in degrees
    //float aspect_ratio: width/height
    //float vert_fov: vertical fov
    float screen_height = std::tan(Radians(get_fov()) / 2.0f)*1.0f*2.0f;//z = -1
    float screen_width = aspect_ratio*screen_height;

    //Calculate the coords in camera space (sample_pos_in_cameraspace)
    float point_x = screen_coord.x*screen_width - 0.5f*screen_width;
    float point_y = screen_coord.y*screen_height - 0.5f*screen_height;
    
    Ray r = Ray();
    r.point = Vec3(0,0,0);//start at the origin
    r.dir = Vec3(point_x,point_y,-1.0f);

    //transform it back into world space.
    r.transform(iview);
    return r;
}
