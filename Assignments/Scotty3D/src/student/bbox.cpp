
#include "../lib/mathlib.h"
#include "debug.h"

bool BBox::hit(const Ray& ray, Vec2& times) const {

    // TODO (PathTracer): Task 3
    // Implement ray - bounding box intersection test
    // Scratchapixel
    bool result = false;

    float tmin, tmax, tymin, tymax, tzmin, tzmax; 
    Vec3 invdir = 1.0f / ray.dir;
    tmin = (invdir.x<0)? ((max.x - ray.point.x) * invdir.x) : ((min.x - ray.point.x) * invdir.x);
    tmax = (1-(invdir.x<0))? ((max.x - ray.point.x) * invdir.x) : ((min.x - ray.point.x) * invdir.x);
    tymin = (invdir.y<0) ? ((max.y - ray.point.y) * invdir.y) : ((min.y - ray.point.y) * invdir.y);
    tymax = (1-(invdir.y<0)) ? ((max.y - ray.point.y) * invdir.y) : ((min.y - ray.point.y) * invdir.y);
 
    if ((tmin > tymax) || (tymin > tmax)) 
        return false; 
    if (tymin > tmin) 
        tmin = tymin; 
    if (tymax < tmax) 
        tmax = tymax; 
 
    tzmin = (invdir.z < 0) ? ((max.z - ray.point.z) * invdir.z):((min.z - ray.point.z) * invdir.z);

    tzmax = (1-(invdir.z < 0)) ? ((max.z - ray.point.z) * invdir.z): ((min.z - ray.point.z) * invdir.z);
 
    if ((tmin > tzmax) || (tzmin > tmax)) 
        return false; 
    if (tzmin > tmin) 
        tmin = tzmin; 
    if (tzmax < tmax) 
        tmax = tzmax; 
 
    // If the ray intersected the bounding box within the range given by
    // [times.x,times.y], update times with the new intersection times.
    //TODO
    if( tmin >= times.x && tmin <= times.y){
        times.x = tmin;
    }
    if(tmax >= times.x && tmax <= times.y){
        times.y = tmax;
    }

    return true;
}
