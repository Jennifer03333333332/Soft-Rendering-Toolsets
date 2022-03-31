
#include "../lib/mathlib.h"
#include "debug.h"

bool BBox::hit(const Ray& ray, Vec2& times) const {

    // TODO (PathTracer): Task 3
    // Implement ray - bounding box intersection test
    // Scratchapixel

    float tmin, tmax, tymin, tymax, tzmin, tzmax; 
    Vec3 invdir = 1.0f / ray.dir;
    int sign[3];
    sign[0] = (invdir.x < 0); 
    sign[1] = (invdir.y < 0); 
    sign[2] = (invdir.z < 0); 
    Vec3 bounds[2];
    bounds[0] = min;
    bounds[1] = max;


    tmin = (bounds[sign[0]].x - ray.point.x) * invdir.x; 
    tmax = (bounds[1-sign[0]].x - ray.point.x) * invdir.x; 
    tymin = (bounds[sign[1]].y - ray.point.y) * invdir.y; 
    tymax = (bounds[1-sign[1]].y - ray.point.y) * invdir.y; 

    // tmin = (invdir.x<0)? ((max.x - ray.point.x) * invdir.x) : ((min.x - ray.point.x) * invdir.x);
    // tmax = (1-(invdir.x<0))? ((max.x - ray.point.x) * invdir.x) : ((min.x - ray.point.x) * invdir.x);
    // tymin = (invdir.y<0) ? ((max.y - ray.point.y) * invdir.y) : ((min.y - ray.point.y) * invdir.y);
    // tymax = (1-(invdir.y<0)) ? ((max.y - ray.point.y) * invdir.y) : ((min.y - ray.point.y) * invdir.y);
 
    if ((tmin > tymax) || (tymin > tmax)) 
        return false; 
    if (tymin > tmin) 
        tmin = tymin; 
    if (tymax < tmax) 
        tmax = tymax; 
    tzmin = (bounds[sign[2]].z - ray.point.z) * invdir.z; 
    tzmax = (bounds[1-sign[2]].z - ray.point.z) * invdir.z; 
    // tzmin = (invdir.z < 0) ? ((max.z - ray.point.z) * invdir.z):((min.z - ray.point.z) * invdir.z);

    // tzmax = (1-(invdir.z < 0)) ? ((max.z - ray.point.z) * invdir.z): ((min.z - ray.point.z) * invdir.z);
 
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
