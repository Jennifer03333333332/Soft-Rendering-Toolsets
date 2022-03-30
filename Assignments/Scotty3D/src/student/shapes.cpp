
#include "../rays/shapes.h"
#include "debug.h"

namespace PT {

const char* Shape_Type_Names[(int)Shape_Type::count] = {"None", "Sphere"};

BBox Sphere::bbox() const {

    BBox box;
    box.enclose(Vec3(-radius));
    box.enclose(Vec3(radius));
    return box;
}

Trace Sphere::hit(const Ray& ray) const {

    // TODO (PathTracer): Task 2
    // Intersect this ray with a sphere of radius Sphere::radius centered at the origin.
    bool result = true;
    float t = 0;
    //Quadratic equation
    float a = (ray.dir).norm_squared();
    float b = 2.0f*dot(ray.point,ray.dir);
    float c = (ray.point).norm_squared() - radius*radius;

    float delta = b*b - 4.0f*a*c;//b^2 - 4ac
    //2 intersection points
    if(delta > 0){
        bool t1_valid = true, t2_valid = true;
        float t1 = ((-2.0f)*dot(ray.point,ray.dir) + sqrt(delta))/(2.0f*(ray.dir).norm_squared());
        float t2 = ((-2.0f)*dot(ray.point,ray.dir) - sqrt(delta))/(2.0f*(ray.dir).norm_squared());
        //check t > 0
        t1_valid = (t1<0)? false:true;
        t2_valid = (t2<0)? false:true;
        //check bounds
        //4 Cases: t1 out of bounds but t2 in the bounds; t1 in the bounds but t2 outof bounds
        float distance1 = std::abs((t1 * ray.dir).norm());
        float distance2 = std::abs((t2 * ray.dir).norm());
        if(distance1 < ray.dist_bounds.x || distance1 > ray.dist_bounds.y) {
            t1_valid = false;
        }
        if(distance2 < ray.dist_bounds.x || distance2 > ray.dist_bounds.y) {
            t2_valid = false;
        }
        // Conclusion
        if( t1_valid && t2_valid ){
            t = std::min(t1,t2);
        }
        else if((!t1_valid) && (!t2_valid)){
            result = false;
        }
        else{
            t = (t1_valid)? t1:t2;
        }

    }
    else if (delta == 0){//1 intersection point: tangent
        t = ((-2.0f)*dot(ray.point,ray.dir))/(2.0f*(ray.dir).norm_squared());
    }
    else{// No intersection
        result = false;
    }
    // If the ray intersects the sphere twice, ret should
    // represent the first intersection, but remember to respect
    // ray.dist_bounds! For example, if there are two intersections,
    // but only the _later_ one is within ray.dist_bounds, you should
    // return that one!

    Trace ret;
    ret.origin = ray.point;
    ret.hit = result;       // was there an intersection?
    if(result){
        ret.distance = std::abs((ray.at(t) - ray.point).norm());   // at what distance did the intersection occur?
        ret.position = ray.at(t); // where was the intersection?
        ret.normal = ray.at(t) - Vec3();   // normals should be out-facing.
    }
    return ret;
}

} // namespace PT
