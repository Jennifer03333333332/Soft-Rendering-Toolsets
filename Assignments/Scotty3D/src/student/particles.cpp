
#include "../scene/particles.h"
#include "../rays/pathtracer.h"

bool Scene_Particles::Particle::update(const PT::Object& scene, float dt, float radius) {

    // TODO(Animation): Task 4

    // Compute the trajectory of this particle for the next dt seconds.

    // (1) Build a ray representing the particle's path if it travelled at constant velocity.

    // (2) Intersect the ray with the scene and account for collisions. Be careful when placing
    // collision points using the particle radius. Move the particle to its next position.
    
    float remain_timestep = dt;
    
    while(remain_timestep > 0){
        Ray r = Ray();
        r.point = pos;//start at the particle pos
        r.dir = velocity;//dir = velocity, should it be normalized?? no. Should it be abs()? no

        PT::Trace trace = scene.hit(r);
        //consider radius
        float cos_t = dot(trace.normal,-1.0f*velocity)/(velocity.norm() * trace.normal.norm());
        Vec3 surface_normal = trace.normal / trace.normal.norm();
        //cos_t == 0: no intersection? cos_t <0: normal is on the other side
        if(cos_t <0){
            cos_t = sqrt(1-cos_t*cos_t);
            surface_normal = -1.0f*surface_normal;
        }

        float interval = std::abs(radius/cos_t);
        float hit_time = (trace.distance - interval)/r.dir.norm();//trace.distance/r.dir.norm(); //(trace.position - trace.origin)/r.dir;
            
        //if not hit || hit_time > remain dt, return
        if(!trace.hit || hit_time > remain_timestep || cos_t == 0){//Not hit in this time step
            //move the particle's a little
            pos = pos + velocity*remain_timestep;
            //velocity
            velocity = velocity + acceleration*remain_timestep;
            break;
        }    
        //If hit
        //1 put particle at the hit point
        pos = trace.position - interval*velocity/velocity.norm();
        //2 reflect the velocity
        //should the surface normal be normalized?
        Vec3 reflect = velocity - 2.0f* dot(velocity,surface_normal)*surface_normal;
        velocity = reflect;
        // (3) Account for acceleration due to gravity.
        //velocity = velocity + acceleration*hit_time;
        remain_timestep -= hit_time;
    }
    // (4) Repeat until the entire time step has been consumed.
    // (5) Decrease the particle's age and return whether it should die.
    age-=dt;
    return age > 0;
}
