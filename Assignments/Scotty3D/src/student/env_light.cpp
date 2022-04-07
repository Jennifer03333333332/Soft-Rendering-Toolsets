
#include "../rays/env_light.h"

#include <limits>

namespace PT {

Vec3 Env_Map::sample() const {

    // TODO (PathTracer): Task 7

    // First, implement Samplers::Sphere::Uniform so the following line works.
    // Second, implement Samplers::Sphere::Image and swap to image_sampler

    //Step one:
    //return uniform_sampler.sample();

    //Steo two:
    return image_sampler.sample();
}

float Env_Map::pdf(Vec3 dir) const {

    // TODO (PathTracer): Task 7


    //Step one:
    //return 1.0f/ (4*PI_F);


    //Steo two:
    float theta = PI_F - std::acos(dir.y);
    theta = std::clamp(theta / PI_F,0.f,1.f); 
    float Jacobian = (float)image.dimension().first * (float)image.dimension().second / (2.f * (float)pow(PI_F,2)*sin(theta));
    return Jacobian * image_sampler.pdf(dir);
}


inline Spectrum lerpSpectrum(float ratio, Spectrum start, Spectrum ends) {
    return (1 - ratio) * start + ratio * ends;
}

Spectrum Env_Map::evaluate(Vec3 dir) const {

    // TODO (PathTracer): Task 7

    // Compute emitted radiance along a given direction by finding the corresponding
    // pixels in the enviornment image. 

    //Vec3 dir -> r, theta, phi -> u,v 

    //Get the coordinate from dir
    float r = dir.norm();// Is dir unit-vector? Yes
    float theta = PI_F - std::acos(dir.y / r);
    float phi = std::atan2(dir.z, dir.x); 
    // You should bi-linearly interpolate the value between the 4 nearest pixels.
    //Theta [0, Π], phi [0,2Π] 
    if (phi < 0) phi = phi + 2.f * PI_F;
    theta = std::clamp(theta / PI_F,0.f,1.f); 
    phi = std::clamp(phi / (2.0f*PI_F),0.f,1.f);  
    
    float h = (float)image.dimension().second;
    float w = (float)image.dimension().first;

    float u = theta * h;//h, u is y
    float v = phi * w;//w, v is x
    //bilinear: (u0,v0),(u0,v1),(u1,v0),(u1,v1)
    float u0 = floor(u), v0 = floor(v);
    float u1,v1;
    if(u - u0<0.5f){
        u1 = u0 - 1;
        std::swap(u0,u1);
    }
    else{
        u1 = u0 + 1.f;
    }
    if(v - v0<0.5f){
        v1 = v0 - 1;
        std::swap(v0,v1);
    }
    else{
        v1 = v0 + 1;
    }
    //Range: >=0 <= w-1 or h-1
    u0 = clamp(u0, 0.f, h - 1.f);
    u1 = clamp(u1, 0.f, h - 1.f);
    v0 = clamp(v0, 0.f, w - 1.f);
    v1 = clamp(v1, 0.f, w - 1.f);
    Spectrum horizontal1 = lerpSpectrum((u-u0)/(u1-u0), image.at((int)v0,(int)u0), image.at((int)v0,(int)u1) );
    Spectrum horizontal2 = lerpSpectrum((u-u0)/(u1-u0), image.at((int)v1,(int)u0), image.at((int)v1,(int)u1) );
    Spectrum vertical =  lerpSpectrum((v-v0)/(v1-v0), horizontal1, horizontal2);
    return vertical;
}

Vec3 Env_Hemisphere::sample() const {
    return sampler.sample();
}

float Env_Hemisphere::pdf(Vec3 dir) const {
    return 1.0f / (2.0f * PI_F);
}

Spectrum Env_Hemisphere::evaluate(Vec3 dir) const {
    if(dir.y > 0.0f) return radiance;
    return {};
}

Vec3 Env_Sphere::sample() const {
    return sampler.sample();
}

float Env_Sphere::pdf(Vec3 dir) const {
    return 1.0f / (4.0f * PI_F);
}

Spectrum Env_Sphere::evaluate(Vec3) const {
    return radiance;
}

} // namespace PT
