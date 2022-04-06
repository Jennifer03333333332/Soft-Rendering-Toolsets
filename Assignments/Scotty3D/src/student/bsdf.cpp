
#include "../rays/bsdf.h"
#include "../util/rand.h"

namespace PT {
//given wi
static Vec3 reflect(Vec3 dir) {

    // TODO (PathTracer): Task 5
    // Return reflection of dir about the surface normal (0,1,0)
    //Vec3 normal(0,1,0); 
    Vec3 in_dir = Vec3((-1.0f)*dir.x,dir.y,(-1.0f)*dir.z);//?? should z be the -z?
    //Vec3 w0 = (-1.0f*dir) + 2.0f * dot(dir,normal)* normal;
    return in_dir;
}

double Schlick_Approximation(float cosine, float ior) {
  float r0 = (1 - ior) / (1 + ior);
  r0 = r0 * r0;
  return r0 + (1 - r0) * pow(1 - cosine, 5);
}

static Vec3 refract(Vec3 out_dir, float index_of_refraction, bool& was_internal) {

    // TODO (PathTracer): Task 5
    // Use Snell's Law to refract out_dir through the surface.
    // Return the refracted direction. Set was_internal to true if
    // refraction does not occur due to total internal reflection,
    // and false otherwise.

    // When dot(out_dir,normal=(0,1,0)) is positive, then out_dir corresponds to a
    // ray exiting the surface into vaccum (ior = 1). However, note that
    // you should actually treat this case as _entering_ the surface, because
    // you want to compute the 'input' direction that would cause this output,
    // and to do so you can simply find the direction that out_dir would refract
    // _to_, as refraction is symmetric.
    Vec3 normal(0,1,0); 
    
    //to test: just let cos_t = dot(out_dir,normal)
    float cos_t = out_dir.y; 
    float sin_t = out_dir.x;//(-1.0f)* reflect??

    // float sin_i = sqrt(1.0f - pow(cos_i,2));
    // float sin_t = sin_i / index_of_refraction;
    // float cos_t = sqrt(1.0f - pow(sin_t,2));
    
    Vec3 indir;
    
    
    float ni,nt;
    if(cos_t> 0){//cos theta_t >0
        //out_dir corresponds to a ray exiting the surface into vaccum (ior = 1)
        //But "entering"
        nt = index_of_refraction;
        ni = 1.0f;
    }
    else{
        nt = 1.0f;
        ni = index_of_refraction;
        
    }
    float sin_i = (nt/ni)*sin_t;
    float cos_i = (float)sqrt(1.0f - pow(sin_i,2));
    was_internal = ((1.0f - pow(ni/nt,2))*(1-pow(cos_i,2)) < 0);//was_internal

    indir.x = (-1.0f)*sin_i;
    indir.y = std::abs(cos_t);
    indir.z = (-1.0f)*sin_i;
    return indir;
}


//Given out_dir:out going(eye and point), generates a ramdom sample for in_dir:incoming(last scatter and point)
//Lambertian: ideal diffuse, outgoing direction are uniformly distributed
Scatter BSDF_Lambertian::scatter(Vec3 out_dir) const {

    // TODO (PathTracer): Task 4

    // Sample the BSDF distribution using the cosine-weighted hemisphere sampler.
    // You can use BSDF_Lambertian::evaluate() to compute attenuation.

    // Scatter ret;
    // ret.direction = Vec3{};
    // ret.attenuation = Spectrum{};
    // return ret;

    Scatter ret;    
    

    //outgoing ray
    ret.direction = sampler.sample();//in bsdf.h has the hemi sampler
    ret.attenuation = evaluate(out_dir, ret.direction);
    return ret;
}

//the ratio of incoming to outgoing radiance = fr a.k.a BRDF * cos theta?
//only defined for continuous BSDFs.
Spectrum BSDF_Lambertian::evaluate(Vec3 out_dir, Vec3 in_dir) const {

    // TODO (PathTracer): Task 4

    // Compute the ratio of reflected/incoming radiance when light from in_dir is reflected through out_dir: albedo * cos(theta).

    //get outgoing angle theta
    //no need for in_dir
    
    Vec3 n(0,1,0);
    float theta = dot(out_dir.unit(), n);//error: .unit() // maybe try in_dir next time
    Spectrum res = albedo * std::cos(theta);//because albedo / PI already
    return res;
}

float BSDF_Lambertian::pdf(Vec3 out_dir, Vec3 in_dir) const {

    // TODO (PathTracer): Task 4

    // Compute the PDF for sampling in_dir from the cosine-weighted hemisphere distribution.
    float theta = dot(out_dir, Vec3(0,1,0)); 
    float cos_theta = std::cos(theta);
    cos_theta = clamp(cos_theta,0.0f,1.0f); //error : clamp
    float pdf = cos_theta/ PI_F; //make sure positive?
    return pdf;
}

Scatter BSDF_Mirror::scatter(Vec3 out_dir) const {

    // TODO (PathTracer): Task 5
    Scatter ret;
    ret.direction = reflect(out_dir);
    ret.attenuation = reflectance;
    return ret;
}

Scatter BSDF_Glass::scatter(Vec3 out_dir) const {

    // TODO (PathTracer): Task 5
    Scatter ret;
    // (1) Compute Fresnel coefficient. Tip: Schlick's approximation.
    // Vec3 normal(0,1,0);
    // float cos_i = dot(out_dir,normal);
    

    
    // (2) Reflect or refract probabilistically based on Fresnel coefficient. Tip: RNG::coin_flip
    bool was_internal = false;
    Vec3 refr = refract(out_dir,index_of_refraction, was_internal);
    //if(was_internal):reflect
    float fresnel = (float)Schlick_Approximation(std::abs(refr.y), index_of_refraction);
    if(RNG::coin_flip(fresnel) || was_internal){
        //reflect
        ret.direction = reflect(out_dir);
        ret.attenuation = reflectance;
    }
    else{
        // refract/transmit
        ret.direction = refr;
        float ratio = (refr.y > 0)? (1.0f/index_of_refraction):index_of_refraction;
        ret.attenuation = transmittance * (float)pow(ratio,2);
    }
    // (3) Compute attenuation based on reflectance or transmittance

    // Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?
    // What happens upon total internal reflection?

    
    // ret.direction = Vec3();
    // ret.attenuation = Spectrum{};
    return ret;
}

Scatter BSDF_Refract::scatter(Vec3 out_dir) const {

    // OPTIONAL (PathTracer): Task 5

    // When debugging BSDF_Glass, it may be useful to compare to a pure-refraction BSDF

    Scatter ret;
    ret.direction = Vec3();
    ret.attenuation = Spectrum{};
    return ret;
}

Spectrum BSDF_Diffuse::emissive() const {
    return radiance;
}

} // namespace PT
