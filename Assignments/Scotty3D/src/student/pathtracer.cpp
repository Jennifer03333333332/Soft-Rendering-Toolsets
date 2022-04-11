
#include "../rays/pathtracer.h"
#include "../rays/samplers.h"
#include "../util/rand.h"
#include "debug.h"

#define TASK_4 0
#define TASK_6 0

namespace PT {


//x,y: screen space coordinates [0,w],[0,h]
Spectrum Pathtracer::trace_pixel(size_t x, size_t y) {

    // TODO (PathTracer): Task 1

    // Generate a ray that uniformly samples pixel (x,y) and return the incoming light.
    // The following code generates a ray at the bottom left of the pixel every time.
    // You'll need to change this so that it gets a new location each time trace_pixel is called for part 3

    // Tip: Use Rect::sample to get a new location each time trace_pixel is called
    // Tip: log_ray is useful for debugging

    Vec2 xy((float)x, (float)y);
    Vec2 wh((float)out_w, (float)out_h);


    Samplers::Rect sampler( Vec2(1.0f,1.0f) );
    Vec2 sample_pos = sampler.sample();//[0,1],[0,1]


    Ray ray = camera.generate_ray((xy + sample_pos) / wh);
    ray.depth = max_depth;
    //debug: 10.0f means show the ray until time become 10.0f
    //if(RNG::coin_flip(0.0005f)) log_ray(ray, 10.0f);
    // Pathtracer::trace() returns the incoming light split into emissive and reflected components.
    auto [emissive, reflected] = trace(ray);
    return emissive + reflected;
}

Spectrum Pathtracer::sample_indirect_lighting(const Shading_Info& hit) {
    // TODO (PathTrace): Task 4

    // This function computes a single-sample Monte Carlo estimate of the _indirect_
    // lighting at our ray intersection point.
    Spectrum radiance;
    // (1) Randomly sample a new ray direction from the BSDF distribution using BSDF::scatter().
    Scatter input = hit.bsdf.scatter(hit.out_dir) ;//
    // (2) Create a new world-space ray and call Pathtracer::trace() to get incoming light. You
    // should modify time_bounds so that the ray does not intersect at time = 0. Remember to
    // set the new depth value.
    Vec3 world_in_dir = hit.object_to_world.rotate(input.direction);//input.direction;//
    //// * error // Mat4::rotate_to(result.normal);
    Ray world_ray(hit.pos,world_in_dir,Vec2(EPS_F,FLT_MAX),hit.depth - 1);//hit.depth - 1 //max_depth

    Spectrum indirect_light = trace(world_ray).second;
    //world_ray.depth--;
    // (3) Add contribution due to incoming light scaled by BSDF attenuation. Whether you
    // compute the BSDF scattering PDF should depend on if the BSDF is a discrete distribution
    // (see BSDF::is_discrete()).
    

    if(hit.bsdf.is_discrete()){// 
        indirect_light *= input.attenuation;//*input.direction;
    }
    else{//continuous: divided by pdf
        float pdf = hit.bsdf.pdf(hit.out_dir,input.direction);
        indirect_light = indirect_light*input.attenuation*(1.0f/pdf); 
    }
    // You should only use the indirect component of incoming light (the second value returned
    // by Pathtracer::trace()), as the direct component will be computed in
    // Pathtracer::sample_direct_lighting().
    radiance += indirect_light;
    return radiance;
}

Spectrum Pathtracer::sample_direct_lighting(const Shading_Info& hit) {

    // This function computes a Monte Carlo estimate of the _direct_ lighting at our ray
    // intersection point by sampling both the BSDF and area lights.

    // Point lights are handled separately, as they cannot be intersected by tracing rays
    // into the scene.
    Spectrum radiance = point_lighting(hit);

    Scatter input = hit.bsdf.scatter(hit.out_dir);
    Vec3 world_in_dir = hit.object_to_world.rotate(input.direction);

    Ray world_ray(hit.pos,world_in_dir,Vec2(EPS_F,FLT_MAX),0);

    Spectrum direct_light = trace(world_ray).first;
    
    float pdf;
    if(hit.bsdf.is_discrete()){
        direct_light *= input.attenuation;
    }
    else{//continuous: divided by pdf
        pdf = hit.bsdf.pdf(hit.out_dir,input.direction);
        direct_light = direct_light*input.attenuation*(1.0f/pdf); 
    }

    radiance += direct_light;
    // TODO (PathTrace): Task 4

    // For task 4, this function should perform almost the same sampling procedure as
    // Pathtracer::sample_indirect_lighting(), but instead accumulates the emissive component of
    // incoming light (the first value returned by Pathtracer::trace()). Note that since we only
    // want emissive, we can trace a ray with depth = 0.

#if TASK_4 == 1
    return radiance;
#endif

    // TODO (PathTrace): Task 6

    // For task 6, we want to upgrade our direct light sampling procedure to also
    // sample area lights using mixture sampling.

    // (1) If the BSDF is discrete, we don't need to bother sampling lights: the behavior
    // should be the same as task 4.
    if(hit.bsdf.is_discrete()){
        return radiance;
    }
    radiance = radiance - direct_light;//?
    // (2) Otherwise, we should randomly choose whether we get our sample from `BSDF::scatter`
    // or `Pathtracer::sample_area_lights`. Note that `Pathtracer::sample_area_lights` returns
    // a world-space direction pointing toward an area light. Choose between the strategies 
    // with equal probability.
    Vec3 random_in_dir;

    Vec3 sample_arealight_inworld = sample_area_lights(hit.pos);//hit.pos ->an area light
    if(RNG::coin_flip(0.5)){
        //use BSDF::Scatter
        random_in_dir = world_in_dir;//world_in_dir is from the task 4's: in_dir in world space
    }else{
        //use sample_area_light
        random_in_dir = sample_arealight_inworld;
    }

    // (3) Create a new world-space ray and call Pathtracer::trace() to get incoming light. You
    // should modify time_bounds so that the ray does not intersect at time = 0. We are again
    // only interested in the emissive component, so the ray depth can be zero.
    Ray world_ray_task6(hit.pos,random_in_dir,Vec2(EPS_F,FLT_MAX),0);

    //if(RNG::coin_flip(0.00005f)) log_ray(world_ray_task6, 5.0f);

    // (4) Add estimate of incoming light scaled by BSDF attenuation. Given a sample,
    // we don't know whether it came from the BSDF or the light, so you should use BSDF::evaluate(),
    // BSDF::pdf(), and Pathtracer::area_lights_pdf() to compute the proper weighting.
    // What is the PDF of our sample, given it could have been produced from either source?
    Spectrum direct_light_task6 = trace(world_ray_task6).first;
    //weighted pdf
    //world space
    float pdf_area_light = area_lights_pdf(hit.pos,sample_arealight_inworld);//area_lights_pdf(hit.pos,world_in_dir);??random_indir hit.world_to_object.rotate(sample_arealight_inworld)
    //local space
    float pdf_task4 = hit.bsdf.pdf(hit.out_dir, hit.world_to_object.rotate(world_in_dir));//bsdf:local
    pdf = (pdf_task4 + pdf_area_light)/2.0f;
    //evaluate: local space
    Spectrum attenuation_task6 = hit.bsdf.evaluate(hit.out_dir,hit.world_to_object.rotate(random_in_dir));//input.attenuation;
    
    direct_light_task6 = direct_light_task6*attenuation_task6*(1.0f/pdf); 
    
    radiance += direct_light_task6;
#if TASK_6 == 1
    return radiance;
#endif

    return radiance;
}
//given a ray, return emissive and reflect
std::pair<Spectrum, Spectrum> Pathtracer::trace(const Ray& ray) {

    // This function orchestrates the path tracing process. For convenience, it
    // returns the incoming light along a ray in two components: emitted from the
    // surface the ray hits, and reflected through that point from other sources.

    // Trace ray into scene.
    Trace result = scene.hit(ray);
    if(!result.hit) {

        // If no surfaces were hit, sample the environemnt map.
        if(env_light.has_value()) {
            return {env_light.value().evaluate(ray.dir), {}};
        }
        return {};
    }

    // If we're using a two-sided material, treat back-faces the same as front-faces
    const BSDF& bsdf = materials[result.material];
    if(!bsdf.is_sided() && dot(result.normal, ray.dir) > 0.0f) {
        result.normal = -result.normal;
    }

    // TODO (PathTracer): Task 4
    // You will want to change the default normal_colors in debug.h, or delete this early out.
    if(debug_data.normal_colors) return {Spectrum::direction(result.normal), {}};

    // If the BSDF is emissive, stop tracing and return the emitted light
    Spectrum emissive = bsdf.emissive();
    if(emissive.luma() > 0.0f) return {emissive, {}};

    // If the ray has reached maximum depth, stop tracing
    if(ray.depth == 0) return {};

    // Set up shading information
    Mat4 object_to_world = Mat4::rotate_to(result.normal);
    Mat4 world_to_object = object_to_world.T();
    Vec3 out_dir = world_to_object.rotate(ray.point - result.position).unit();

    Shading_Info hit = {bsdf,    world_to_object, object_to_world, result.position,
                        out_dir, result.normal,   ray.depth};

    // Sample and return light reflected through the intersection
    return {emissive, sample_direct_lighting(hit) + sample_indirect_lighting(hit)};// Spectrum{}
}

} // namespace PT
