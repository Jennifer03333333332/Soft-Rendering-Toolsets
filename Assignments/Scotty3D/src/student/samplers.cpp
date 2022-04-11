
#include "../rays/samplers.h"
#include "../util/rand.h"

namespace Samplers {

Vec2 Rect::sample() const {

    // TODO (PathTracer): Task 1

    // Generate a uniformly random point on a rectangle of size size.x * size.y
    // Tip: RNG::unit()
    //how to guarantee there's no repeat? RNG::unit() did it?
    return Vec2{RNG::unit()*size.x, RNG::unit()*size.y};
}

Vec3 Sphere::Uniform::sample() const {

    // TODO (PathTracer): Task 7

    // Generate a uniformly random point on the unit sphere.
    // Tip: start with Hemisphere::Uniform
    Hemisphere::Uniform u;
    Vec3 hemi_uni = u.sample();
    return hemi_uni;
}
//Calculate index of _pdf, _cdf
inline size_t image_coor_to_index(size_t row,size_t col,size_t w){
    return row * w + col;
}
inline Vec2 index_to_image_coor(size_t index,size_t w,size_t h){
    Vec2 xy;//col-row
    xy.x = (float)(int)(index % w);//col
    xy.y = (float)(int)(index / w);//row
    return xy;
}
Sphere::Image::Image(const HDR_Image& image) {

    // TODO (PathTracer): Task 7

    // Set up importance sampling data structures for a spherical environment map image.
    // You may make use of the _pdf, _cdf, and total members, or create your own.
    const auto [_w, _h] = image.dimension();
    w = _w;
    h = _h;
    _pdf.resize(w*h,0);
    _cdf.resize(w*h,0);

    for(size_t row = 0;row<h;row++){// y
        for(size_t col = 0;col<w;col++){// x
            //Problem: which theta? the one in a texture map[pi, 0], or the one we get by calculation [0,pi]? 
            float theta = PI_F - PI_F * (float)row / (float)h;  //try PI_F -  then
            size_t index = image_coor_to_index(row,col,w);
            float Lsintheta = image.at(col,row).luma() * std::abs(std::sin(theta));

            _pdf[index] = Lsintheta;
            total+= Lsintheta;
            _cdf[index] = total;
        }
    }

    for(auto& p: _pdf){
        p /= total;
    }
    for(auto& c: _cdf){
        c /= total;
    }
    // float pre_pdf = 0;
    // for(size_t i = 0;i < w*h;i++){
    //     _pdf[i] = _pdf[i] / total;//normalize
    //     _cdf[i] = _pdf[i] + pre_pdf;
    //     pre_pdf += _pdf[i];
    // }
    //row-major
    //_pdf, _cdf, and total

    //Luminosity*sin theta, theta nees to be the center of pixel

}

Vec3 Sphere::Image::sample() const {//importance sampler

    // TODO (PathTracer): Task 7

    // Use your importance sampling data structure to generate a sample direction.
    // Tip: std::upper_bound
    // Generate a weighted random phi and theta coordinate pair, convert to xyz coordinates (Vec3) 
    float unit_sample = RNG::unit();
    
    auto upper = std::upper_bound(_cdf.begin(),_cdf.end(),unit_sample);
    size_t index = upper - _cdf.begin();
    Vec2 col_row = index_to_image_coor(index,w,h);
    // for(size_t i = 0;i < w*h ; i++){
    // //for(size_t i = w*h - 1;i >= 0; i--){//should use some algorhitm, nvm
    //     if(_cdf[i] > unit_sample){
    //         row_col = index_to_image_coor(i,w,h);
    //         break;
    //     }
    // }
    //now I have x,y(in screen coords), how to get xyz?
    //Problem again: should phi in [0,2pi], or [-pi,pi]?
    Vec3 dir;
    float theta = PI_F - PI_F * col_row.y / (float)h;//Change 2
    float phi = (2.0f*PI_F) * col_row.x / (float)w;
    //if(phi >= PI_F)phi = phi - 2.f * PI_F;

    theta = std::clamp(theta, 0.f, PI_F); 
    phi = std::clamp(phi, 0.f, 2.0f*PI_F);  

    dir.x = sin(theta)*cos(phi);
    dir.y = cos(theta);
    dir.z = sin(theta)*sin(phi);
    return dir;
}

float Sphere::Image::pdf(Vec3 dir) const {

    // TODO (PathTracer): Task 7

    // What is the PDF of this distribution at a particular direction?
    //use Jacobian
    //float r = dir.norm();// Is dir unit-vector? Yes
    float theta = PI_F - std::acos(dir.y );/// r
    float phi = std::atan2(dir.z, dir.x); 
    // You should bi-linearly interpolate the value between the 4 nearest pixels.
    // theta [0, Π], phi [0,2Π] 
    if (phi < 0) phi = phi + 2.f * PI_F;
    theta = std::clamp(theta / PI_F,0.f,1.f); 
    phi = std::clamp(phi / (2.0f * PI_F),0.f,1.f);  
    size_t y = (size_t)(theta * (float)h);//h, u is y  Get the cood y, now I need the index y
    size_t x = (size_t)(phi * (float)w);//w, v is x
    if(y > h - 1)y = h - 1;
    if(x > w - 1)x = w - 1;
    //the right-top index should be (w-1,h-1)
    float Jacobian = (float)w * (float)h / (2.f * PI_F * PI_F *sin(theta));
    return Jacobian * _pdf[image_coor_to_index(y,x,w)];
}

Vec3 Point::sample() const {
    return point;
}

Vec3 Triangle::sample() const {
    float u = std::sqrt(RNG::unit());
    float v = RNG::unit();
    float a = u * (1.0f - v);
    float b = u * v;
    return a * v0 + b * v1 + (1.0f - a - b) * v2;
}

Vec3 Hemisphere::Uniform::sample() const {

    float Xi1 = RNG::unit();
    float Xi2 = RNG::unit();

    float theta = std::acos(Xi1);
    float phi = 2.0f * PI_F * Xi2;

    float xs = std::sin(theta) * std::cos(phi);
    float ys = std::cos(theta);
    float zs = std::sin(theta) * std::sin(phi);

    return Vec3(xs, ys, zs);
}

Vec3 Hemisphere::Cosine::sample() const {

    float phi = RNG::unit() * 2.0f * PI_F;
    float cos_t = std::sqrt(RNG::unit());

    float sin_t = std::sqrt(1 - cos_t * cos_t);
    float x = std::cos(phi) * sin_t;
    float z = std::sin(phi) * sin_t;
    float y = cos_t;

    return Vec3(x, y, z);
}

} // namespace Samplers
