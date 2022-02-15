#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors (JUST AS A PLACEHOLDER)
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 6: Implement nearest neighbour interpolation
  
  // return magenta for invalid level
  if( level < 0 || level >= tex.mipmap.size()) return Color(1,0,1,1);//this is magenta
  
  //u, v belongs to [0,1], map it to texture's width and height
  //try clamp the edges
  int sx = (int)floor(clamp(u, 0.0f, 0.9999f) * tex.mipmap[level].width);
  int sy = (int)floor(clamp(v, 0.0f, 0.9999f) * tex.mipmap[level].height);
  
  float r = tex.mipmap[level].texels[4 * (sx + sy * tex.mipmap[level].width)] / 255.0;//error forget/255.0
  float g = tex.mipmap[level].texels[4 * (sx + sy * tex.mipmap[level].width) + 1] / 255.0;
  float b = tex.mipmap[level].texels[4 * (sx + sy * tex.mipmap[level].width) + 2] / 255.0;
  float a = tex.mipmap[level].texels[4 * (sx + sy * tex.mipmap[level].width) + 3] / 255.0;
  return Color(r,g,b,a);
}

template<typename T>
inline Color lerpColor(T ratio, Color start, Color ends){
  return start + ratio * ( ends + (-1)*start);
}

inline Color GetColorFromTexure(Texture& tex, const int& level, const int& x, const int& y){
  float r = tex.mipmap[level].texels[4 * (x + y * tex.mipmap[level].width)]/ 255.0;
  float g = tex.mipmap[level].texels[4 * (x + y * tex.mipmap[level].width) + 1]/ 255.0;
  float b = tex.mipmap[level].texels[4 * (x + y * tex.mipmap[level].width) + 2]/ 255.0;
  float a = tex.mipmap[level].texels[4 * (x + y * tex.mipmap[level].width) + 3]/ 255.0;
  return Color(r,g,b,a);
}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 6: Implement bilinear filtering

  // return magenta for invalid level
    if (level < 0 || level >= tex.mipmap.size()) return Color(1, 0, 1, 1);
    float su = clamp(u, 0.0f, 0.99999f) * tex.mipmap[level].width;
    float sv = clamp(v, 0.0f, 0.99999f) * tex.mipmap[level].height;
     
    float u0 = floor(su)+0.5f; float v0 = floor(sv) + 0.5f;
    //Need other 3 near texels
    //if u,v is on the edge??
    float u1, v1;
    if (round(su) == (int)su) {//su-(int)su <0.5
        u1 = clamp<float>(floor(su) - 1.0 + 0.5f, 0.0, tex.mipmap[level].width);//if out of range : constrain
        swap(u1, u0);
    }
    else { u1 = clamp<float>(ceil(su)+0.5f, 0.0, tex.mipmap[level].width); }
    if (round(sv) == floor(sv) - 1.0 + 0.5f) {
        v1 = clamp<float>((int)sv - 1, 0.0, tex.mipmap[level].height);
        swap(v1, v0);
    }
    else { v1 = clamp<float>(ceil(sv) + 0.5f, 0, tex.mipmap[level].height); }
    //Error maybe I should use texel's center?
    


    //get 4 coordinate: (u0,v0)(u0,v1)(u1,v0)(u1,v1)
    // Color lerpHorizontal1 = lerpColor((su-u0)/(u1-u0), GetColorFromTexure(tex, level, u0, v0), GetColorFromTexure(tex, level, u1, v0));
    // Color lerpHorizontal2 = lerpColor((su-u0)/(u1-u0), GetColorFromTexure(tex, level, u0, v1), GetColorFromTexure(tex, level, u1, v1));
    // Color lerpVertical = lerpColor((sv-v0)/(v1-v0), lerpHorizontal1, lerpHorizontal2);
    
    //try use the bilinear identical 
    float t = (su-u0)/(u1-u0), s = (sv-v0)/(v1-v0);
    Color bilinear = (1-t)*((1-s)*GetColorFromTexure(tex, level, u0, v0)+s*GetColorFromTexure(tex, level, u1, v0))
      +t*((1-s)*GetColorFromTexure(tex, level, u0, v1)+s*GetColorFromTexure(tex, level, u1, v1));

    return bilinear;
}



Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Task 7: Implement trilinear filtering

  // return magenta for invalid level
    //if (level < 0 || level >= tex.mipmap.size()) return Color(1, 0, 1, 1);//this is magenta


}

} // namespace CMU462
