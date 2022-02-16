#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

    //self-defined functions
    template<typename T>
    inline Color lerpColor(T ratio, Color start, Color ends) {
        //return start + ratio * ( ends + (-1)*start);
        return (1 - ratio) * start + ratio * ends;
    }
    //x,y are texture's coordinates
    inline Color GetColorFromTexture(Texture& tex, const int& level, const int& x, const int& y) {
        float r = tex.mipmap[level].texels[4 * (x + y * tex.mipmap[level].width)] / 255.0f;
        float g = tex.mipmap[level].texels[4 * (x + y * tex.mipmap[level].width) + 1] / 255.0f;
        float b = tex.mipmap[level].texels[4 * (x + y * tex.mipmap[level].width) + 2] / 255.0f;
        float a = tex.mipmap[level].texels[4 * (x + y * tex.mipmap[level].width) + 3] / 255.0f;
        return Color(r, g, b, a);
    }
    inline void SetColorToTexture(Texture& tex, const int& level, const int& x, const int& y, const Color &c) {
        tex.mipmap[level].texels[4 * (x + y * tex.mipmap[level].width)] = c.r*255.0f;
        tex.mipmap[level].texels[4 * (x + y * tex.mipmap[level].width) + 1] = c.g*255.0f;
        tex.mipmap[level].texels[4 * (x + y * tex.mipmap[level].width) + 2] = c.b*255.0f;
        tex.mipmap[level].texels[4 * (x + y * tex.mipmap[level].width) + 3] = c.a*255.0f;
        return;
    }
    inline Color average4Colors() {

    }

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
  //Initialize each mipmap level
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
  //Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {//for each level of mipmap
    MipLevel& mip = tex.mipmap[i];
    //for every texel in this level of mipmap
    for (int x = 0; x < mip.width; x++) {
        for (int y = 0; y < mip.height; y++) {
            Color sum = Color(0, 0, 0, 0);
            //Concentrate from the i-1 level. 4 -> 1
            for (int m = 0; m < 2; m++) {
                for (int n = 0; n < 2; n++) {
                    Color c = GetColorFromTexture(tex, i - 1, 2 * x + m, 2 * y + n);
                    sum += Color(c.r * c.a, c.g * c.a, c.b * c.a, c.a);//sum += c;
                }
            }
            sum *= 0.25f;
             //if (sum.a != 0) {
             //    sum.r /= sum.a;
             //    sum.g /= sum.a;
             //    sum.b /= sum.a;
             //}
            //SetColorToTexture(tex, i, x, y, sum);
            float_to_uint8(&mip.texels[4 * (x + y * mip.width)], &sum.r);
        }
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



//u,v are [0,1] coordinates
Color Sampler2DImp::sample_bilinear(Texture& tex, float u, float v, int level) {
    if (level < 0 || level >= tex.mipmap.size()) return Color(1, 0, 1, 1);
    float su = clamp(u, 0.0f, 0.99999f) * tex.mipmap[level].width;//clamp the edges
    float sv = clamp(v, 0.0f, 0.99999f) * tex.mipmap[level].height;
    //Calculate the nearest 4 texel's center:
    float u0 = floor(su) + 0.5f; float v0 = floor(sv) + 0.5f;
    float u1, v1;
    if (su - (int)su < 0.5f) {//
        u1 = clamp<float>(u0 - 1, 0.0, tex.mipmap[level].width);//if out of range : constrain
        swap(u1, u0);
    }
    else { u1 = clamp<float>(u0 + 1, 0.0, tex.mipmap[level].width); }
    if (sv - (int)sv < 0.5f) {
        v1 = clamp<float>(v0 - 1, 0.0, tex.mipmap[level].height);
        swap(v1, v0);
    }
    else { v1 = clamp<float>(v0+1, 0, tex.mipmap[level].height); }
    //Bilinear interpolate the color
    Color lerpHorizontal1 = lerpColor((su - u0) / (u1 - u0), GetColorFromTexture(tex, level, u0, v0), GetColorFromTexture(tex, level, u1, v0));
    Color lerpHorizontal2 = lerpColor((su - u0) / (u1 - u0), GetColorFromTexture(tex, level, u0, v1), GetColorFromTexture(tex, level, u1, v1));
    Color lerpVertical = lerpColor((sv - v0) / (v1 - v0), lerpHorizontal1, lerpHorizontal2);
    return lerpVertical;
}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Task 7: Implement trilinear filtering

  // return magenta for invalid level
  // 1 Calculate the level. u_scale v_scale means dx, dy; tex.width, tex.height means du,dv
    float Lsquare_x = pow(tex.width / u_scale,2) + pow(tex.height / u_scale, 2);
    float Lsquare_y = pow(tex.width / v_scale, 2) + pow(tex.height / v_scale, 2);
    //float L = max(tex.width/u_scale,tex.height/v_scale);
    float level = log2f(sqrt(max(Lsquare_x,Lsquare_y)));
    if (level < 0) level = 0.0f;
    if (level >= tex.mipmap.size()) return Color(1, 0, 1, 1);//this is magenta
    //Trilinear interpolation
    int lowLevel = (int)floor(level);
    int highLevel = lowLevel + 1;
    if(highLevel>= tex.mipmap.size())return sample_bilinear(tex, u, v, tex.mipmap.size()-1);

    Color lowLevelColor = sample_bilinear(tex, u, v, lowLevel);
    Color highLevelColor = sample_bilinear(tex, u, v, highLevel);
    return lerpColor((level - (int)level), lowLevelColor, highLevelColor);
}

} // namespace CMU462
