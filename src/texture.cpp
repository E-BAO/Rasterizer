#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>

namespace CGL {

    Color MipLevel::get_texel(int tx, int ty){
        if(tx < 0 || tx >= width || ty < 0 || ty >= height)
            return Color(1,1,1);
        unsigned char rgb[3];
        int offset = (ty * width + tx) * 3 ;
        if(offset + 2 >= texels.size())
            cout<<"error: over texels size = "<<offset<<endl;
        rgb[0] = texels[offset + 0];
        rgb[1] = texels[offset + 1];
        rgb[2] = texels[offset + 2];
        return Color(rgb);
    }
    
Color Texture::sample(const SampleParams &sp) {
  // Parts 5 and 6: Fill this in.
  // Should return a color sampled based on the psm and lsm parameters given
    if(sp.lsm == L_ZERO){
        if(sp.psm == P_NEAREST)
            return sample_nearest(sp.p_uv, 0);
        else
            return sample_bilinear(sp.p_uv, 0);
    }else if(sp.lsm == L_NEAREST){
        int level = get_level(sp);
        if(level >= mipmap.size())
            level = mipmap.size() - 1;
        if(level < 0)
            level = 0;
        if(sp.psm == P_NEAREST)
            return sample_nearest(sp.p_uv, level);
        else
            return sample_bilinear(sp.p_uv, level);
    }else if(sp.lsm == L_LINEAR){
        double levelD = get_level(sp);
//        cout<<"levelD"<<levelD<<endl;
        int level0 = floor(levelD);
        int level1 = level0 + 1;

        if(level0 >= mipmap.size())
            level0 = max((double)mipmap.size() - 2.0,0.0);
        if(level0 < 0)
            level0 = 0;
    
        double s = (double)level1 - levelD;
        
        if(level1 >= mipmap.size() || level1 < 0){
            if(sp.psm == P_NEAREST)
                return sample_nearest(sp.p_uv, level0);
            else  //trilinear
                return sample_bilinear(sp.p_uv, level0);
        }
        
        if(level1 >= mipmap.size() || level1 < 0)
            cout<<"error level1"<<endl;
        if(level0 >= mipmap.size() || level0 < 0)
            cout<<"error level0"<<endl;
        
        if(sp.psm == P_NEAREST){
//            Color c1 = sample_nearest(sp.p_uv, level0);
            Color c2 = s * sample_nearest(sp.p_uv, level0) + (1.0 - s) * sample_nearest(sp.p_uv, level1);
            if(c2.r <=1e-3 && c2.g <=1e-3 && c2.b <=1e-3){
                cout<<"error black"<<endl;
                c2.g = 1;
            }
            return c2;
        }else{//trilinear
            Color c1 = sample_bilinear(sp.p_uv, level0);
            Color c2 = sample_bilinear(sp.p_uv, level1);
            Color c3 = s * c1 + (1 - s) *c2;
            if(c1.r <=1e-3 && c1.g <=1e-3 && c1.b <=1e-3){
                c3.g = 1;
                cout<<mipmap.size()<<level0 << " c1 error black"<<endl;
            }
            if(c2.r <=1e-3 && c2.g <=1e-3 && c2.b <=1e-3){
                c3.b = 1;
                cout<<mipmap.size()<<level1 << " c2 error black"<<endl;
            }
            
            return c3;
        }
    }
    cout<<"error exit"<<endl;
    return Color(0,0,0);
}

float Texture::get_level(const SampleParams &sp) {
  // Optional helper function for Parts 5 and 6
    Vector2D dx_uv = sp.p_dx_uv - sp.p_uv;
    Vector2D dy_uv = sp.p_dy_uv - sp.p_uv;
    dx_uv.x *= width;
    dy_uv.x *= width;
    dx_uv.y *= height;
    dy_uv.y *= height;
    double L = max(sqrt(pow(dx_uv.x, 2.0) + pow(dx_uv.y, 2.0)), sqrt(pow(dy_uv.x, 2.0) + pow(dy_uv.y, 2.0)));
    return log2(L);
}

// Returns the nearest sample given a particular level and set of uv coords
Color Texture::sample_nearest(Vector2D uv, int level) {
  // Optional helper function for Parts 5 and 6
  // Feel free to ignore or create your own
//    if(level >= mipmap.size())
//        level = mipmap.size() - 1;
//    else if(level < 0)
//        level = 0;
    level = max(0, min((int)mipmap.size() - 1, level));
    MipLevel * ml = &mipmap[level];
    int tx = uv.x * ml->width;
    int ty = uv.y * ml->height;
    return ml->get_texel(tx,ty);
}
    
// Returns the bilinear sample given a particular level and set of uv coords
Color Texture::sample_bilinear(Vector2D uv, int level) {
  // Optional helper function for Parts 5 and 6
  // Feel free to ignore or create your own
    if(level >= mipmap.size())
        level = mipmap.size() - 1;
    else if(level < 0)
        level = 0;
    MipLevel * ml = &mipmap[level];
    double posX = uv.x * ml->width;
    double posY = uv.y * ml->height;
    
    double gridX = floor(posX);
    double gridY = floor(posY);
    
    int x0, y0, x1, y1;
    if(posX - gridX < 0.5){
        x1 = gridX;
        x0 = x1 - 1;
    }else{
        x0 = gridX;
        x1 = x0 + 1;
    }
    double s =  posX - 0.5 - (double)x0;

    if(posY - gridY < 0.5){
        y1 = gridY;
        y0 = y1 - 1;
    }else{
        y0 = gridY;
        y1 = y0 + 1;
    }
    double t = posY - 0.5 - (double)y0;
    
    if(t > 1 || s > 1)
        cout<<"error  t = "<<t<<"  s = "<<s<<endl;
    
    Color c0 = ml->get_texel(x0, y0) + s * (ml->get_texel(x1, y0) - ml->get_texel(x0, y0));
    Color c1 = ml->get_texel(x0, y1) + s * (ml->get_texel(x1, y1) - ml->get_texel(x0, y1));
    
    return c0 + t * (c1 - c0);
}




/****************************************************************************/



inline void uint8_to_float(float dst[3], unsigned char *src) {
  uint8_t *src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
}

inline void float_to_uint8(unsigned char *dst, float src[3]) {
  uint8_t *dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
  dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
  dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
}

void Texture::generate_mips(int startLevel) {

  // make sure there's a valid texture
  if (startLevel >= mipmap.size()) {
    std::cerr << "Invalid start level";
  }

  // allocate sublevels
  int baseWidth = mipmap[startLevel].width;
  int baseHeight = mipmap[startLevel].height;
  int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  mipmap.resize(startLevel + numSubLevels + 1);

  int width = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel &level = mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width = max(1, width / 2);
    //assert (width > 0);
    height = max(1, height / 2);
    //assert (height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(3 * width * height);
  }

  // create mips
  int subLevels = numSubLevels - (startLevel + 1);
  for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
       mipLevel++) {

    MipLevel &prevLevel = mipmap[mipLevel - 1];
    MipLevel &currLevel = mipmap[mipLevel];

    int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
    int currLevelPitch = currLevel.width * 3; // 32 bit RGB

    unsigned char *prevLevelMem;
    unsigned char *currLevelMem;

    currLevelMem = (unsigned char *)&currLevel.texels[0];
    prevLevelMem = (unsigned char *)&prevLevel.texels[0];

    float wDecimal, wNorm, wWeight[3];
    int wSupport;
    float hDecimal, hNorm, hWeight[3];
    int hSupport;

    float result[3];
    float input[3];

    // conditional differentiates no rounding case from round down case
    if (prevLevel.width & 1) {
      wSupport = 3;
      wDecimal = 1.0f / (float)currLevel.width;
    } else {
      wSupport = 2;
      wDecimal = 0.0f;
    }

    // conditional differentiates no rounding case from round down case
    if (prevLevel.height & 1) {
      hSupport = 3;
      hDecimal = 1.0f / (float)currLevel.height;
    } else {
      hSupport = 2;
      hDecimal = 0.0f;
    }

    wNorm = 1.0f / (2.0f + wDecimal);
    hNorm = 1.0f / (2.0f + hDecimal);

    // case 1: reduction only in horizontal size (vertical size is 1)
    if (currLevel.height == prevLevel.height) {
      //assert (currLevel.height == 1);

      for (int i = 0; i < currLevel.width; i++) {
        wWeight[0] = wNorm * (1.0f - wDecimal * i);
        wWeight[1] = wNorm * 1.0f;
        wWeight[2] = wNorm * wDecimal * (i + 1);

        result[0] = result[1] = result[2] = 0.0f;

        for (int ii = 0; ii < wSupport; ii++) {
          uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
          result[0] += wWeight[ii] * input[0];
          result[1] += wWeight[ii] * input[1];
          result[2] += wWeight[ii] * input[2];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (3 * i), result);
      }

      // case 2: reduction only in vertical size (horizontal size is 1)
    } else if (currLevel.width == prevLevel.width) {
      //assert (currLevel.width == 1);

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        result[0] = result[1] = result[2] = 0.0f;
        for (int jj = 0; jj < hSupport; jj++) {
          uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
          result[0] += hWeight[jj] * input[0];
          result[1] += hWeight[jj] * input[1];
          result[2] += hWeight[jj] * input[2];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (currLevelPitch * j), result);
      }

      // case 3: reduction in both horizontal and vertical size
    } else {

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          // convolve source image with a trapezoidal filter.
          // in the case of no rounding this is just a box filter of width 2.
          // in the general case, the support region is 3x3.
          for (int jj = 0; jj < hSupport; jj++)
            for (int ii = 0; ii < wSupport; ii++) {
              float weight = hWeight[jj] * wWeight[ii];
              uint8_to_float(input, prevLevelMem +
                                        prevLevelPitch * (2 * j + jj) +
                                        3 * (2 * i + ii));
              result[0] += weight * input[0];
              result[1] += weight * input[1];
              result[2] += weight * input[2];
            }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
        }
      }
    }
  }
}

}
