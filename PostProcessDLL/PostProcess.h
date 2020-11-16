﻿#pragma once

#include <cstdint>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <iomanip>
#include <omp.h>
#include <cmath>
#include <algorithm> 

#include "../../HydraAPI/hydra_api/LiteMath.h"

////////////////////////////////////////////////////////////////////////////////////////////////////

using LiteMath::float3;
using LiteMath::float4;


float Min3(const float a, const float b, const float c)
{  
  if      (a <= b && a <= c) return a;
  else if (b <= a && a <= c) return b;
  else                       return c;
}

float Max3(const float a, const float b, const float c)
{
  if      (a >= b && a >= c) return a;
  else if (b >= a && a >= c) return b;
  else                       return c;
}

float Clamp(const float u, const float a, const float b) { const float r = fmax(a, u); return fmin(r, b); }

void ClampMinusToZero(float3& data)
{
  data.x = fmax(data.x, 0.0F);
  data.y = fmax(data.y, 0.0F);
  data.z = fmax(data.z, 0.0F);
}

float Luminance(const float3 data)
{
  return dot(data, float3(0.2126F, 0.7152F, 0.0722F));  
}


float Mean(const float3 a_data)
{
  return (a_data.x + a_data.y + a_data.z) / 3.0F;
}

float Distance(const float x1, const float y1, const float x2, const float y2)
{
  float differenceX = x1 - x2;
  float differenceY = y1 - y2;

  float distance = sqrt(differenceX * differenceX + differenceY * differenceY);
  return distance;
}

float Distance(const float3 point1, const float point2)
{
  float3 difference(point1.x - point2, point1.y - point2, point1.z - point2);

  float distance = sqrt(difference.x * difference.x + difference.y * difference.y + difference.z * difference.z);
  return distance;
}

void SpectralToRGB(float &r, float &g, float &b, const float l) // RGB <0,1> <- lambda l <400,700> [nm]
{
  float t;
  r = 0.0;
  g = 0.0;
  b = 0.0;

  if ((l >= 400.0f) && (l < 410.0f)) { t = (l - 400.0f) / (410.0f - 400.0f); r = +(0.33f * t) - (0.20f * t * t); }
  else if ((l >= 410.0f) && (l < 475.0f)) { t = (l - 410.0f) / (475.0f - 410.0f); r = 0.14f - (0.13f * t * t); }
  else if ((l >= 545.0f) && (l < 595.0f)) { t = (l - 545.0f) / (595.0f - 545.0f); r = +(1.98f * t) - (t * t); }
  else if ((l >= 595.0f) && (l < 650.0f)) { t = (l - 595.0f) / (650.0f - 595.0f); r = 0.98f + (0.06f * t) - (0.40f * t * t); }
  else if ((l >= 650.0f) && (l < 700.0f)) { t = (l - 650.0f) / (700.0f - 650.0f); r = 0.65f - (0.84f * t) + (0.20f * t * t); }
  if ((l >= 415.0f) && (l < 475.0f)) { t = (l - 415.0f) / (475.0f - 415.0f); g = +(0.80f * t * t); }
  else if ((l >= 475.0f) && (l < 590.0f)) { t = (l - 475.0f) / (590.0f - 475.0f); g = 0.8f + (0.76f * t) - (0.80f * t * t); }
  else if ((l >= 585.0f) && (l < 639.0f)) { t = (l - 585.0f) / (639.0f - 585.0f); g = 0.84f - (0.84f * t); }
  if ((l >= 400.0f) && (l < 475.0f)) { t = (l - 400.0f) / (475.0f - 400.0f); b = +(2.20f * t) - (1.50f * t * t); }
  else if ((l >= 475.0f) && (l < 560.0f)) { t = (l - 475.0f) / (560.0f - 475.0f); b = 0.7f - (t)+(0.30f * t * t); }
}

void ConvertRgbToHsv(float3& a_data)
{
  float r = a_data.x;
  float g = a_data.y;
  float b = a_data.z;

  r *= 255.0F;
  g *= 255.0F;
  b *= 255.0F;

  float rgb_min = Min3(r, g, b);
  float rgb_max = Max3(r, g, b);

  float H = 0.0F;
  float S = 0.0F;
  float V = rgb_max;

  // norm. val 1
  if (V != 0.0F)
  {
    r /= V;
    g /= V;
    b /= V;
  }

  rgb_min = Min3(r, g, b);
  rgb_max = Max3(r, g, b);

  S = rgb_max - rgb_min;

  if (S != 0.0F)
  {
    // Saturation 
    r = (r - rgb_min) / S;
    g = (g - rgb_min) / S;
    b = (b - rgb_min) / S;
  }

  //rgb_min = Min3(r, g, b);
  //rgb_max = Max3(r, g, b);

  if (rgb_max == r)
  {
    H = 0.0F + 60.0F*(g - b);
  }
  else if (rgb_max == g)
  {
    H = 120.0F + 60.0F*(b - r);
  }
  else /* rgb_max == *b */
  {
    H = 240.0F + 60.0F*(r - g);
  }

  if (H < 0.0F)   H += 360.0F;
  if (H > 360.0F) H -= 360.0F;

  V /= 255.0F;

  a_data.x = H;
  a_data.y = S;
  a_data.z = V;
}

void ConvertRgbToXyz65(float4* data)
{
  float var_R = data->x; //var_R from 0 to 1
  float var_G = data->y; //var_G from 0 to 1
  float var_B = data->z; //var_B from 0 to 1

  if (var_R > 0.04045f) var_R = pow((var_R + 0.055f) / 1.055f, 2.4f);
  else                  var_R = var_R / 12.92f;

  if (var_G > 0.04045f) var_G = pow((var_G + 0.055f) / 1.055f, 2.4f);
  else                  var_G = var_G / 12.92f;

  if (var_B > 0.04045f) var_B = pow((var_B + 0.055f) / 1.055f, 2.4f);
  else                  var_B = var_B / 12.92f;

  var_R *= 100.0f;
  var_G *= 100.0f;
  var_B *= 100.0f;

  //Chromatic adaptation. Observer. = 2°, Illuminant = D65

  data->x = var_R * 0.4124564f + var_G * 0.3575761f + var_B * 0.1804375f;
  data->y = var_R * 0.2126729f + var_G * 0.7151522f + var_B * 0.0721750f;
  data->z = var_R * 0.0193339f + var_G * 0.1191920f + var_B * 0.9503041f;
}
void ConvertSrgbToXyz(float3& a_data)
{
  const float R = a_data.x;
  const float G = a_data.y;
  const float B = a_data.z;

  a_data.x = R * 0.4124564F + G * 0.3575761F + B * 0.1804375F;
  a_data.y = R * 0.2126729F + G * 0.7151522F + B * 0.0721750F;
  a_data.z = R * 0.0193339F + G * 0.1191920F + B * 0.9503041F;
}


void ConvertXyz65ToRgb(float4* data)
{
  float var_X = data->x / 100.0f;  //X from 0 to  95.047      (Observer = 2°, Illuminant = D65)
  float var_Y = data->y / 100.0f;  //Y from 0 to 100.000
  float var_Z = data->z / 100.0f;  //Z from 0 to 108.883

                                   //Chromatic adaptation.

  float R = var_X * 3.2404542f + var_Y * -1.5371385f + var_Z * -0.4985314f;
  float G = var_X * -0.9692660f + var_Y * 1.8760108f + var_Z * 0.0415560f;
  float B = var_X * 0.0556434f + var_Y * -0.2040259f + var_Z * 1.0572252f;

  if (R > 0.0031308f) R = 1.055f * (pow(R, 0.4166666666f)) - 0.055f; // 1.0f / 2.4f = 0.4166666666f
  else                R = 12.92f * R;
  if (G > 0.0031308f) G = 1.055f * (pow(G, 0.4166666666f)) - 0.055f;
  else                G = 12.92f * G;
  if (B > 0.0031308f) B = 1.055f * (pow(B, 0.4166666666f)) - 0.055f;
  else                B = 12.92f * B;

  data->x = R;
  data->y = G;
  data->z = B;
}

void ConvertXyzToSrgb(float3& a_data)
{
  const float X = a_data.x;
  const float Y = a_data.y;
  const float Z = a_data.z;

  a_data.x = X *  3.2404542F + Y * -1.5371385F + Z * -0.4985314F;
  a_data.y = X * -0.9692660F + Y *  1.8760108F + Z *  0.0415560F;
  a_data.z = X *  0.0556434F + Y * -0.2040259F + Z *  1.0572252F;
}


void ConvertXyzToLab(float3& a_data)
{
  float var_X = a_data.x / 95.0470F;   //Observer= 2°, Illuminant= D65 
  float var_Y = a_data.y / 100.000F;   
  float var_Z = a_data.z / 108.883F;   

  if (var_X > 0.008856F) var_X = pow(var_X, 0.33333333F);           //1.0f / 3.0f    = 0.33333333f
  else                   var_X = (7.787037F * var_X) + 0.13793103F; //16.0f / 116.0f = 0.13793103f
  if (var_Y > 0.008856F) var_Y = pow(var_Y, 0.33333333F);
  else                   var_Y = (7.787037F * var_Y) + 0.13793103F;
  if (var_Z > 0.008856F) var_Z = pow(var_Z, 0.33333333F);
  else                   var_Z = (7.787037F * var_Z) + 0.13793103F;

  a_data.x = 116.0F *  var_Y - 16.0F;  // CIE-L
  a_data.y = 500.0F * (var_X - var_Y); // CIE-a
  a_data.z = 200.0F * (var_Y - var_Z); // CIE-b
}

void ConvertHsvToRgb(float4* data)
{
  float h = data->x;
  float s = data->y;
  float v = data->z;

  float hh, p, q, t, ff, r, g, b;
  int i;

  if (s <= 0.0f) // < is bogus, just shuts up warnings
  {
    r = v;
    g = v;
    b = v;
  }

  hh = h;

  if (hh >= 360.0f) hh = 0.0f;
  if (hh < 0.0f)    hh = 359.9999f;

  hh /= 60.0f;
  i = (int)hh;
  ff = hh - i;
  p = v * (1.0f - s);
  q = v * (1.0f - (s * ff));
  t = v * (1.0f - (s * (1.0f - ff)));

  switch (i)
  {
    case 0:
      r = v;
      g = t;
      b = p;
      break;
    case 1:
      r = q;
      g = v;
      b = p;
      break;
    case 2:
      r = p;
      g = v;
      b = t;
      break;
    case 3:
      r = p;
      g = q;
      b = v;
      break;
    case 4:
      r = t;
      g = p;
      b = v;
      break;
    case 5:
    default:
      r = v;
      g = p;
      b = q;
      break;
  }

  data->x = r;
  data->y = g;
  data->z = b;
}

void ConvertLabToXyz(float4* data)
{
  float var_Y = (data->x + 16.0f) / 116.0f;
  float var_X = var_Y + data->y / 500.0f;
  float var_Z = var_Y - data->z / 200.0f;

  float pow3var_X = var_X * var_X * var_X;
  float pow3var_Y = var_Y * var_Y * var_Y;
  float pow3var_Z = var_Z * var_Z * var_Z;

  if (pow3var_X > 0.008856f) var_X = pow3var_X;
  else                       var_X = (var_X - 0.13793103f) / 7.787f;  //16.0f / 116.0f = 0.13793103f
  if (pow3var_Y > 0.008856f) var_Y = pow3var_Y;
  else                       var_Y = (var_Y - 0.13793103f) / 7.787f;
  if (pow3var_Z > 0.008856f) var_Z = pow3var_Z;
  else                       var_Z = (var_Z - 0.13793103f) / 7.787f;

  data->x = var_X * 95.047f;   // Observer= 2°, Illuminant= D65
  data->y = var_Y * 100.0f;    //
  data->z = var_Z * 108.883f;  //
}

void MatrixCat02(float3& a_data)
{
  // Conversion to cone responses. Mcat02 in CAM02. 
  const float R =  0.7328F * a_data.x + 0.4296F * a_data.y + -0.1624F * a_data.z;
  const float G = -0.7036F * a_data.x + 1.6975F * a_data.y +  0.0061F * a_data.z;
  const float B =  0.0030F * a_data.x + 0.0136F * a_data.y +  0.9834F * a_data.z;

  a_data.x = R;
  a_data.y = G;
  a_data.z = B;
}

void InverseMatrixCat02(float3& a_data)
{
  // Inverse Mcat02 in CAM02. 
  const float X =  1.096124F * a_data.x + -0.278869F * a_data.y + 0.182745F * a_data.z;
  const float Y =  0.454369F * a_data.x +  0.473533F * a_data.y + 0.072098F * a_data.z;
  const float Z = -0.009628F * a_data.x + -0.005698F * a_data.y + 1.015326F * a_data.z;

  a_data.x = X;
  a_data.y = Y;
  a_data.z = Z;
}

void MatrixHpe(float3& a_data)
{
  // Conversion to cone responses. Matrix HPE in CAM02. (Hunt-Pointer-Estevez)
  const float L =  0.38971F * a_data.x + 0.68898F * a_data.y + -0.07868F * a_data.z; //Kuo modify −0.07869
  const float M = -0.22981F * a_data.x + 1.18340F * a_data.y +  0.04641F * a_data.z; //            0.04642f ?
  const float S = a_data.z;

  a_data.x = L;
  a_data.y = M;
  a_data.z = S;
}

void InverseMatrixHpe(float3& a_data)
{
  // Inverse matrix HPE in CAM02. 
  const float R = 1.910197F * a_data.x + -1.112124F * a_data.y +  0.201908F * a_data.z;
  const float G = 0.370950F * a_data.x +  0.629054F * a_data.y + -0.000008F * a_data.z;
  const float B = a_data.z;

  a_data.x = R;
  a_data.y = G;
  a_data.z = B;
}

void ConvertXyzToLmsVonKries(float3& a_data)
{
  const float L =  0.4002400F * a_data.x + 0.7076000F * a_data.y + -0.0808100F * a_data.z;
  const float M = -0.2263000F * a_data.x + 1.1653200F * a_data.y +  0.0457000F * a_data.z;
  const float S =  0.9182200F * a_data.z;

  a_data.x = L;
  a_data.y = M;
  a_data.z = S;
}

void ConvertLmsToXyzVonKries(float3& a_data)
{
  const float L = a_data.x;
  const float M = a_data.y;
  const float S = a_data.z;

  a_data.x = 1.8599364F * L + -1.1293816F * M +  0.2198974F * S;
  a_data.y = 0.3611914F * L +  0.6388125F * M + -0.0000064F * S;
  a_data.z = 1.0890636F * S;
}

void ConvertXyzToLms(float3& a_data)
{
  float L =  0.4002f * a_data.x + 0.7075f * a_data.y + -0.0807f * a_data.z;
  float M = -0.2280f * a_data.x + 1.1500f * a_data.y +  0.0612f * a_data.z;
  float S =  0.9184f * a_data.z;

  //const float a = 0.43f;
  //
  //if      (L >= 0.0f) L =  pow( L, a);
  //else if (L <  0.0f) L = -pow(-L, a);
  //if      (M >= 0.0f) M =  pow( M, a);
  //else if (M <  0.0f) M = -pow(-M, a);
  //if      (S >= 0.0f) S =  pow( S, a);
  //else if (S <  0.0f) S = -pow(-S, a);

  a_data.x = L;
  a_data.y = M;
  a_data.z = S;
}

void ConvertXyzToLmsPower(float3& a_data, const float a_power)
{
  float L =  0.4002F * a_data.x + 0.7075F * a_data.y + -0.0807F * a_data.z;
  float M = -0.2280F * a_data.x + 1.1500F * a_data.y +  0.0612F * a_data.z;
  float S =  0.9184F * a_data.z;

  const float a = a_power;

  if      (L >= 0.0F) L =  pow( L, a);
  else if (L <  0.0F) L = -pow(-L, a);
  if      (M >= 0.0F) M =  pow( M, a);
  else if (M <  0.0F) M = -pow(-M, a);
  if      (S >= 0.0F) S =  pow( S, a);
  else if (S <  0.0F) S = -pow(-S, a);

  a_data.x = L;
  a_data.y = M;
  a_data.z = S;
}

void ConvertLmsToIpt(float3& a_data)
{
  const float I = 0.4000F * a_data.x +  0.4000F * a_data.y +  0.2000F * a_data.z;
  const float P = 4.4550F * a_data.x + -4.8510F * a_data.y +  0.3960F * a_data.z;
  const float T = 0.8056F * a_data.x +  0.3572F * a_data.y + -1.1628F * a_data.z;

  a_data.x = I;
  a_data.y = P;
  a_data.z = T;
}

void ConvertLmsToXyz(float3& a_data)
{
  float L = a_data.x;
  float M = a_data.y;
  float S = a_data.z;

  a_data.x = 1.8493F * L + -1.1383F * M +  0.2381F * S;
  a_data.y = 0.3660F * L +  0.6444F * M + -0.0100F * S;
  a_data.z = 1.0893F * S;
}

void ConvertLmsToXyzPower(float3& a_data)
{
  float L = a_data.x;
  float M = a_data.y;
  float S = a_data.z;

  const float a = 1.0F / 0.43F;
  
  if      (L >= 0.0F) L =  pow( L, a);
  else if (L <  0.0F) L = -pow(-L, a);
  if      (M >= 0.0F) M =  pow( M, a);
  else if (M <  0.0F) M = -pow(-M, a);
  if      (S >= 0.0F) S =  pow( S, a);
  else if (S <  0.0F) S = -pow(-S, a);

  a_data.x = 1.8493F * L + -1.1383F * M +  0.2381F * S;
  a_data.y = 0.3660F * L +  0.6444F * M + -0.0100F * S;
  a_data.z = 1.0893F * S;
}

void ConvertIptToLms(float3& a_data)
{
  const float L = 0.9999F * a_data.x +  0.0970F * a_data.y +  0.2053F * a_data.z;
  const float M = 0.9999F * a_data.x + -0.1138F * a_data.y +  0.1332F * a_data.z;
  const float S = 0.9999F * a_data.x +  0.0325F * a_data.y + -0.6768F * a_data.z;

  a_data.x = L;
  a_data.y = M;
  a_data.z = S;
}

void ChromAdaptIcam(float3& a_data, const float3 a_dataW, const float3 a_d65, const float a_D)
{
  float3 dataTemp(a_data.x, a_data.y, a_data.z);

  MatrixCat02(dataTemp);

  const float Rw = fmax(a_dataW.x, 1e-6F);
  const float Gw = fmax(a_dataW.y, 1e-6F);
  const float Bw = fmax(a_dataW.z, 1e-6F);

  const float Rc = (a_d65.x * a_D / Rw + 1.0F - a_D) * dataTemp.x;
  const float Gc = (a_d65.y * a_D / Gw + 1.0F - a_D) * dataTemp.y;
  const float Bc = (a_d65.z * a_D / Bw + 1.0F - a_D) * dataTemp.z;

  dataTemp.x = Rc;
  dataTemp.y = Gc;
  dataTemp.z = Bc;

  InverseMatrixCat02(dataTemp);

  a_data.x = dataTemp.x;
  a_data.y = dataTemp.y;
  a_data.z = dataTemp.z;
}

void InverseChromAdaptIcam(float3& a_data, const float D)
{
  static const float3 xyz_d65 = { 0.9505F, 1.0F, 1.0888F };
  static const float Xd65     =  0.8562F * xyz_d65.x +  0.3372F * xyz_d65.y + -0.1934F * xyz_d65.z;  
  static const float Yd65     = -0.8360F * xyz_d65.x +  1.8327F * xyz_d65.y +  0.0033F * xyz_d65.z;
  static const float Zd65     =  0.0357F * xyz_d65.x + -0.0469F * xyz_d65.y +  1.0112F * xyz_d65.z;
    
  const float R =  0.8562F * a_data.x +  0.3372F * a_data.y + -0.1934F * a_data.z; 
  const float G = -0.8360F * a_data.x +  1.8327F * a_data.y +  0.0033F * a_data.z;
  const float B =  0.0357F * a_data.x + -0.0469F * a_data.y +  1.0112F * a_data.z;

  const float Rc = (D * Xd65 / Xd65 + 1.0F - D) * R;
  const float Gc = (D * Yd65 / Yd65 + 1.0F - D) * G;
  const float Bc = (D * Zd65 / Zd65 + 1.0F - D) * B;

  const float Xadapt =  0.9873F * Rc + -0.1768F * Gc + 0.1894F * Bc; 
  const float Yadapt =  0.4504F * Rc +  0.4649F * Gc + 0.0846F * Bc;
  const float Zadapt = -0.0139F * Rc +  0.0278F * Gc + 0.9861F * Bc;

  a_data.x = Xadapt;
  a_data.y = Yadapt;
  a_data.z = Zadapt;
}

float3 Mediana(const float4* inData, const int size, const int a_width, const int a_height, const int x, const int y)
{
  int a = 0;
  std::vector<float> medianaR(size*size);
  std::vector<float> medianaG(size*size);
  std::vector<float> medianaB(size*size);

  const int halfSize = (size - 1) / 2.0f;

  for (int i = -halfSize; i <= halfSize; ++i)
  {
    for (int j = -halfSize; j <= halfSize; ++j)
    {
      int Y = y + i;
      int X = x + j;

      if (Y < 0)         Y = 1;
      else if (Y >= a_height) Y = a_height - 1;
      if (X < 0)         X = 1;
      else if (X >= a_width)  X = a_width - 1;

      medianaR[a] = inData[Y * a_width + X].x;
      medianaG[a] = inData[Y * a_width + X].y;
      medianaB[a] = inData[Y * a_width + X].z;

      a++;
    }
  }
  sort(medianaR.begin(), medianaR.end());
  sort(medianaG.begin(), medianaG.end());
  sort(medianaB.begin(), medianaB.end());

  const int a2 = size * size / 2.0f + 0.5f;
  return { medianaR[a2], medianaG[a2], medianaB[a2] };
}
void Mediana(float inData[], const int radius, const int a_width, const int a_height, const int sizeImage)
{
  std::vector<float> tmpData(sizeImage);

  int sizeMediana = (radius * 2 + 1) * (radius * 2 + 1);
  std::vector<float> mediana(sizeMediana);

  if (radius == 1) sizeMediana = 10;

  for (int y = 0; y < a_height; ++y)
  {
    for (int x = 0; x < a_width; ++x)
    {
      int a = 0;

      for (int i = -radius; i <= radius; ++i)
      {
        for (int j = -radius; j <= radius; ++j)
        {
          int Y = y + i;
          int X = x + j;

          if (Y < 0)         Y = 1;
          else if (Y >= a_height) Y = a_height - 1;
          if (X < 0)         X = 1;
          else if (X >= a_width)  X = a_width - 1;

          mediana[a] = inData[Y * a_width + X];
          a++;
        }
      }
      sort(mediana.begin(), mediana.end());

      tmpData[y * a_width + x] = mediana[sizeMediana / 2];
    }
  }

#pragma omp parallel for
  for (int i = 0; i < sizeImage; ++i)
    inData[i] = tmpData[i];
}

void Blend(float& inData1, const float& inData2, const float amount) // 0 - data1, 1 - data2
{
  inData1 = inData1 + (inData2 - inData1) * amount;
}

template <typename T>
void Normalize(T& data, const float inMin, const float inMax, const float outMin, const float outMax)
{
  data = (data - inMin) * (outMax - outMin) / (inMax - inMin) + outMin;
}
void OffsetCenter(float& data, const float inMin, const float inMax, const float coef)
{
  Normalize(data, inMin, inMax, 0, 1);
  data = pow(data, coef);
  Normalize(data, 0, 1, inMin, inMax);
}
void ContrastField(float& data, const float amount)
{
  if (data > 0.0f && data < 1.0f)
  {
    const float a = data * 1.5707f;
    const float outData = sin(a) * sin(a);
    Blend(data, outData, amount);
  }
}
void Compress(float& data, const float maxRgb)
{
  data = (data * (1.0f + data / (maxRgb * maxRgb))) / (1.0f + data);
}
void CalculateHistogram(const float& data, float histogram[], const int histogramBin)
{
  if (data <= 1.0f)
  {
    // Current bin.      
    const int currentHistogramBin = (data * float(histogramBin - 1)) + 0.5f;

    // Increase count currents bin.
    if (currentHistogramBin >= 0 && currentHistogramBin < histogramBin)
      histogram[currentHistogramBin] += 1.0f;
  }
}
void UniformHistogram(float histogram[], const int histogramBin)
{
  // Recalculate histogramm
  for (int i = 1; i < histogramBin; i++)
  {
    histogram[i] = histogram[i - 1] + histogram[i]; // cummulative uniform
  }
}
void NormalizeHistogram(float histogram[], const int histogramBin)
{
  float histMin = 999999.9f;
  float histMax = 0.0f;

  // Calculate min and max value
//#pragma omp parallel for
  for (int i = 0; i < histogramBin; ++i)
  {
    if (histogram[i] < histMin && histogram[i] != 0)
      histMin = histogram[i];

    if (histogram[i] > histMax)
      histMax = histogram[i];
  }

  // Normalize to 0 - 1
#pragma omp parallel for
  for (int i = 0; i < histogramBin; ++i)
  {
    Normalize(histogram[i], histMin, histMax, 0.0f, 1.0f);
  }
}
void UniformContrastRgb(float* data, const int sizeImage, const int histogramBin)
{
  float histMin = 9999999.9f;
  float histMax = 0.0f;

  const float iterrHistogramBin = 1.0f / (float)(sizeImage * 3);

  std::vector<float> histogram(histogramBin);

  // --- Calculate the histogram ---
  //#pragma omp parallel for
  for (int i = 0; i < sizeImage * 3; i++)
  {
    if (data[i] <= 1.0f)
    {
      const int currentHistogramBin = (int)round(data[i] * float(histogramBin - 1));

      if (currentHistogramBin >= 0 && currentHistogramBin < histogramBin)
        histogram[currentHistogramBin] += iterrHistogramBin;
    }
  }

  // Recalculate the histogram (redistribution)
  for (int i = 1; i < histogramBin; i++)
    histogram[i] = histogram[i - 1] + histogram[i];

  //#pragma omp parallel for
  for (int i = 0; i < histogramBin; i++)
  {
    if (histogram[i] < histMin && histogram[i] != 0)
      histMin = histogram[i];

    if (histogram[i] > histMax)
      histMax = histogram[i];
  }

#pragma omp parallel for
  for (int i = 0; i < histogramBin; i++)
    Normalize(histogram[i], histMin, histMax, 0.0f, 1.0f);

  // Assign a brightness histogram.
#pragma omp parallel for
  for (int i = 0; i < sizeImage * 3; i++)
  {
    if (data[i] <= 1.0f)
    {
      const int currentHistogramBin = (int)round(data[i] * float(histogramBin - 1));

      if (currentHistogramBin >= 0 && currentHistogramBin < histogramBin)
        data[i] = pow(histogram[currentHistogramBin], 2.2f);
      else
        data[i] = 0.0f;
    }
  }
}
void AssignHistogramToData(float& data, const float* histogram, const int histogramBin)
{
  if (data <= 1.0f)
  {
    // Assign the histogram to the data->
    const int currentHistogramBin = (int)round(data * float(histogramBin - 1));

    if (currentHistogramBin >= 0 && currentHistogramBin < histogramBin)
      data = pow(histogram[currentHistogramBin], 2.2f);
    else
      data = 0.0f;
  }
}
void UniformContrast(float data[], const int sizeImage, const int histogramBin)
{
#define MAX3(x,y,z)  ((y) >= (z) ? ((x) >= (y) ? (x) : (y)) : ((x) >= (z) ? (x) : (z)))

  float histMin = 10000000.0f;
  float histMax = 0.0f;

  const float iterrHistogramBin = 1.0f / (float)(sizeImage);

  std::vector<float> histogram(histogramBin);

  // --- Calculate the histogram ---
  //#pragma omp parallel for
  for (int i = 0; i < sizeImage; i++)
  {
    data[i] = data[i] / (1.0f + data[i]);

    const int currentHistogramBin = (int)round(data[i] * float(histogramBin - 1));

    if (currentHistogramBin >= 0 && currentHistogramBin < histogram.size())
      histogram[currentHistogramBin] += iterrHistogramBin;
  }

  // Recalculate the histogram (redistribution)
  for (int i = 1; i < histogramBin; i++)
    histogram[i] = histogram[i - 1] + histogram[i];


#pragma omp parallel for
  for (int i = 0; i < histogramBin; i++)
  {
    if (histogram[i] < histMin && histogram[i] != 0)
      histMin = histogram[i];

    if (histogram[i] > histMax)
      histMax = histogram[i];
  }

#pragma omp parallel for
  for (int i = 0; i < histogramBin; i++)
    Normalize(histogram[i], histMin, histMax, 0.0f, 1.0f);



  // Assign a brightness histogram.
#pragma omp parallel for
  for (int i = 0; i < sizeImage; i++)
  {
    const int currentHistogramBin = (int)round(data[i] * float(histogramBin - 1));

    if (currentHistogramBin >= 0 && currentHistogramBin < histogram.size())
      data[i] = pow(histogram[currentHistogramBin], 2.2f);
    else
      data[i] = 0.0f;
  }
}

void ComputeWhitePoint(const float3 summRgb, float3* whitePoint, const int sizeImage)
{
  whitePoint->x = summRgb.x / (float)sizeImage;
  whitePoint->y = summRgb.y / (float)sizeImage;
  whitePoint->z = summRgb.z / (float)sizeImage;
}

void GetWhitePointForWhiteBalance(const bool a_autoWhiteBalance, const float3 a_summRgb, float3& a_whitePointColor, const int a_sizeImage, float3& a_d65)
{
  if (a_autoWhiteBalance)
    ComputeWhitePoint(a_summRgb, &a_whitePointColor, a_sizeImage);

  ConvertSrgbToXyz(a_whitePointColor);
  MatrixCat02(a_whitePointColor);
  MatrixCat02(a_d65);

  const float lum = Luminance(a_whitePointColor);

  if (lum > 0.0F)
  {
    a_whitePointColor.x /= lum;
    a_whitePointColor.y /= lum;
    a_whitePointColor.z /= lum;
  }
}

void WhiteBalance(float3& a_data, const float3 a_whitePointColor, const float3 a_d65, const float a_whiteBalance)
{
  ConvertSrgbToXyz(a_data);
  ChromAdaptIcam  (a_data, a_whitePointColor, a_d65, a_whiteBalance);
  ConvertXyzToSrgb(a_data);
}


void Saturation(float3& a_data, const float a_saturation)
{
  const float lum = Luminance(a_data);
  a_data          = lum + (a_data - lum) * a_saturation;
  ClampMinusToZero(a_data);
}


void VibranceIPT(float3& a_dataIPT, const float3 a_dataRGB, const float a_vibrance) // saturation of unsaturated colors. 
{  
  float3 normColor     = normalize(a_dataRGB);
  ConvertSrgbToXyz    (normColor);
  ConvertXyzToLmsPower(normColor, 0.43F);
  ConvertLmsToIpt     (normColor);

  const float  satVal  = sqrt(normColor.y * normColor.y + normColor.z * normColor.z) / 0.7852F; // 0.7852F - max sat. in normalize [0-1] IPT.
  const float  satMult = satVal + (1.0F - satVal) * a_vibrance;

  a_dataIPT.y         *= satMult;
  a_dataIPT.z         *= satMult;
}


void CompressIPT(float3& a_dataIPT, const float a_compress)
{
  // Global compress.
  float knee           = 10.0F;
  Blend(knee, 2.0F, pow(a_compress, 0.175F)); // lower = softer
  const float antiKnee = 1.0F / knee;

  const float compLum  = a_dataIPT.x / pow((1.0F + pow(a_dataIPT.x, knee)), antiKnee);
  const float multSat  = pow(compLum / fmax(a_dataIPT.x, 1e-6F), 1.5F);

  a_dataIPT.x          = compLum;
  a_dataIPT.y         *= multSat;
  a_dataIPT.z         *= multSat;
}


void ContrastIPT(float3& a_dataIPT, const float a_contrast)
{
  float lumContr        = a_dataIPT.x;
  ContrastField(lumContr, a_contrast - 1.0F);
  const float multColor = sqrt(lumContr / fmax(a_dataIPT.x, 1e-6F));

  a_dataIPT.x           = lumContr;
  //a_data.y           *= multColor;
  a_dataIPT.z          *= multColor; // compress only blue-yellow axis for fix purple in blue spectrum.
}




void Bilinear(const float inData, float outData[], float newPosX, float newPosY, const int m_width, const int m_height, const int sizeImage)
{
  // |------|------|
  // | dxy3 | dxy4 |
  // |------|------|
  // | dxy1 | dxy2 |
  // |------|------|

  const int floorY = floor(newPosY);
  const int floorX = floor(newPosX);

  const int     dxy1 = floorY * m_width + floorX;
  int           dxy2 = dxy1 + 1;
  if (dxy1 < 0) dxy2 = dxy1 - 1;
  const int     dxy3 = (floorY + 1) * m_width + floorX;
  int           dxy4 = dxy3 + 1;
  if (dxy3 < 0) dxy4 = dxy3 - 1;

  float dx = newPosX - (int)newPosX;
  float dy = newPosY - (int)newPosY;

  if (dx < 0.0f) dx = 1.0f - abs(dx);
  if (dy < 0.0f) dy = 1.0f - abs(dy);

  float multBrightness1 = (1.0f - dx) * (1.0f - dy);
  float multBrightness2 = dx * (1.0f - dy);
  float multBrightness3 = dy * (1.0f - dx);
  float multBrightness4 = dx * dy;

  if (floorY >= 0.0f && floorX >= 0.0f && floorY < m_height - 1 && floorX < m_width - 1)
  {
    outData[dxy1] += inData * multBrightness1;
    outData[dxy2] += inData * multBrightness2;
    outData[dxy3] += inData * multBrightness3;
    outData[dxy4] += inData * multBrightness4;
  }
}
void Bilinear3(const float3 inData, float3 outData[], float newPosX, float newPosY, const int m_width, const int m_height, const int sizeImage)
{
  // |------|------|
  // | dxy3 | dxy4 |
  // |------|------|
  // | dxy1 | dxy2 |
  // |------|------|

  const int floorY = floor(newPosY);
  const int floorX = floor(newPosX);

  const int     dxy1 = floorY * m_width + floorX;
  int           dxy2 = dxy1 + 1;
  if (dxy1 < 0) dxy2 = dxy1 - 1;
  const int     dxy3 = (floorY + 1) * m_width + floorX;
  int           dxy4 = dxy3 + 1;
  if (dxy3 < 0) dxy4 = dxy3 - 1;

  float dx = newPosX - (int)newPosX;
  float dy = newPosY - (int)newPosY;

  if (dx < 0.0f) dx = 1.0f - abs(dx);
  if (dy < 0.0f) dy = 1.0f - abs(dy);

  const float multBrightness1 = (1.0f - dx) * (1.0f - dy);
  const float multBrightness2 = dx * (1.0f - dy);
  const float multBrightness3 = dy * (1.0f - dx);
  const float multBrightness4 = dx * dy;

  if (floorY >= 0.0f && floorX >= 0.0f && floorY < m_height - 1 && floorX < m_width - 1)
  {
    outData[dxy1].x += inData.x * multBrightness1;
    outData[dxy2].x += inData.x * multBrightness2;
    outData[dxy3].x += inData.x * multBrightness3;
    outData[dxy4].x += inData.x * multBrightness4;

    outData[dxy1].y += inData.y * multBrightness1;
    outData[dxy2].y += inData.y * multBrightness2;
    outData[dxy3].y += inData.y * multBrightness3;
    outData[dxy4].y += inData.y * multBrightness4;

    outData[dxy1].z += inData.z * multBrightness1;
    outData[dxy2].z += inData.z * multBrightness2;
    outData[dxy3].z += inData.z * multBrightness3;
    outData[dxy4].z += inData.z * multBrightness4;
  }
}


void Vignette(float4* data4, const float a_vignette, const float centerImageX, const float centerImageY, const float radiusImage, 
  const int x, const int y, const int i)
{
  const float distanceFromCenter = Distance((float)x, (float)y, centerImageX, centerImageY) / radiusImage;

  float vignetteIntensity = sqrt(1.0F - distanceFromCenter * a_vignette);
  if (vignetteIntensity < 0.0f) vignetteIntensity = 0;

  data4[i].x *= vignetteIntensity;
  data4[i].y *= vignetteIntensity;
  data4[i].z *= vignetteIntensity;
}

void SummValueOnField(const float4 data[], float3* summRgb, const int i)
{
#pragma omp atomic
  summRgb->x += data[i].x;
#pragma omp atomic
  summRgb->y += data[i].y;
#pragma omp atomic
  summRgb->z += data[i].z;
}
void ChrommAberr(const float4 inData[], float outData1[], float outData2[], const int m_width, const int m_height, const int sizeImage, const float a_chromAberr, const int x, const int y, const int i)
{
  // Генерируем карту велосити.
  // (-1, -1) - влево вниз, (0, 0) - нет движения. (1, 1) - вправо вверх.  Красный > зелёный > синий.
  const float stepVelocityX = (float)x / m_width;
  const float stepVelocityY = (float)y / m_height;
  const float velocityX = stepVelocityX * 2.0f - 1.0f; // -1.0f to 1.0f
  const float velocityY = stepVelocityY * 2.0f - 1.0f; // -1.0f to 1.0f

  const float newPosXr = x + velocityX * a_chromAberr;
  const float newPosYr = y + velocityY * a_chromAberr;
  const float newPosXg = x + velocityX * 0.5f * a_chromAberr;
  const float newPosYg = y + velocityY * 0.5f * a_chromAberr;

  // Anti-aliasing
  Bilinear(inData[i].x, outData1, newPosXr, newPosYr, m_width, m_height, sizeImage);
  Bilinear(inData[i].y, outData2, newPosXg, newPosYg, m_width, m_height, sizeImage);
}


void MinMaxRgb(float4* data, float& minRgb, float& maxRgb, const float thresholdFloor, const int i)
{
  omp_set_lock;
  if (data[i].x < minRgb && data[i].x > 0.0f) minRgb = data[i].x;
  if (data[i].x > maxRgb)                     maxRgb = data[i].x;

  if (data[i].y < minRgb && data[i].y > 0.0f) minRgb = data[i].y;
  if (data[i].y > maxRgb)                     maxRgb = data[i].y;

  if (data[i].z < minRgb && data[i].z > 0.0f) minRgb = data[i].z;
  if (data[i].z > maxRgb)                     maxRgb = data[i].z;

  if (minRgb < thresholdFloor) minRgb = thresholdFloor;

  if (data[i].x < thresholdFloor || std::isnan(data[i].x)) data[i].x = thresholdFloor;
  if (data[i].y < thresholdFloor || std::isnan(data[i].y)) data[i].y = thresholdFloor;
  if (data[i].z < thresholdFloor || std::isnan(data[i].z)) data[i].z = thresholdFloor;
}

void MinMaxArray(float& data, float& min, float& max, const float thresholdFloor, const float thresholdTop)
{
  if (data < min && data > 0.0f)  min = data;
  else if (data > max)                 max = data;

  if (min < thresholdFloor) min = thresholdFloor;
  if (max > thresholdTop)   max = thresholdTop;

  if (data < thresholdFloor) data = thresholdFloor;
  else if (data > thresholdTop)   data = thresholdTop;
}
void DrawColumn(float4 data[], const int width, const int height, const int offsetX, const int ImgWidth, const float3 color, const bool transparent)
{
  const int endA = height * ImgWidth + width + offsetX;
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
    {
      const int a = y * ImgWidth + x + offsetX;

      if (transparent)
      {
        data[a].x = data[a].x / 2.0f + color.x * (pow((float)a / endA, 4) + 0.2f);
        data[a].y = data[a].y / 2.0f + color.y * (pow((float)a / endA, 4) + 0.2f);
        data[a].z = data[a].z / 2.0f + color.z * (pow((float)a / endA, 4) + 0.2f);
      }
      else
      {
        data[a].x = color.x;
        data[a].y = color.y;
        data[a].z = color.z;
      }
    }
  }
}

void ViewHistorgam(float4 data[], float3* histogram, const int m_width, const int m_height, const int histogramBin)
{
  const int countCol = m_width;
  const int step = histogramBin / countCol + 0.5f;
  int width = m_width / countCol + 0.5f;
  if (width < 1) width = 1;

  float height;
  float gamma = 0.455f;

  float3 color;
  color.x = 0;
  color.y = 0;
  color.z = 0;

  // Draw black histogram
  for (int i = 0; i < countCol; ++i)
  {
    const int j = i * step;

    histogram[j].x = pow(histogram[j].x, gamma);
    histogram[j].y = pow(histogram[j].y, gamma);
    histogram[j].z = pow(histogram[j].z, gamma);

    height = histogram[j].x * m_height;
    DrawColumn(data, width, (int)height, i * width, m_width, color, 0);

    height = histogram[j].y * m_height;
    DrawColumn(data, width, (int)height, i * width, m_width, color, 0);

    height = histogram[j].z * m_height;
    DrawColumn(data, width, (int)height, i * width, m_width, color, 0);
  }

  // Draw color histogram
  for (int i = 0; i < countCol; ++i)
  {
    const int j = i * step;

    color.x = 0.0f;
    color.y = 0.0f;
    color.z = 0.5f;
    height = histogram[j].z * m_height;
    DrawColumn(data, width, (int)height, i * width, m_width, color, 1);

    color.x = 0.5f;
    color.y = 0.0;
    color.z = 0.0;
    height = histogram[j].x * m_height;
    DrawColumn(data, width, (int)height, i * width, m_width, color, 1);

    color.x = 0.0f;
    color.y = 0.5f;
    color.z = 0.0f;
    height = histogram[j].y * m_height;
    DrawColumn(data, width, (int)height, i * width, m_width, color, 1);

  }

}

void MinMaxHistBin(const float* histogram, float& minHistBin, float& maxHistBin, const int sizeImage, const int histogramBin)
{
  int i = 0;
  const float floor = (float)sizeImage * 0.00001f; // 0.001% from image resolution

  while (i < histogramBin && histogram[i] <= floor)
  {
    minHistBin = (float)i;
    i++;
  }

  i = histogramBin - 1;

  while (i >= 0 && histogram[i] <= floor)
  {
    maxHistBin = (float)i;
    i--;
  }

  minHistBin /= histogramBin;
  maxHistBin /= histogramBin;
}


int FloatToInt(const float inData)
{
  return int(inData + 0.5f);
}

void Resize(float data[], const float sourceWidth, const float sourceHeight, const float sourceSizeImage, const float resize)
{
  // Resize works correctly only if the resize is a multiple of 2.
  const int newSizeImage = sourceSizeImage * (resize * resize);
  const int newHeight = sourceHeight * resize;
  const int newWidth = sourceWidth * resize;

  float* resizeArray = new float[newSizeImage];

  // Downsample
  if (resize < 1.0f)
  {
#pragma omp parallel for
    for (int y = 0; y < newHeight; ++y)
    {
      for (int x = 0; x < newWidth; ++x)
      {
        const int i = y * newWidth + x;
        const int j = (int)(y / resize) * sourceWidth + x / resize;
        if (j >= sourceSizeImage) continue;
        resizeArray[i] = data[j];
      }
    }

#pragma omp parallel for
    for (int y = 0; y < newHeight; ++y)
    {
      for (int x = 0; x < newWidth; ++x)
      {
        const int i = y * sourceWidth + x;
        const int j = y * newWidth + x;

        data[i] = resizeArray[j];
      }
    }
  }
  else // Upsample
  {
#pragma omp parallel for
    for (int y = 0; y < newHeight; ++y)
    {
      for (int x = 0; x < newWidth; ++x)
      {
        const int i = y * newWidth + x;
        const int j = (int)(y / resize) * newWidth + x / resize;
        resizeArray[i] = data[j];
      }
    }

#pragma omp parallel for
    for (int i = 0; i < newSizeImage; ++i)
      data[i] = resizeArray[i];
  }

  delete[] resizeArray;
}

void Sharp(float4 image4out[], const float lumForSharp[], const float a_sharpness, const int a_width, const int a_height)
{
#pragma omp parallel for
  for (int y = 0; y < a_height; ++y)
  {
    for (int x = 0; x < a_width; ++x)
    {
      const int i = y * a_width + x;

      if (image4out[i].x < 0.0f) image4out[i].x = 0.0f;
      if (image4out[i].y < 0.0f) image4out[i].y = 0.0f;
      if (image4out[i].z < 0.0f) image4out[i].z = 0.0f;

      float mean = 0.0f;

      for (int i = -1; i <= 1; ++i)
      {
        for (int j = -1; j <= 1; ++j)
        {
          int Y = y + i;
          int X = x + j;

          if (Y < 0)         Y = 1;
          else if (Y >= a_height) Y = a_height - 1;
          if (X < 0)         X = 1;
          else if (X >= a_width)  X = a_width - 1;

          mean += lumForSharp[Y * a_width + X];
        }
      }
      mean /= 9.0f;

      float dispers = 0.0f;
      for (int i = -1; i <= 1; i++)
      {
        for (int j = -1; j <= 1; j++)
        {
          int Y = y + i;
          int X = x + j;

          if (Y < 0)         Y = 1;
          else if (Y >= a_height) Y = a_height - 1;
          if (X < 0)         X = 1;
          else if (X >= a_width)  X = a_width - 1;

          float a = lumForSharp[Y * a_width + X] - mean;
          dispers += abs(a);
        }
      }

      dispers /= (1.0f + dispers);
      const float lumCompr = lumForSharp[i] / (1.0f + lumForSharp[i]);
      const float meanCompr = mean / (1.0f + mean);
      float hiPass = (lumForSharp[i] - mean) * 5 + 1.0f;
      hiPass = pow(hiPass, 2);

      float sharp = 1.0f;
      Blend(sharp, hiPass, a_sharpness * (1.0f - dispers));

      image4out[i].x *= sharp;
      image4out[i].y *= sharp;
      image4out[i].z *= sharp;
    }
  }
}

std::vector<float> createGaussKernelWeights1D_HDRImage2(const int blurRadius)
{
  std::vector<float> gKernel(blurRadius);

  // sum is for normalization
  float sum = 0.0f;

  for (int x = 0; x < blurRadius; ++x)
  {
    gKernel[x] = exp(-x * x / ((blurRadius * blurRadius) / 6.28318530f));
    sum += gKernel[x];
  }

  // normalize the Kernel
  for (int i = 0; i < blurRadius; ++i)
    gKernel[i] /= (2.0f * sum);

  return gKernel;
}

void Blur(float data[], int blurRadius, const int m_width, const int m_height, const int sizeImage)
{
  const float diagonalImage = Distance(0, 0, m_width, m_height);
  float radiusImage = diagonalImage / 2.0f;

  // Resize works correctly only if the resize is a multiple of 2.
  float resize = 1.0f;
  if (blurRadius >= (radiusImage / 4.0f)) resize = 0.5f;
  else if (blurRadius >= (radiusImage / 2.0f)) resize = 0.25f;
  blurRadius *= resize;

  const int resizeWidth     = m_width * resize;
  const int resizeHeight    = m_height * resize;
  const int resizeSizeImage = sizeImage * (resize * resize);

  float* blurPass = new float[resizeSizeImage];

  std::vector<float> weights = createGaussKernelWeights1D_HDRImage2(blurRadius + 1);

  // Downsize for speed up blur
  if (resize < 1.0f) Resize(data, m_width, m_height, sizeImage, resize);

  //------------ Blur ------------

  float summ;
  float multRadius = 1.0f;
  if (blurRadius >= resizeWidth)
    multRadius = (float)(resizeWidth - 1) / blurRadius;

  // All row.  

// cannot be parallelize.
  for (int y = 0; y < resizeHeight; ++y)
  {
    for (int x = 0; x < resizeWidth; ++x)
    {
      summ = 0.0f;

      for (int row = -blurRadius; row <= blurRadius; ++row)
      {
        int currX = x + row;
        if (currX < 0)            currX = 0 + abs(row) * multRadius;
        else if (currX >= resizeWidth) currX = resizeWidth - abs(row) * multRadius;

        summ += data[y * m_width + currX] * weights[abs(row)];
      }
      blurPass[y * resizeWidth + x] = summ;
    }
  }

  // All col.

  multRadius = 1.0f;
  if (blurRadius >= resizeHeight)
    multRadius = (float)(resizeHeight - 1) / blurRadius;

#pragma omp parallel for
  for (int x = 0; x < resizeWidth; ++x)
  {
    for (int y = 0; y < resizeHeight; ++y)
    {
      summ = 0.0f;

      for (int col = -blurRadius; col <= blurRadius; ++col)
      {
        int currY = y + col;
        if (currY < 0)             currY = 0 + abs(col) * multRadius;
        else if (currY >= resizeHeight) currY = resizeHeight - abs(col) * multRadius;

        summ += blurPass[currY * resizeWidth + x] * weights[abs(col)];
      }
      data[y * m_width + x] = summ;
    }
  }

  // Upsize
  resize = 1.0f / resize;
  if (resize > 1.0f) Resize(data, resizeWidth, resizeHeight, resizeSizeImage, resize);

  delete[] blurPass;
}


float ConvertAngleToRad(int angle)
{
  return angle * 0.01745329252f;
}
void DiffractionStars(float4* inData, float3 diffrStars[], const float a_sizeStar, const int a_numRay, const int a_RotateRay, const float a_randomAngle,
                      const float a_sprayRay, const int m_width, const int m_height, const int sizeImage, const float radiusImage, const int x, const int y, const int i)
{
  const float twoPI = 6.283185f;
  const float angle = twoPI / a_numRay;
  const float waveLenghtR = 0.000700f;
  const float waveLenghtY = 0.000575f;
  const float waveLenghtG = 0.000530f;
  const float waveLenghtC = 0.000485f;
  const float waveLenghtB = 0.000450f;
  const float waveLenghtV = 0.000400f;

  float lum = Luminance(float3(inData[i].x, inData[i].y, inData[i].z));
  lum /= (1.0f + lum);

  // Ray
  for (int numRay = 0; numRay < a_numRay; ++numRay)
  {
    float jitter = (rand() % 100) / 100.0f + 0.5f;
    float nextAngle = angle * numRay - ConvertAngleToRad(a_RotateRay);
    if ((rand() % 1000) / 1000.0f <= a_randomAngle)
    {
      nextAngle = (rand() % 628) / 100.0f;
    }
    const float sprayAngle = ((float)(rand() % (int)(angle * 100)) / 100.0f) * a_sprayRay;
    nextAngle += (sprayAngle - angle / 2.0f);

    const float sizeStar = a_sizeStar * radiusImage / 100.0f;
    const float sizeStarFromLum = sizeStar * lum * jitter / 5.0f;

    // Aperture big   - dist ray small, wave small
    // Aperture small - dist ray big, wave big
    const float diafragma = 1.0f / (sizeStarFromLum * 5.0f * jitter);

    for (float sample = 1.0f; sample < sizeStarFromLum; sample += 0.1f)
    {
      const float angleFi = sample;
      const float temp = diafragma * angleFi;
      const float multDist = 1.0f - (float)sample / sizeStarFromLum;
      const float u = temp / waveLenghtG + 0.000001f;

      const float newX = x + cos(nextAngle) * sample;
      const float newY = y + sin(nextAngle) * sample;

      if (newX > 0 && newX < m_width && newY > 0 && newY < m_height)
      {
        const float I = pow((sinf(u) / u), 2) * multDist;

        const float3 color ={ inData[i].x * I, inData[i].y * I ,inData[i].z * I };

        Bilinear3(color, diffrStars, newX, newY, m_width, m_height, sizeImage);
      }
    }
  }
}

void UniformContrastRGBFilter(float4 image4out[], const int sizeImage, const int histogramBin, const float amount)
{
  std::vector<float> rgbArray(sizeImage * 3);

  // Convert 3 field RGB to linear array.
#pragma omp parallel for
  for (int i = 0; i < sizeImage; ++i)
  {
    rgbArray[i] = image4out[i].x;
    rgbArray[sizeImage + i] = image4out[i].y;
    rgbArray[sizeImage * 2 + i] = image4out[i].z;
  }

  UniformContrastRgb(&rgbArray[0], sizeImage, histogramBin);

  // Return to main array
#pragma omp parallel for
  for (int i = 0; i < sizeImage; ++i)
  {
    const float meanRGBsource = (image4out[i].x + image4out[i].y + image4out[i].z) / 3.0F;
    const float meanRGB_UC    = (rgbArray[i] + rgbArray[i + sizeImage] + rgbArray[i + sizeImage * 2]) / 3.0F;

    float diff                = meanRGB_UC / fmax(meanRGBsource, 1e-6F);
    float3 rgbDiff;
    rgbDiff.x                 = image4out[i].x;
    rgbDiff.y                 = image4out[i].y;
    rgbDiff.z                 = image4out[i].z;

    Blend(rgbDiff.x, image4out[i].x * diff, amount);
    Blend(rgbDiff.y, image4out[i].y * diff, amount);
    Blend(rgbDiff.z, image4out[i].z * diff, amount);

    Blend(image4out[i].x, rgbArray[i], amount);
    Blend(image4out[i].y, rgbArray[i + sizeImage], amount);
    Blend(image4out[i].z, rgbArray[i + sizeImage * 2], amount);

    Blend(image4out[i].x, rgbDiff.x, 0.5F);
    Blend(image4out[i].y, rgbDiff.y, 0.5F);
    Blend(image4out[i].z, rgbDiff.z, 0.5F);
  }
}

void NormalizeFilter(float4 image4out[], float3 histogram[], const int sizeImage, const int histogramBin, float3& minHistBin,
                     float3& maxHistBin, float& minAllHistBin, float& maxAllHistBin, const float amount)
{
#pragma omp parallel for
  for (int i = 0; i < sizeImage; ++i)
  {
    CalculateHistogram(image4out[i].x, &histogram[0].x, histogramBin);
    CalculateHistogram(image4out[i].y, &histogram[0].y, histogramBin);
    CalculateHistogram(image4out[i].z, &histogram[0].z, histogramBin);
  }

  // Calculate min/max histogram for a_normalize
  MinMaxHistBin(&histogram[0].x, minHistBin.x, maxHistBin.x, sizeImage, histogramBin);
  MinMaxHistBin(&histogram[0].y, minHistBin.y, maxHistBin.y, sizeImage, histogramBin);
  MinMaxHistBin(&histogram[0].z, minHistBin.z, maxHistBin.z, sizeImage, histogramBin);

  minAllHistBin = Min3(minHistBin.x, minHistBin.y, minHistBin.z);
  maxAllHistBin = Max3(maxHistBin.x, maxHistBin.y, maxHistBin.z);

  // Normalize 
#pragma omp parallel for
  for (int i = 0; i < sizeImage; ++i)
  {
    float3 dataNorm;
    dataNorm.x = image4out[i].x;
    dataNorm.y = image4out[i].y;
    dataNorm.z = image4out[i].z;

    Normalize(dataNorm.x, minAllHistBin, maxAllHistBin, 0.0f, 1.0f);
    Normalize(dataNorm.y, minAllHistBin, maxAllHistBin, 0.0f, 1.0f);
    Normalize(dataNorm.z, minAllHistBin, maxAllHistBin, 0.0f, 1.0f);

    Blend(image4out[i].x, dataNorm.x, amount);
    Blend(image4out[i].y, dataNorm.y, amount);
    Blend(image4out[i].z, dataNorm.z, amount);
  }
}

void ViewHistorgamFilter(float4 image4out[], float3 histogram[], const int sizeImage, const int histogramBin, const int a_width,
                         const int a_height)
{
#pragma omp parallel for
  for (int i = 0; i < sizeImage; ++i)
  {
    CalculateHistogram(image4out[i].x, &histogram[0].x, histogramBin);
    CalculateHistogram(image4out[i].y, &histogram[0].y, histogramBin);
    CalculateHistogram(image4out[i].z, &histogram[0].z, histogramBin);
  }

  NormalizeHistogram(&histogram[0].x, histogramBin);
  NormalizeHistogram(&histogram[0].y, histogramBin);
  NormalizeHistogram(&histogram[0].z, histogramBin);

  ViewHistorgam(image4out, &histogram[0], a_width, a_height, histogramBin);
}

