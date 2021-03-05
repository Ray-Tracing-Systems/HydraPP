#include "Functions.h"

//--------------------------------------------------------------------

float Min(const float value1, const float value2, const float value3)
{
  return fmin(fmin(value1, value2), value3);
}

float Min(const float3 value1, const float value2)
{
  return fmin(Min(value1.x, value1.y, value1.z), value2);
}

float Max(const float value1, const float value2, const float value3)
{
  return fmax(fmax(value1, value2), value3);
}

float Distance(const float x1, const float y1, const float x2, const float y2)
{
  float differenceX = x1 - x2;
  float differenceY = y1 - y2;

  float distance = sqrt(differenceX * differenceX + differenceY * differenceY);
  return distance;
}

float Distance(const float a, const float b)
{   
  return sqrt(a * a + b * b);
}

float Distance(const float3& point1, const float point2)
{
  float3 difference(point1.x - point2, point1.y - point2, point1.z - point2);

  float distance = sqrt(difference.x * difference.x + difference.y * difference.y + difference.z * difference.z);
  return distance;
}

void ClampMinusToZero(float3& a_data)
{
  a_data.x = fmax(a_data.x, 0.0F);
  a_data.y = fmax(a_data.y, 0.0F);
  a_data.z = fmax(a_data.z, 0.0F);
}

float Clamp(const float u, const float a, const float b) { const float r = fmax(a, u); return fmin(r, b); }

[[nodiscard]] float Luminance(const float3& a_data)
{  
  return dot(a_data, float3(0.2126F, 0.7152F, 0.0722F));
}


[[nodiscard]] float Mean(const float3& a_data)
{
  return (a_data.x + a_data.y + a_data.z) / 3.0F;
}

void SpectralToRGB(float &r, float &g, float &b, const float l) // RGB <0,1> <- lambda l <400,700> [nm]
{
  float t;
  r = 0.0;
  g = 0.0;
  b = 0.0;

  if      ((l >= 400.0f) && (l < 410.0f)) { t = (l - 400.0f) / (410.0f - 400.0f); r = +(0.33f * t) - (0.20f * t * t); }
  else if ((l >= 410.0f) && (l < 475.0f)) { t = (l - 410.0f) / (475.0f - 410.0f); r = 0.14f - (0.13f * t * t); }
  else if ((l >= 545.0f) && (l < 595.0f)) { t = (l - 545.0f) / (595.0f - 545.0f); r = +(1.98f * t) - (t * t); }
  else if ((l >= 595.0f) && (l < 650.0f)) { t = (l - 595.0f) / (650.0f - 595.0f); r = 0.98f + (0.06f * t) - (0.40f * t * t); }
  else if ((l >= 650.0f) && (l < 700.0f)) { t = (l - 650.0f) / (700.0f - 650.0f); r = 0.65f - (0.84f * t) + (0.20f * t * t); }
  if      ((l >= 415.0f) && (l < 475.0f)) { t = (l - 415.0f) / (475.0f - 415.0f); g = +(0.80f * t * t); }
  else if ((l >= 475.0f) && (l < 590.0f)) { t = (l - 475.0f) / (590.0f - 475.0f); g = 0.8f + (0.76f * t) - (0.80f * t * t); }
  else if ((l >= 585.0f) && (l < 639.0f)) { t = (l - 585.0f) / (639.0f - 585.0f); g = 0.84f - (0.84f * t); }
  if      ((l >= 400.0f) && (l < 475.0f)) { t = (l - 400.0f) / (475.0f - 400.0f); b = +(2.20f * t) - (1.50f * t * t); }
  else if ((l >= 475.0f) && (l < 560.0f)) { t = (l - 475.0f) / (560.0f - 475.0f); b = 0.7f - (t) + (0.30f * t * t); }
}

void ConvertRgbToHsv(float3& a_data)
{
  float r = a_data.x;
  float g = a_data.y;
  float b = a_data.z;

  r *= 255.0F;
  g *= 255.0F;
  b *= 255.0F;

  float rgb_max = Max(r, g, b);

  float H = 0.0F;
  float V = rgb_max;

  // norm. val 1
  if (V != 0.0F)
  {
    r /= V;
    g /= V;
    b /= V;
  }

  float rgb_min = Min(r, g, b);
  rgb_max       = Max(r, g, b);

  const float S = rgb_max - rgb_min;

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
  float var_R = data->x; //var_R From 0 to 1
  float var_G = data->y; //var_G From 0 to 1
  float var_B = data->z; //var_B From 0 to 1

  if (var_R > 0.04045F) var_R = pow((var_R + 0.055F) / 1.055F, 2.4F);
  else                  var_R = var_R / 12.92F;

  if (var_G > 0.04045F) var_G = pow((var_G + 0.055F) / 1.055F, 2.4F);
  else                  var_G = var_G / 12.92F;

  if (var_B > 0.04045F) var_B = pow((var_B + 0.055F) / 1.055F, 2.4F);
  else                  var_B = var_B / 12.92F;

  var_R *= 100.0F;
  var_G *= 100.0F;
  var_B *= 100.0F;

  //Chromatic adaptation. Observer. = 2°, Illuminant = D65

  data->x = var_R * 0.4124564F + var_G * 0.3575761F + var_B * 0.1804375F;
  data->y = var_R * 0.2126729F + var_G * 0.7151522F + var_B * 0.0721750F;
  data->z = var_R * 0.0193339F + var_G * 0.1191920F + var_B * 0.9503041F;
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


void ConvertSrgbToXyzIcam06(float4* data)
{
  const float R = data->x;
  const float G = data->y;
  const float B = data->z;

  data->x = R * 0.412424f + G * 0.2126560f + B * 0.0193324f;
  data->y = R * 0.357579f + G * 0.7151580f + B * 0.1191900f;
  data->z = R * 0.180464f + G * 0.0721856f + B * 0.9504440f;
}

void ConvertXyz65ToRgb(float3& a_data)
{
  float var_X = a_data.x / 100.0F;  //X From 0 to  95.047      (Observer = 2°, Illuminant = D65)
  float var_Y = a_data.y / 100.0F;  //Y From 0 to 100.000
  float var_Z = a_data.z / 100.0F;  //Z From 0 to 108.883

  //Chromatic adaptation.

  float R = var_X *  3.2404542F + var_Y * -1.5371385F + var_Z * -0.4985314F;
  float G = var_X * -0.9692660F + var_Y *  1.8760108F + var_Z *  0.0415560F;
  float B = var_X *  0.0556434F + var_Y * -0.2040259F + var_Z *  1.0572252F;

  if (R > 0.0031308F) R = 1.055F * (pow(R, 0.4166666666F)) - 0.055F; // 1.0F / 2.4F = 0.4166666666F
  else                R = 12.92F * R;
  if (G > 0.0031308F) G = 1.055F * (pow(G, 0.4166666666F)) - 0.055F;
  else                G = 12.92F * G;
  if (B > 0.0031308F) B = 1.055F * (pow(B, 0.4166666666F)) - 0.055F;
  else                B = 12.92F * B;

  a_data.x = R;
  a_data.y = G;
  a_data.z = B;
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

void ConvertXyzToSrgbIcam06(float3& a_data)
{
  const float X = a_data.x;
  const float Y = a_data.y;
  const float Z = a_data.z;

  a_data.x = X *  3.2407F + Y * -0.9693F + Z *  0.0556F;
  a_data.y = X * -1.5373F + Y *  1.8760F + Z * -0.2040F;
  a_data.z = X * -0.4986F + Y *  0.0416F + Z *  1.0571F;
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


void ConvertXyzToZlab(float4* data, const float3& whitePoint, const float a_whiteBalance)
{
  const float X = data->x;
  const float Y = data->y;
  const float Z = data->z;

  const float Xw = whitePoint.x;
  const float Yw = whitePoint.y;
  const float Zw = whitePoint.z;

  // ----- chromatic adaptation -----

  // convert to sharpened cone responses using Бредфорский расчёт. 

  const float R =  0.8951f * (X / Y) +  0.2664f + -0.1614f * (Z / Y);
  const float G = -0.7502f * (X / Y) +  1.7135f +  0.0367f * (Z / Y);
  const float B =  0.0389f * (X / Y) + -0.0685f +  1.0296f * (Z / Y);

  const float Rw =  0.8951f * (Xw / Yw) +  0.2664f + -0.1614f * (Zw / Yw);
  const float Gw = -0.7502f * (Xw / Yw) +  1.7135f +  0.0367f * (Zw / Yw);
  const float Bw =  0.0389f * (Xw / Yw) + -0.0685f +  1.0296f * (Zw / Yw);

  const float p = pow(Bw, 0.0834f);
  const float D = a_whiteBalance;

  const float Rc = (D * (1.0f / Rw) + 1.0f - D) * R;
  const float Gc = (D * (1.0f / Gw) + 1.0f - D) * G;
  float       Bc = (D * (1.0f / pow(Bw, p)) + 1.0f - D) * pow(abs(B), p);

  if (B < 0.0f) Bc = -abs(B);

  // Back all to XYZ with inverse M matrix

  const float Xc = 0.9871f * Rc * Y + -0.1469f * Gc * Y + 0.1597f * Bc * Y;
  const float Yc = 0.4324f * Rc * Y + 0.5184f * Gc * Y + 0.0491f * Bc * Y;
  const float Zc = -0.0083f * Rc * Y + 0.0402f * Gc * Y + 0.9687f * Bc * Y;

  // ----- Appearance Correlates -----


  const float Lz = 100.0f *  pow(Yc / 100.0f, 0.50025f);
  const float Az = 500.0f * (pow(Xc / 100.0f, 0.345f) - pow(Yc / 100.0f, 0.345f));
  const float Bz = 200.0f * (pow(Yc / 100.0f, 0.345f) - pow(Zc / 100.0f, 0.345f));

  data->x = Lz;
  data->y = Az;
  data->z = Bz;
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
void ConvertXyzToRlab(float4* data, float power)
{
  const float var_X = data->x / 95.0470f;   //Observer= 2°, Illuminant= D65 
  const float var_Y = data->y / 100.000f;   //
  const float var_Z = data->z / 108.883f;   //

  const float p = 1.0f / (2.3f * power + 0.001f);

  data->x = 100.0f *  pow(var_Y, p);                   // L
  data->y = 430.0f * (pow(var_X, p) - pow(var_Y, p));  // a
  data->z = 170.0f * (pow(var_Y, p) - pow(var_Z, p));  // b
}
void ConvertRlabToXyz(float4* data)
{
  const float var_Y = pow(data->x / 100.0f, 2.3f);
  const float var_X = pow(pow(var_Y, 1.0f / 2.3f) + data->y / 430.0f, 2.3f);
  const float var_Z = pow(pow(var_Y, 1.0f / 2.3f) - data->z / 170.0f, 2.3f);

  data->x = var_X * 95.047f;   // Observer= 2°, Illuminant= D65
  data->y = var_Y * 100.0f;    //
  data->z = var_Z * 108.883f;  //
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
void ConvertZlabToXyz(float4* data, const float3& whitePoint)
{
  const float Lz = data->x;
  const float Az = data->y;
  const float Bz = data->z;

  const float Rw = whitePoint.x;
  const float Gw = whitePoint.y;
  const float Bw = whitePoint.z;

  // ----- Appearance Correlates -----

  const float Yc = pow(Lz, 1.999f);
  const float Xc = pow(pow(Yc, 0.345f) + Az / 5.0f, 2.8985f);
  const float Zc = pow(pow(Yc, 0.345f) - Bz / 2.0f, 2.8985f);

  // ----- chromatic adaptation -----

  // M matrix

  const float Rc = 0.8951f * Xc + 0.2664f * Yc + -0.1614f * Zc;
  const float Gc = -0.7502f * Xc + 1.7135f * Yc + 0.0367f * Zc;
  const float Bc = 0.0389f * Xc + -0.0685f * Yc + 1.0296f * Zc;

  const float D = 0.0f;
  const float p = pow(Bw, 0.0834f);

  const float R = Rc / (D * (1.0f / Rw) + 1.0f - D);
  const float G = Gc / (D * (1.0f / Gw) + 1.0f - D);
  const float B = pow(abs(Bc / (D * (1.0f / pow(Bw, p)) + 1.0f - D)), 1 / p);

  // inverse M matrix

  const float X = 0.9871f * R + -0.1469f * G + 0.1597f * B;
  const float Y = 0.4324f * R + 0.5184f * G + 0.0491f * B;
  const float Z = -0.0083f * R + 0.0402f * G + 0.9687f * B;

  data->x = Xc * 100.0f;
  data->y = Yc * 100.0f;
  data->z = Zc * 100.0f;
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


void ExponentIcam(float3& a_data, const float a_dataYblur)
{
  float L = a_data.x;
  float M = a_data.y;
  float S = a_data.z;

  const float La = 0.2F * a_dataYblur;

  const float Fl = 0.2F * pow(1.0F / (5.0F * La + 1.0F), 4.0F)        *    (5.0F * La) +
        0.1F * pow(1.0F - pow(1.0F / (5.0F * La + 1.0F), 4.0F), 2.0F) * pow(5.0F * La, 0.3333333F);
      
  const float a = fmax(0.43F * Fl, 0.3F);

  if      (L >= 0.0F) a_data.x =  pow( L, a);
  else if (L <  0.0F) a_data.x = -pow(-L, a);

  if      (M >= 0.0F) a_data.y =  pow( M, a);
  else if (M <  0.0F) a_data.y = -pow(-M, a);

  if      (S >= 0.0F) a_data.z =  pow( S, a);
  else if (S <  0.0F) a_data.z = -pow(-S, a);
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

void ConvertLmsToXyzPower(float3& a_data, const float a_power)
{
  float L = a_data.x;
  float M = a_data.y;
  float S = a_data.z;  

  if      (L >= 0.0F) L =  pow( L, a_power);
  else if (L <  0.0F) L = -pow(-L, a_power);
  if      (M >= 0.0F) M =  pow( M, a_power);
  else if (M <  0.0F) M = -pow(-M, a_power);
  if      (S >= 0.0F) S =  pow( S, a_power);
  else if (S <  0.0F) S = -pow(-S, a_power);

  a_data.x = 1.8502F * L + -1.1383F * M +  0.2384F * S;
  a_data.y = 0.3668F * L +  0.6438F * M + -0.0106F * S;
  a_data.z = 1.0888F * S;
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

void ChromAdaptIcam06(float3& a_data, const float3& a_dataW, const float3& a_d65, const float a_D)
{
  float3 dataTemp(a_data.x, a_data.y, a_data.z);
  MatrixCat02(dataTemp);

  const float Rw = fmax(a_dataW.x, 1e-6F);
  const float Gw = fmax(a_dataW.y, 1e-6F);
  const float Bw = fmax(a_dataW.z, 1e-6F);

  const float fakeBlurLum = a_data.y;
  //const float La = 0.2F * fakeBlurLum;
  const float F  = a_D; 
  //const float D = /*0.3F **/ F* (1.0F - 1.0F / 3.6F * exp(-(La - 42.0F) / 92.0F));  
  const float D = F * WeightWhiteBalance(a_data.y);
  
  const float Rc = (a_d65.x * D / Rw + (1.0F - D)) * dataTemp.x;
  const float Gc = (a_d65.y * D / Gw + (1.0F - D)) * dataTemp.y;
  const float Bc = (a_d65.z * D / Bw + (1.0F - D)) * dataTemp.z;

  dataTemp.x     = Rc;
  dataTemp.y     = Gc;
  dataTemp.z     = Bc;

  InverseMatrixCat02(dataTemp);

  a_data.x       = dataTemp.x;
  a_data.y       = dataTemp.y;
  a_data.z       = dataTemp.z;
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

void ChromAdaptCam02(float3& a_data, const float3& a_dataW, const float a_D)
{
  float3 dataTemp(a_data.x, a_data.y, a_data.z);

  MatrixCat02(dataTemp);

  // D65 offset: 95.047f / 100.0f = 0.95047f
  // D65 offset: 108.883 / 100.0f = 1.08883f
  const float Rc = (0.95047F * a_D / a_dataW.x + (1.0F - a_D)) * dataTemp.x;
  const float Gc = (           a_D / a_dataW.y + (1.0F - a_D)) * dataTemp.y;
  const float Bc = (1.08883F * a_D / a_dataW.z + (1.0F - a_D)) * dataTemp.z;

  dataTemp.x = Rc;
  dataTemp.y = Gc;
  dataTemp.z = Bc;
  
  InverseMatrixCat02(dataTemp);

  a_data.x = dataTemp.x;
  a_data.y = dataTemp.y;
  a_data.z = dataTemp.z;
}

void InverseChromAdaptCam02(float3& a_data, const float3& a_dataW, const float a_Fl, const float a_D, const float a_c)
{
  float3 dataTemp(a_data.x, a_data.y, a_data.z);

  MatrixCat02(dataTemp);

  const float Rc = dataTemp.x;
  const float Gc = dataTemp.y;
  const float Bc = dataTemp.z;

  const float R = Rc / (100.0F * a_D / a_dataW.x + 1.0F - a_D);
  const float G = Gc / (100.0F * a_D / a_dataW.y + 1.0F - a_D);
  const float B = Bc / (100.0F * a_D / a_dataW.z + 1.0F - a_D);

  dataTemp.x = R;
  dataTemp.y = G;
  dataTemp.z = B;

  InverseMatrixCat02(dataTemp); // Return XYZ

  a_data.x = dataTemp.x;
  a_data.y = dataTemp.y;
  a_data.z = dataTemp.z;
}



void GetWhitePointForWhiteBalance(const bool a_autoWhiteBalance, const float3& a_summRgb, float3& a_whitePointColor, const int a_sizeImage, float3& a_d65)
{
  if (a_autoWhiteBalance)
    ComputeWhitePoint(a_summRgb, a_whitePointColor, a_sizeImage);

  
  if (float lum = Luminance(a_whitePointColor); lum > 0.0F)
  {
    a_whitePointColor.x /= lum;
    a_whitePointColor.y /= lum;
    a_whitePointColor.z /= lum;
  }

  ConvertSrgbToXyz(a_whitePointColor);
  MatrixCat02(a_whitePointColor);
  MatrixCat02(a_d65);
}

void WhiteBalance(float3& a_data, const float3& a_whitePointColor, const float3& a_d65, const float a_whiteBalance)
{
  ConvertSrgbToXyz(a_data);
  ChromAdaptIcam06(a_data, a_whitePointColor, a_d65, a_whiteBalance);
  ConvertXyzToSrgb(a_data);
}


void Saturation(float3& a_data, const float a_saturation)
{
  const float lum = Luminance(a_data);
  a_data          = lum + (a_data - lum) * a_saturation;
  ClampMinusToZero(a_data);
}


void VibranceIPT(float3& a_dataIPT, const float3& a_dataRGB, const float a_vibrance) // saturation of unsaturated colors. 
{    
  float3 normColor     = normalize(a_dataRGB);
  ConvertSrgbToXyz    (normColor);
  ConvertXyzToLmsPower(normColor, 0.43F);
  ConvertLmsToIpt     (normColor);

  const float satVal  = sqrt(normColor.y * normColor.y + normColor.z * normColor.z) / 0.7852F; // 0.7852F - max sat. in normalize [0-1] IPT.
  
  if (satVal < 0.01F) // fix falloff some desaturated images to black.
    return;

  const float satMult = satVal + (1.0F - satVal) * a_vibrance;

  a_dataIPT.y         *= satMult;
  a_dataIPT.z         *= satMult;
}


void CompressIPT(float3& a_dataIPT, const float a_compress)
{
  // Global compress.
  float compLum = a_dataIPT.x;
  CompressWithKnee(compLum, a_compress);
  const float diff     = pow(compLum / fmax(a_dataIPT.x, 1e-6F), 3.0F - a_compress * 2.0F);

  a_dataIPT.x          = compLum;
  a_dataIPT.y         *= diff;
  a_dataIPT.z         *= diff;
  
  MoreCompressColor_IPT(a_dataIPT);
}


void ContrastIPT(float3& a_dataIPT, const float a_contrast)
{
  if (a_dataIPT.x > 0.0F && a_dataIPT.x < 1.0F)
  {
    //const float weight   = 1.0F - pow(fabs(2.0F * (a_dataIPT.x - 0.5F)), 5.0F); // save deep shadow and lights
    const float lumContr =  ContrastField(a_dataIPT.x, (a_contrast - 1.0F)/* * weight*/);

    ChangeColorWithNewLuminance_IPT(lumContr, a_dataIPT);
    MoreCompressColor_IPT(a_dataIPT);
  }

  // cylinder space.
  //float saturation = sqrt(a_dataIPT.y * a_dataIPT.y + a_dataIPT.z * a_dataIPT.z);
  //const float hue  = atan2f(a_dataIPT.z, a_dataIPT.y);

  //a_dataIPT.x      = lumContr;
  //a_dataIPT.y      = saturation * diff * cos(hue);
  //a_dataIPT.z      = saturation * diff * sin(hue);  
}


void MoreCompressColor_IPT(float3& a_dataIPT)
{
  const float saturation = sqrt(a_dataIPT.y * a_dataIPT.y + a_dataIPT.z * a_dataIPT.z);
  const float compSat    = tanh(saturation);
  const float colorDiff  = compSat / fmax(saturation, 1e-6F);

  a_dataIPT.y           *= colorDiff;
  a_dataIPT.z           *= colorDiff;
}

void ChangeColorWithNewLuminance_IPT(const float & a_newLum, float3& a_dataIPT)
{
  const float diff                     = a_newLum / fmax(a_dataIPT.x, 1e-6F);
  
  const float multColor                = fmax(diff, 1.0F);
  const float multColBlue              = fmin(diff, 1.0F); // we compress blue only when the brightness decreases, to prevent the appearance of purple.

  a_dataIPT.x                          = a_newLum;
  a_dataIPT.y                         *= multColor;               // red-green  
  if (a_dataIPT.z < 0.0F) a_dataIPT.z *= multColor * multColBlue; // blue 
  else                    a_dataIPT.z *= multColor;               // yellow
}


void ContrastRGB(float3& a_data, const float a_contrast)
{
  const float weight = a_contrast - 1.0F;
  const float gamma  = 1.0F / 2.2F; 
  a_data.x = ContrastField(pow(a_data.x, gamma), weight);
  a_data.y = ContrastField(pow(a_data.y, gamma), weight);
  a_data.z = ContrastField(pow(a_data.z, gamma), weight);

  a_data.x = pow(a_data.x, 1.0F / gamma);
  a_data.y = pow(a_data.y, 1.0F / gamma);
  a_data.z = pow(a_data.z, 1.0F / gamma);
}


float3 Mediana(const float4* inData, const int size, const int a_width, const int a_height, const int x, const int y)
{
  size_t a          = 0;
  const int sizeWin = size * size;

  std::vector<float> medianaR((size_t)sizeWin);
  std::vector<float> medianaG((size_t)sizeWin);
  std::vector<float> medianaB((size_t)sizeWin);

  const auto halfSize = (int)((float)(size - 1) / 2.0F);

  for (int i = -halfSize; i <= halfSize; ++i)
  {
    for (int j = -halfSize; j <= halfSize; ++j)
    {
      int Y = y + i;
      int X = x + j;

      if      (Y < 0)         Y = 1;
      else if (Y >= a_height) Y = a_height - 1;
      if      (X < 0)         X = 1;
      else if (X >= a_width)  X = a_width - 1;

      const int a2 = Y * a_width + X;

      medianaR[a] = inData[(size_t)a2].x;
      medianaG[a] = inData[(size_t)a2].y;
      medianaB[a] = inData[(size_t)a2].z;

      a++;
    }
  }
  sort(medianaR.begin(), medianaR.end());
  sort(medianaG.begin(), medianaG.end());
  sort(medianaB.begin(), medianaB.end());

  const auto a2 = (int)((float)(size * size) / 2.0F + 0.5F);
  return { medianaR[(size_t)a2], medianaG[(size_t)a2], medianaB[(size_t)a2] };
}

void Mediana(float inData[], const int radius, const int a_width, const int a_height, const int sizeImage)
{  
  std::vector<float> tmpData((size_t)sizeImage);

  int sizeMediana = (radius * 2 + 1) * (radius * 2 + 1);
  std::vector<float> mediana((size_t)sizeMediana);

  if (radius == 1) sizeMediana = 10;

  for (int y = 0; y < a_height; ++y)
  {
    for (int x = 0; x < a_width; ++x)
    {
      size_t a = 0;

      for (int i = -radius; i <= radius; ++i)
      {
        for (int j = -radius; j <= radius; ++j)
        {
          int Y = y + i;
          int X = x + j;
                    
          if      (Y < 0)         Y = 1;
          else if (Y >= a_height) Y = a_height - 1;
          if      (X < 0)         X = 1;
          else if (X >= a_width)  X = a_width - 1;

          const int a2 = Y * a_width + X;

          mediana[a] = inData[(size_t)a2];
          a++;
        }
      }
      sort(mediana.begin(), mediana.end());

      const int a3 = y * a_width + x;
      tmpData[(size_t)(a3)] = mediana[(size_t)(sizeMediana / 2)];
    }    
  }

#pragma omp parallel for
  for (int i = 0; i < sizeImage; ++i)
    inData[(size_t)i] = tmpData[(size_t)i];
}


void Blend(float& inData1, const float inData2, const float amount) // 0 - data1, 1 - data2
{
  inData1 = inData1 + (inData2 - inData1) * amount;
}

void Blend(float3& inData1, const float inData2, const float amount) // 0 - data1, 1 - data2
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


float ContrastField(const float a_data, const float a_amount)
{
  float result = a_data;

  if (a_data > 0.0F && a_data < 1.0F)
  {
    const float a       = a_data * 1.5707F;
    const float outData = sin(a) * sin(a);
    Blend(result, outData, a_amount);
  }

  return result;
}


void CompressWithMax(float& data, const float maxRgb)
{
  data = (data * (1.0F + data / (maxRgb * maxRgb))) / (1.0F + data);
}


void CompressWithKnee(float& a_data, const float a_compress)
{
  float knee = 10.0F;
  Blend(knee, 2.0F, pow(a_compress, 0.175F)); // lower = softer
  const float antiKnee = 1.0F / knee;  
  
  a_data = a_data / pow((1.0F + pow(a_data, knee)), antiKnee);
}


void CalculateHistogram(const float& a_data, std::vector<float>& a_histogram, const size_t a_histogramBin, const float a_iterrHistBin)
{
  if (a_data <= 1.0F)
  {
    const auto currentHistogramBin = (size_t)round(a_data * float(a_histogramBin - 1));

    if (currentHistogramBin < a_histogramBin)
      a_histogram[currentHistogramBin] += a_iterrHistBin;
  }
}


void CalculateHistogram(const float4& a_data, std::vector<float3>& a_histogram, const float a_iterrHistBin)
{
  const size_t sizeHist = a_histogram.size();

  if (a_data.x <= 1.0F)
  {
    const auto currentHistogramBin = (size_t)round(a_data.x * float(sizeHist - 1));

    if (currentHistogramBin < sizeHist)
      a_histogram[currentHistogramBin] += a_iterrHistBin;
  }

  if (a_data.y <= 1.0F)
  {
    const auto currentHistogramBin = (size_t)round(a_data.y * float(sizeHist - 1));

    if (currentHistogramBin < sizeHist)
      a_histogram[currentHistogramBin] += a_iterrHistBin;
  }

  if (a_data.z <= 1.0F)
  {
    const auto currentHistogramBin = (size_t)round(a_data.z * float(sizeHist - 1));

    if (currentHistogramBin < sizeHist)
      a_histogram[currentHistogramBin] += a_iterrHistBin;
  }
}


void CummulativeHistogram(std::vector<float>& a_histogram)
{
  // Recalculate histogramm
  for (size_t i = 1; i < a_histogram.size(); ++i)
    a_histogram[i] = a_histogram[i - 1] + a_histogram[i]; // cummulative uniform
}


void NormalizeHistogram(std::vector<float>& a_histogram)
{
  float histMin = 1e8F;
  float histMax = 0.0F;

  // Calculate min and max value (not parallel)
  for (const auto i : a_histogram)
  {
    histMin = fmin(i, histMin);
    histMax = fmax(i, histMax);
  }

  // Normalize to 0 - 1
#pragma omp parallel for
  for (int i = 0; i < (int)a_histogram.size(); ++i)
    Normalize(a_histogram[(size_t)i], histMin, histMax, 0.0F, 1.0F);
}


void NormalizeHistogram(std::vector<float3>& a_histogram)
{
  float histMin = 1e8F;
  float histMax = 0.0F;

  // Calculate min and max value (not parallel)
  for (const float3 i : a_histogram)      
    MinMaxRgb(i, histMin, histMax);       

  // Normalize to 0 - 1
#pragma omp parallel for
  for (int i = 0; i < (int)a_histogram.size(); ++i)  
    Normalize(a_histogram[(size_t)i], histMin, histMax, 0.0F, 1.0F);
}


void UniformContrastRgb(float* a_data, const size_t a_sizeImage, const size_t a_histogramBin)
{
  const size_t sizeImg3ch          = a_sizeImage * 3;
  const float iterrHistogramBin = 1.0F / (float)(sizeImg3ch);

  std::vector<float> histogram(a_histogramBin);
    
  //#pragma omp parallel for
  for (size_t i = 0; i < sizeImg3ch; ++i)
    CalculateHistogram(fmin(a_data[i], 1.0F), histogram, a_histogramBin, iterrHistogramBin);


  CummulativeHistogram(histogram);
  NormalizeHistogram  (histogram);


  // Assign a brightness histogram.
#pragma omp parallel for
  for (int i = 0; i < (int)a_sizeImage * 3; ++i)
    AssignHistogramToData(a_data[(size_t)i], histogram, a_histogramBin);
}


void AssignHistogramToData(float& a_data, const std::vector<float>& a_histogram, const size_t a_histogramBin)
{
  if (a_data <= 1.0F)
  {
    // Assign the histogram to the data->
    const auto currentHistogramBin = (size_t)round(a_data * float(a_histogramBin - 1));

    if (currentHistogramBin < a_histogramBin)
      a_data = a_histogram[currentHistogramBin];
    else
      a_data = 0.0F;
  }
}

void UniformContrast(float* a_data, const size_t a_sizeImage, const size_t a_histogramBin)
{
  const float iterrHistogramBin = 1.0F / (float)(a_sizeImage);
  std::vector<float> histogram(a_histogramBin);
  
  // not parallelized
  for (size_t i = 0; i < a_sizeImage; ++i)
    CalculateHistogram(a_data[i], histogram, a_histogramBin, iterrHistogramBin);


  CummulativeHistogram(histogram);
  NormalizeHistogram(histogram);


  // 1D blur
  const auto a_blurRadius     = (int)((float)a_histogramBin * 0.3F);
  std::vector<float> weights = CreateGaussKernelWeights1D_HDRImage2(a_blurRadius + 1);
  Blur1D(histogram, (int)a_histogramBin, a_blurRadius, weights);


  // Assign a brightness histogram.
#pragma omp parallel for
  for (int i = 0; i < (int)a_sizeImage; ++i)
    AssignHistogramToData(a_data[(size_t)i], histogram, a_histogramBin);
}


void Blur1D(std::vector<float>& a_array, const int a_sizeArray, const int a_blurRadius, std::vector<float>& a_weights)
{
  for (int i = 1; i < a_sizeArray; i++)
  {
    float summ = 0.0F;

    for (int row = -a_blurRadius; row <= a_blurRadius; ++row)
    {
      int currX = i + row;
      if      (currX < 0)            currX = int(0.0F + abs((float)row));
      else if (currX >= a_sizeArray) currX = int((float)a_sizeArray - abs((float)row));

      summ += a_array[(size_t)currX] * a_weights[(size_t)abs(row)];
    }
    a_array[(size_t)i] = summ;
  }
}


void Bilinear(const float a_inData, float* a_outData, const float a_newPosX, const float a_newPosY, const int a_width, const int a_height)
{
  // |------|------|
  // | dxy3 | dxy4 |
  // |------|------|
  // | dxy1 | dxy2 |
  // |------|------|

  const auto floorY = (int)a_newPosY;
  const auto floorX = (int)a_newPosX;

  const int     dxy1 = floorY * a_width + floorX;
  int           dxy2 = dxy1 + 1;
  if (dxy1 < 0) dxy2 = dxy1 - 1;  
  const int     dxy3 = (floorY + 1) * a_width + floorX;
  int           dxy4 = dxy3 + 1;
  if (dxy3 < 0) dxy4 = dxy3 - 1;

  float dx = a_newPosX - floor(a_newPosX);
  float dy = a_newPosY - floor(a_newPosY);

  if (dx < 0.0F) dx = 1.0F - abs(dx);
  if (dy < 0.0F) dy = 1.0F - abs(dy);

  float multBrightness1 = (1.0F - dx) * (1.0F - dy);
  float multBrightness2 = dx * (1.0F - dy);
  float multBrightness3 = dy * (1.0F - dx);
  float multBrightness4 = dx * dy;

  if (floorY >= 0 && floorX >= 0 && floorY < a_height - 1 && floorX < a_width - 1)
  {
    a_outData[(size_t)dxy1] += a_inData * multBrightness1;
    a_outData[(size_t)dxy2] += a_inData * multBrightness2;
    a_outData[(size_t)dxy3] += a_inData * multBrightness3;
    a_outData[(size_t)dxy4] += a_inData * multBrightness4;
  }
}


void Bilinear3(const float3& a_inData, float3* a_outData, const float a_newPosX, const float a_newPosY, const int a_width, const int a_height)
{
  // |------|------|
  // | dxy3 | dxy4 |
  // |------|------|
  // | dxy1 | dxy2 |
  // |------|------|

  const auto floorY = (int)a_newPosY;
  const auto floorX = (int)a_newPosX;

  const int     dxy1 = floorY * a_width + floorX;
  int           dxy2 = dxy1 + 1;
  if (dxy1 < 0) dxy2 = dxy1 - 1;
  const int     dxy3 = (floorY + 1) * a_width + floorX;
  int           dxy4 = dxy3 + 1;
  if (dxy3 < 0) dxy4 = dxy3 - 1;

  float dx = a_newPosX - floor(a_newPosX);
  float dy = a_newPosY - floor(a_newPosY);

  if (dx < 0.0F) dx = 1.0F - abs(dx);
  if (dy < 0.0F) dy = 1.0F - abs(dy);

  float multBrightness1 = (1.0F - dx) * (1.0F - dy);
  float multBrightness2 = dx * (1.0F - dy);
  float multBrightness3 = dy * (1.0F - dx);
  float multBrightness4 = dx * dy;

  if (floorY >= 0 && floorX >= 0 && floorY < a_height - 1 && floorX < a_width - 1)
  {
    a_outData[(size_t)dxy1].x += a_inData.x * multBrightness1;
    a_outData[(size_t)dxy2].x += a_inData.x * multBrightness2;
    a_outData[(size_t)dxy3].x += a_inData.x * multBrightness3;
    a_outData[(size_t)dxy4].x += a_inData.x * multBrightness4;
              
    a_outData[(size_t)dxy1].y += a_inData.y * multBrightness1;
    a_outData[(size_t)dxy2].y += a_inData.y * multBrightness2;
    a_outData[(size_t)dxy3].y += a_inData.y * multBrightness3;
    a_outData[(size_t)dxy4].y += a_inData.y * multBrightness4;
                     
    a_outData[(size_t)dxy1].z += a_inData.z * multBrightness1;
    a_outData[(size_t)dxy2].z += a_inData.z * multBrightness2;
    a_outData[(size_t)dxy3].z += a_inData.z * multBrightness3;
    a_outData[(size_t)dxy4].z += a_inData.z * multBrightness4;
  }
}


void Vignette(float3& a_data, const float a_vignette, const float a_centerImageX, const float a_centerImageY, 
  const float a_radiusImage, const int a_x, const int a_y)
{
  const float distanceFromCenter = Distance((float)a_x, (float)a_y, a_centerImageX, a_centerImageY) / a_radiusImage;
  const float vignetteIntensity  = sqrt(fmax(1.0F - distanceFromCenter * a_vignette, 0.0F));

  a_data *= vignetteIntensity;
}


void CalculateAutoWhitePoint(const float3& a_data, float3& a_autoWhitePoint)
{
  const float lum    = Luminance(a_data);  
  const float weight = WeightWhiteBalance(lum);

  a_autoWhitePoint.x += a_data.x * weight;
  a_autoWhitePoint.y += a_data.y * weight;
  a_autoWhitePoint.z += a_data.z * weight;
}


void ChrommAberr(const float3& a_inData, float* a_outData1, float* a_outData2, const int a_width, const int a_height,
  const float a_chromAberr, const int x, const int y)
{
  // Generating a velocity map.
  // (-1, -1) - left down, (0, 0) - no movement. (1, 1) - right up. Red > green > blue.

  const float stepVelocityX = (float)x / (float)a_width;
  const float stepVelocityY = (float)y / (float)a_height;
  const float velocityX = stepVelocityX * 2.0F - 1.0F; // -1.0f to 1.0f
  const float velocityY = stepVelocityY * 2.0F - 1.0F; // -1.0f to 1.0f

  const float newPosXr = (float)x + velocityX * a_chromAberr;
  const float newPosYr = (float)y + velocityY * a_chromAberr;
  const float newPosXg = (float)x + velocityX * 0.5F * a_chromAberr;
  const float newPosYg = (float)y + velocityY * 0.5F * a_chromAberr;
  
  // Anti-aliasing
  Bilinear(a_inData.x, a_outData1, newPosXr, newPosYr, a_width, a_height);
  Bilinear(a_inData.y, a_outData2, newPosXg, newPosYg, a_width, a_height);
}


void ComputeWhitePoint(const float3& a_summRgb, float3& a_whitePoint, const int a_sizeImage)
{
  a_whitePoint.x = a_summRgb.x / (float)a_sizeImage;
  a_whitePoint.y = a_summRgb.y / (float)a_sizeImage;
  a_whitePoint.z = a_summRgb.z / (float)a_sizeImage;
}


void MinMaxRgb(const float3& a_data, float& a_minRgb, float& a_maxRgb)
{
#pragma omp critical
  {
    a_minRgb = fmin(Min(a_data.x, a_data.y, a_data.z), a_minRgb);
    a_maxRgb = fmax(Max(a_data.x, a_data.y, a_data.z), a_maxRgb);
  }
}


void DrawColumn(float4* a_data, const int width, const int height, const int offsetX, const int ImgWidth, const float3& a_color, const bool transparent)
{
  const int endA = height * ImgWidth + width + offsetX;

  for (size_t y = 0; y < (size_t)height; ++y)
  {
    for (size_t x = 0; x < (size_t)width; ++x)
    {
      const size_t a = y * (size_t)ImgWidth + x + (size_t)offsetX;
      
      const float multTransp = powf((float)a / (float)endA, 4) + 0.2F;

      if (transparent)
      {
        a_data[a].x /= (2.0F + a_color.x * multTransp);
        a_data[a].y /= (2.0F + a_color.y * multTransp);
        a_data[a].z /= (2.0F + a_color.z * multTransp);
      }
      else
      {
        a_data[a].x = a_color.x;
        a_data[a].y = a_color.y;
        a_data[a].z = a_color.z;
      }
    }
  }
}


void ViewHistorgam(float4* a_data, std::vector<float3>& a_histogram, const int a_width, const int a_height)
{
  const auto  countCol = (size_t)a_width;
  const auto  sizeHist = a_histogram.size();

  const auto  step     = (size_t)((float)sizeHist / (float)countCol + 0.5F);
  const auto  width    = (int)fmax((float)a_width / (float)countCol + 0.5F, 1);  
  const float gamma    = 0.455F;
  int         height   = 0;

  float3 color(0, 0, 0);
    

  // Draw black histogram
  for (size_t i = 0; i < countCol; ++i)
  {
    const size_t j = i * step;
  
    a_histogram[j].x = powf(a_histogram[j].x, gamma);
    a_histogram[j].y = powf(a_histogram[j].y, gamma);
    a_histogram[j].z = powf(a_histogram[j].z, gamma);
  
    height         = (int)(a_histogram[j].x * (float)a_height);
    DrawColumn(a_data, width, height, (int)i * width, a_width, color, false);
  
    height         = (int)(a_histogram[j].y * (float)a_height);
    DrawColumn(a_data, width, height, (int)i * width, a_width, color, false);
  
    height         = (int)(a_histogram[j].z * (float)a_height);
    DrawColumn(a_data, width, height, (int)i * width, a_width, color, false);
  }

  // Draw color histogram
  for (size_t i = 0; i < countCol; ++i)
  {
    const size_t j = i * step;

    color          = { 0.0F, 0.0F, 0.5F };        
    height         = (int)(a_histogram[j].z * (float)a_height);
    DrawColumn(a_data, width, height, (int)i * width, a_width, color, true);
    
    color          = { 0.5F, 0.0F, 0.0F };    
    height         = (int)(a_histogram[j].x * (float)a_height);
    DrawColumn(a_data, width, height, (int)i * width, a_width, color, true);

    color          = { 0.0F, 0.5F, 0.0F };        
    height         = (int)(a_histogram[j].y * (float)a_height);
    DrawColumn(a_data, width, height, (int)i * width, a_width, color, true);
  }  
}


void MinMaxHistBin(const float* a_histogram, float& minHistBin, float& maxHistBin, const size_t a_sizeImage, const size_t a_histogramBin)
{
  size_t i = 0;
  const float floor = (float)a_sizeImage * 0.00001F; // 0.001% from image resolution

  while (i < a_histogramBin && a_histogram[i] <= floor)
  {
    minHistBin = (float)i;
    i++;
  }

  i = a_histogramBin - 1;

  while (a_histogram[i] <= floor)
  {
    maxHistBin = (float)i;
    i--;
  }

  minHistBin /= (float)a_histogramBin;
  maxHistBin /= (float)a_histogramBin;
}

void MinMaxHistBin(const std::vector<float3>& a_histogram, float& a_minAllHistBin, float& a_maxAllHistBin, const size_t a_sizeImage)
{
  size_t i              = 0;
  const float floor     = (float)a_sizeImage * 0.00001F; // 0.001% from image resolution
  const size_t sizeHist = a_histogram.size();


  while (i < sizeHist && (a_histogram[i].x <= floor || a_histogram[i].y <= floor || a_histogram[i].z <= floor))
  {
    a_minAllHistBin = (float)i;
    i++;
  }

  i = sizeHist - 1;

  while (a_histogram[i].x <= floor || a_histogram[i].y <= floor || a_histogram[i].z <= floor)
  {
    a_maxAllHistBin = (float)i;
    i--;
  }

  a_minAllHistBin /= (float)sizeHist;
  a_maxAllHistBin /= (float)sizeHist;
}


int FloatToInt(const float inData)
{   
  return int(inData + 0.5f);
}

void Resize(float* a_data, const int a_sourceWidth, const int a_sourceHeight, const int a_sourceSizeImage, const float a_resize)
{
  // Resize works correctly only if the resize is a multiple of 2.
  const auto newSizeImage = (size_t)((float)a_sourceSizeImage * (a_resize * a_resize));
  const auto newHeight    =    (int)((float)a_sourceHeight    * a_resize);
  const auto newWidth     =    (int)((float)a_sourceWidth     * a_resize);

  std::vector<float> resizeArray(newSizeImage);

  // Downsample
  if (a_resize < 1.0F)
  {
#pragma omp parallel for
    for (int y = 0; y < newHeight; ++y)
    {
      for (size_t x = 0; x < (size_t)newWidth; ++x)
      {
        const auto j = (int)((float)y / a_resize * (float)a_sourceWidth + (float)x / a_resize);
        
        if (j >= a_sourceSizeImage)
          continue;

        const size_t i = (size_t)y * (size_t)newWidth + x;
        resizeArray[i] = a_data[(size_t)j];
      }
    }


#pragma omp parallel for
    for (int y = 0; y < newHeight; ++y)
    {
      for (size_t x = 0; x < (size_t)newWidth; ++x)
      {
        const size_t i = (size_t)y * (size_t)a_sourceWidth + x;
        const size_t j = (size_t)y * (size_t)newWidth      + x;

        a_data[i] = resizeArray[j];
      }
    }
  }
  else // Upsample
  {
#pragma omp parallel for
    for (int y = 0; y < newHeight; ++y)
    {
      for (size_t x = 0; x < (size_t)newWidth; ++x)
      {
        const size_t i = (size_t)y * (size_t)newWidth + x;
        const auto   j = (int)((float)y / a_resize * (float)newWidth + (float)x / a_resize);
        resizeArray[i] = a_data[(size_t)j];
      }
    }


#pragma omp parallel for
  for (int i = 0; i < (int)newSizeImage; ++i)
    a_data[(size_t)i] = resizeArray[(size_t)i];
  }
}

void Sharp(float4* a_image4out, const std::vector<float>& a_lumForSharp, const float a_sharpness, const int a_width, const int a_height)
{
#pragma omp parallel for
  for (int y = 0; y < a_height; ++y)
  {
    for (size_t x = 0; x < (size_t)a_width; ++x)
    {
      const size_t a = (size_t)y * (size_t)a_width + x;

      a_image4out[a].x = fmax(a_image4out[a].x, 0.0F);
      a_image4out[a].y = fmax(a_image4out[a].y, 0.0F);
      a_image4out[a].z = fmax(a_image4out[a].z, 0.0F);
            
      float mean = 0.0F;

      for (int i = -1; i <= 1; ++i)
      {
        for (int j = -1; j <= 1; ++j)
        {
          int Y =      y + i;
          int X = (int)x + j;

          if      (Y < 0)         Y = 1;
          else if (Y >= a_height) Y = a_height - 1;
          if      (X < 0)         X = 1;
          else if (X >= a_width)  X = a_width - 1;

          const size_t a2 = (size_t)Y * (size_t)a_width + (size_t)X;

          mean += a_lumForSharp[a2];
        }
      }

      mean /= 9.0F;

      float dispers = 0.0F;

      for (int i = -1; i <= 1; i++)
      {
        for (int j = -1; j <= 1; j++)
        {
          int Y =      y + i;
          int X = (int)x + j;

          if      (Y < 0)         Y = 1;
          else if (Y >= a_height) Y = a_height - 1;
          if      (X < 0)         X = 1;
          else if (X >= a_width)  X = a_width - 1;

          const size_t a2 = (size_t)Y * (size_t)a_width + (size_t)X;

          const float val = a_lumForSharp[a2] - mean;
          dispers        += abs(val);
        }
      }

      dispers            /= (1.0F + dispers);
      float hiPass        = (a_lumForSharp[a] - mean) * 5 + 1.0F;
      hiPass              = powf(hiPass, 2);
                        
      float sharp         = 1.0F;
      Blend(sharp, hiPass, a_sharpness * (1.0F - dispers));

      a_image4out[a].x     *= sharp;
      a_image4out[a].y     *= sharp;
      a_image4out[a].z     *= sharp;
    }
  }
}

float WeightWhiteBalance(const float a_lum)
{
  return (100.0F + a_lum) / (100.0F + a_lum * a_lum);
}


std::vector<float> CreateGaussKernelWeights1D_HDRImage2(const int blurRadius)
{
  std::vector<float> gKernel(blurRadius);

  // sum is for normalization
  float sum = 0.0F;

  for (size_t x = 0; x < (size_t)blurRadius; ++x)
  {
    gKernel[x] = exp(-(float)x * (float)x / ((float)blurRadius * (float)blurRadius / 6.28318530F));
    sum += gKernel[x];
  }

  // normalize the Kernel
  for (int i = 0; i < blurRadius; ++i)
    gKernel[(size_t)i] /= (2.0F * sum);

  return gKernel;
}


void Blur(float* a_data, int a_blurRadius, const int a_width, const int a_height, const int a_sizeImage)
{
  const float diagonalImage = Distance(0.0F, 0.0F, (float)a_width, (float)a_height);
  const float radiusImage   = diagonalImage / 2.0F;

  // Resize works correctly only if the resize is a multiple of 2.
  float resize = 1.0F;
  if      (a_blurRadius >= (int)(radiusImage / 4.0F)) resize = 0.5F;
  else if (a_blurRadius >= (int)(radiusImage / 2.0F)) resize = 0.25F;
  a_blurRadius = (int)((float)a_blurRadius * resize);

  const int resizeWidth      = int((float)a_width     * resize);
  const int resizeHeight     = int((float)a_height    * resize);
  const int resizeSizeImage  = int((float)a_sizeImage * (resize * resize));

  std::vector<float> blurPass(resizeSizeImage);
  std::vector<float> weights = CreateGaussKernelWeights1D_HDRImage2(a_blurRadius+1);

  // Downsize for speed up blur
  if (resize < 1.0F) 
    Resize(a_data, a_width, a_height, a_sizeImage, resize);

  //------------ Blur ------------

  float summ;
  float multRadius = 1.0F;
  if (a_blurRadius >= resizeWidth)
    multRadius = (float)(resizeWidth - 1) / (float)a_blurRadius;
  
  // All row.  
  
// cannot be parallelize.
  for (int y = 0; y < resizeHeight; ++y)
  {
    for (int x = 0; x < resizeWidth; ++x)
    {
      summ = 0.0F;      

      for (int row = -a_blurRadius; row <= a_blurRadius; ++row)
      {
        int currX = x + row;

        if      (currX < 0)            currX = 0           + int(abs((float)row) * multRadius);
        else if (currX >= resizeWidth) currX = resizeWidth - int(abs((float)row) * multRadius);
      
        const size_t a = (size_t)y * (size_t)a_width + (size_t)currX;

        summ += a_data[a] * weights[(size_t)abs(row)];
      }
      const size_t a = (size_t)y * (size_t)resizeWidth + (size_t)x;
      blurPass[a] = summ;
    }
  }

  // All col.

  multRadius = 1.0F;
  if (a_blurRadius >= resizeHeight)
    multRadius = (float)(resizeHeight-1) / (float)a_blurRadius;

#pragma omp parallel for
  for (int x = 0; x < resizeWidth; ++x)
  {
    for (int y = 0; y < resizeHeight; ++y)
    {
      summ = 0.0F;

      for (int col = -a_blurRadius; col <= a_blurRadius; ++col)
      {
        int currY = y + col;        

        if      (currY < 0)             currY = 0            + (int)(abs((float)col) * multRadius);
        else if (currY >= resizeHeight) currY = resizeHeight - (int)(abs((float)col) * multRadius);
        
        const size_t a = (size_t)currY * (size_t)resizeWidth + (size_t)x;
        summ += blurPass[a] * weights[abs(col)];
      }
      const size_t a = (size_t)y * (size_t)a_width + (size_t)x;
      a_data[a] = summ;
    }
  }    
  
  // Upsize
  resize = 1.0F / resize;
  if (resize > 1.0F) Resize(a_data, resizeWidth, resizeHeight, resizeSizeImage, resize);
}


[[nodiscard]] float ConvertAngleToRad(const int angle)
{
  return (float)angle * 0.01745329252F;
}

void DiffractionStars(const float3& a_data, float3& a_diffrStars, const float a_sizeStar, const int a_numRay, const int a_RotateRay,
 const float a_randomAngle,  const float a_sprayRay, const int m_width, const int m_height, const float radiusImage, const int x, const int y)
{  
  //const float PI          = 3.141592F;
  const float twoPI       = 6.283185F;
  const float angle       = twoPI / (float)a_numRay;
  //const float waveLenghtR = 0.000700F; 
  //const float waveLenghtY = 0.000575F; 
  const float waveLenghtG = 0.000530F; 
  //const float waveLenghtC = 0.000485F; 
  //const float waveLenghtB = 0.000450F;   
  //const float waveLenghtV = 0.000400F;   

  float lum                 = Luminance(a_data);
  lum                      /= (1.0F + lum);   
  
  // Ray
  for (int numRay = 0; numRay < a_numRay; ++numRay)
  {
    float jitter                = (rand() % 100) / 100.0F + 0.5F;
    float nextAngle             = angle * numRay - ConvertAngleToRad(a_RotateRay);

    if ((rand() % 1000) / 1000.0F <= a_randomAngle)    
      nextAngle = (rand() % 628) / 100.0F;
    
    const float sprayAngle      = ((float)(rand() % (int)(angle * 100)) / 100.0F) * a_sprayRay;
    nextAngle                  += (sprayAngle - angle / 2.0F);

    const float sizeStar        = a_sizeStar * radiusImage / 100.0F;
    const float sizeStarFromLum = sizeStar * lum * jitter / 5.0F;

    // Aperture big   - dist ray small, wave small
    // Aperture small - dist ray big, wave big
    const float diafragma = 1.0F / (sizeStarFromLum * 5.0F * jitter);

    for (float sample = 1.0F; sample < sizeStarFromLum; sample += 0.1F)
    {
      const float angleFi  = sample;
      const float temp     = diafragma * angleFi;
      const float multDist = 1.0F - sample / sizeStarFromLum;
      const float u        = temp / waveLenghtG + 0.000001F;

      const float newX     = (float)x + cos(nextAngle) * sample;
      const float newY     = (float)y + sin(nextAngle) * sample;

      if (newX > 0 && newX < (float)m_width && newY > 0 && newY < (float)m_height)
      {
        const float I = powf((sinf(u) / u), 2) * multDist;
        Bilinear3(a_data * I, &a_diffrStars, newX, newY, m_width, m_height);
      }
    }
  }
}

void UniformContrastRGBFilter(float4* a_data, const size_t a_sizeImage, const size_t a_histogramBin, const float a_amount)
{  
  const size_t size3image = a_sizeImage * 3;
  std::vector<float> rgbArray(size3image);

  // Convert 3 field RGB to linear array.
#pragma omp parallel for
  for (int i = 0; i < (int)a_sizeImage; ++i)
  {
    const auto a           = (size_t)i;
    const size_t indColorG = a + a_sizeImage;
    const size_t indColorB = a + a_sizeImage * 2;

    rgbArray[a]            = a_data[a].x;
    rgbArray[indColorG]    = a_data[a].y;
    rgbArray[indColorB]    = a_data[a].z;
  }


  UniformContrastRgb(&rgbArray[0], a_sizeImage, a_histogramBin);


  // Return to main array
#pragma omp parallel for
  for (int i = 0; i < (int)a_sizeImage; ++i)
  {
    const auto   a            = (size_t)i;
    const size_t indColorG    = a + a_sizeImage;
    const size_t indColorB    = a + a_sizeImage * 2;

    rgbArray[a]               = pow(rgbArray[a], 2.2F);
    rgbArray[indColorG]       = pow(rgbArray[indColorG], 2.2F);
    rgbArray[indColorB]       = pow(rgbArray[indColorB], 2.2F);

    const float meanRGBsource = (a_data[a].x + a_data[a].y + a_data[a].z) / 3.0F;
    const float meanRGB_UC    = (rgbArray[a] + rgbArray[indColorG] + rgbArray[indColorB]) / 3.0F;

    const float diff          = meanRGB_UC / fmax(meanRGBsource, 1e-6F);

    float3 rgbDiff(a_data[a].x, a_data[a].y, a_data[a].z);    

    Blend(rgbDiff.x, a_data[a].x * diff, a_amount);
    Blend(rgbDiff.y, a_data[a].y * diff, a_amount);
    Blend(rgbDiff.z, a_data[a].z * diff, a_amount);

    Blend(a_data[a].x, rgbArray[a]        , a_amount);
    Blend(a_data[a].y, rgbArray[indColorG], a_amount);
    Blend(a_data[a].z, rgbArray[indColorB], a_amount);

    Blend(a_data[a].x, rgbDiff.x, 0.5F);
    Blend(a_data[a].y, rgbDiff.y, 0.5F);
    Blend(a_data[a].z, rgbDiff.z, 0.5F);
  }
}

void UniformContrastIPTFilter(float4* a_data, const size_t a_sizeImage, const size_t a_histogramBin, const float a_amount)
{
  std::vector<float3> IPT_array    (a_sizeImage);
  std::vector<float>  IPT_lum_array(a_sizeImage);

#pragma omp parallel for
  for (int i = 0; i < (int)a_sizeImage; ++i)
  {
    const auto a       = (size_t)i;
    IPT_array[a]       = float3(a_data[a].x, a_data[a].y, a_data[a].z);
    ConvertSrgbToXyz    (IPT_array[a]);
    ConvertXyzToLmsPower(IPT_array[a], 0.43F);
    ConvertLmsToIpt     (IPT_array[a]);
    IPT_lum_array[a]   = IPT_array[a].x;
  }


  UniformContrast(&IPT_lum_array[0], a_sizeImage, a_histogramBin);
  

  // Return to main array
#pragma omp parallel for
  for (int i = 0; i < (int)a_sizeImage; ++i)
  {
    const auto a               = (size_t)i;
    const float diff           = IPT_lum_array[a] / fmax(IPT_array[a].x, 1e-6F);
    const float diffCompress   = tanh(diff - 1.0F); // hyperbolic tangent for smooth limit shadow and light.
    const float lum            = IPT_array[a].x * (1.0F + diffCompress);

    ChangeColorWithNewLuminance_IPT(lum, IPT_array[a]);
    MoreCompressColor_IPT(IPT_array[a]);

    ConvertIptToLms     (IPT_array[a]);
    ConvertLmsToXyzPower(IPT_array[a], 1.0F / 0.43F);
    ConvertXyzToSrgb    (IPT_array[a]);
    ClampMinusToZero    (IPT_array[a]);

    // little compress in RGB
    CompressWithKnee(IPT_array[a].x, 0.01F);
    CompressWithKnee(IPT_array[a].y, 0.01F);
    CompressWithKnee(IPT_array[a].z, 0.01F);    

    Blend(a_data[a].x, IPT_array[a].x, a_amount);
    Blend(a_data[a].y, IPT_array[a].y, a_amount);
    Blend(a_data[a].z, IPT_array[a].z, a_amount);
  }
}


void NormalizeFilter(float4* a_data, std::vector<float3>& a_histogram, const size_t a_sizeImage, const size_t a_histogramBin,
 float3& minHistBin, float3& maxHistBin, float& minAllHistBin, float& maxAllHistBin, const float amount)
{
// not parallelized
  for (size_t i = 0; i < a_sizeImage; ++i)
    CalculateHistogram(a_data[i], a_histogram, 1.0F);


  // Calculate min/max a_histogram for a_normalize
  MinMaxHistBin(a_histogram, minAllHistBin, maxAllHistBin, a_sizeImage);
  

  // Normalize 
#pragma omp parallel for
  for (int i = 0; i < (int)a_sizeImage; ++i)
  {
    const auto a    = (size_t)i;
    float3 dataNorm = float3(a_data[a].x, a_data[a].y, a_data[a].z);

    Normalize(dataNorm, minAllHistBin, maxAllHistBin, 0.0F, 1.0F);
    
    Blend(a_data[a].x, dataNorm.x, amount);
    Blend(a_data[a].y, dataNorm.y, amount);
    Blend(a_data[a].z, dataNorm.z, amount);
  }
}

void ViewHistorgamFilter(float4* a_data, std::vector<float3>& a_histogram, const int a_width, const int a_height, const bool a_viewHistogram)
{
  if (!a_viewHistogram)
    return;

  const int sizeImage    = a_width * a_height;
  const size_t histSize  = a_histogram.size();

#pragma omp parallel for
  for (int i = 0; i < sizeImage; ++i)
    CalculateHistogram(a_data[(size_t)i], a_histogram, 1.0F);


  NormalizeHistogram(a_histogram);
  ViewHistorgam(a_data, a_histogram, a_width, a_height);
}


void CalculateStatisticsImage(const float4* a_data, float& a_minRgb, float& a_maxRgb, std::vector<float3>& a_histogram, const int a_width, const int a_height)
{
  const auto histogramBin = (int)a_histogram.size();
  const int sizeImage    = a_width * a_height;
  a_minRgb               = 1e6F;
  a_maxRgb               = 0.0F;

// not parallelized
  for (size_t i = 0; i < (size_t)sizeImage; ++i)
  {
    const float3 color(a_data[i].x, a_data[i].y, a_data[i].z);
    MinMaxRgb(color, a_minRgb, a_maxRgb);
    CalculateHistogram(a_data[i], a_histogram, 1.0F);
  }
}


void PrintStatistics(const bool a_statistics, const float a_endTime, const float a_minRgb, const float a_maxRgb, const float3& a_whitePointColor)
{
  if (!a_statistics)
    return;

  std::cout << std::endl;
  std::cout << "endTime = "    << a_endTime / 1000.0F << std::endl;
  std::cout << "maxThreads = " << omp_get_num_procs() << std::endl;
  std::cout << "minRgb = "     << a_minRgb            << std::endl;
  std::cout << "maxRgb = "     << a_maxRgb            << std::endl;
  std::cout << "whitePoint = " << a_whitePointColor.x << " " << a_whitePointColor.y << " " << a_whitePointColor.z << std::endl;
  std::cout << std::endl;
}
