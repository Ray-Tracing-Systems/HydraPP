#include "PostProcess.h"

#include "../../HydraAPI/hydra_api/pugixml.hpp"
#include "../../HydraAPI/hydra_api/HydraXMLHelpers.h"

#include "../hydra_pp/HydraPostProcessCommon.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class PostProcessHydra1 : public IFilter2D
{
public: 
  void Release() override;   // COM like destructor; 
  bool Eval()    override;   //
};


extern "C" __declspec(dllexport) IFilter2D* CreateFilter(const wchar_t* a_filterName)
{
  std::wstring inFilterName(a_filterName);

  if (inFilterName == L"post_process_hydra1")
    return new PostProcessHydra1;
  else
    return nullptr;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PostProcessHydra1::Release()
{   
  // nothing to destroy here, we have simple implementation
}

void ExecutePostProcessHydra1(
  const float* a_input, float* a_output, unsigned int a_numThreads,
  const float a_exposure, const float a_compress, const float a_contrast,
  const float a_saturation, const float a_whiteBalance, float3 a_whitePointColor,
  const float a_uniformContrast, const float a_normalize, const float a_vignette,
  const float a_chromAberr, const float a_sharpness, const float a_sizeStar,
  const int a_numRay, const int a_rotateRay, const float a_randomAngle,
  const float a_sprayRay, unsigned int a_width, unsigned int a_height); // forward declaraion

bool PostProcessHydra1::Eval()
{
  // read presets
  //
  pugi::xml_document docXml;
  docXml.load_string(GetSettingsStr());
  pugi::xml_node settings = docXml.child(L"settings");

  // read input and output arguments
  //
  unsigned int numThreads = settings.attribute(L"numThreads").as_int();

  float exposure          = settings.attribute(L"exposure").as_float();
  float compress          = settings.attribute(L"compress").as_float();
  float contrast          = settings.attribute(L"contrast").as_float();
  float saturation        = settings.attribute(L"saturation").as_float();
  float whiteBalance      = settings.attribute(L"whiteBalance").as_float();
  float3 whitePointColor  = HydraXMLHelpers::ReadFloat3(settings.attribute(L"whitePointColor"));
  float uniformContrast   = settings.attribute(L"uniformContrast").as_float();
  float normalize         = settings.attribute(L"normalize").as_float();
  float vignette          = settings.attribute(L"vignette").as_float();
  float chromAberr        = settings.attribute(L"chromAberr").as_float();
  float sharpness         = settings.attribute(L"sharpness").as_float();

  // Diffraction stars
  float sizeStar          = settings.attribute(L"diffStars_sizeStar").as_float(); // 0-100
  int   numRay            = settings.attribute(L"diffStars_numRay").as_int();  // 0-16
  int   rotateRay         = settings.attribute(L"diffStars_rotateRay").as_int();  // 0-360
  float randomAngle       = settings.attribute(L"diffStars_randomAngle").as_float();   // 0-1
  float sprayRay          = settings.attribute(L"diffStars_sprayRay").as_float();   // 0-1 

  int w1, h1, bpp1;
  int w2, h2, bpp2;

  const float* input      = GetInputByName(L"in_color"  , &w1, &h1, &bpp1);
  float*       output     = GetOutputByName(L"out_color", &w2, &h2, &bpp2);

  // check we have correct in and out arguments
  //
  if (exposure < 0.0f || compress < 0.0f || contrast < 0.0f ||
    saturation < 0.0f || whiteBalance < 0.0f || whitePointColor.x < 0.0f || 
    whitePointColor.y < 0.0f || whitePointColor.z < 0.0f || uniformContrast < 0.0f || 
    normalize < 0.0f || vignette < 0.0f || chromAberr < 0.0f || sharpness < 0.0f ||
    sizeStar < 0.0f || numRay < 0.0f || rotateRay < 0.0f || randomAngle < 0.0f || 
    sprayRay < 0.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; Arguments must be greater than zero.", ERR_MSG_SIZE);
    return false;
  }

  if (exposure > 10.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; exposure should be in the range 0 - 10. Default 1. This will be limited to 10.", ERR_MSG_SIZE);
    exposure = 10.0f;
  }
  if (compress > 1.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; compress should be in the range 0 - 1. Default 0. This will be limited to 1.", ERR_MSG_SIZE);
    compress = 1.0f;
  }
  if (contrast < 1.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; contrast should be in the range 1 - 2. Default 1. This will be limited to 1.", ERR_MSG_SIZE);
    contrast = 1.0f;
  }
  else if (contrast > 2.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; contrast should be in the range 1 - 2. Default 1. This will be limited to 2.", ERR_MSG_SIZE);
    contrast = 2.0f;
  }
  if (saturation > 2.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; saturation should be in the range 0 - 2. Default 1. This will be limited to 2.", ERR_MSG_SIZE);
    saturation = 2.0f;
  }
  if (whiteBalance > 1.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; whiteBalance should be in the range 0 - 1. Default 0. This will be limited to 1.", ERR_MSG_SIZE);
    whiteBalance = 1.0f;
  }
  if (uniformContrast > 1.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; uniformContrast should be in the range 0 - 1. Default 0. This will be limited to 1.", ERR_MSG_SIZE);
    uniformContrast = 1.0f;
  }
  if (normalize > 1.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; normalize should be in the range 0 - 1. Default 0. This will be limited to 1.", ERR_MSG_SIZE);
    normalize = 1.0f;
  }
  if (vignette > 1.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; vignette should be in the range 0 - 1. Default 0. This will be limited to 1.", ERR_MSG_SIZE);
    vignette = 1.0f;
  }
  if (chromAberr > 1.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; chromAberr should be in the range 0 - 1. This will be limited to 1.", ERR_MSG_SIZE);
    chromAberr = 1.0f;
  }
  if (sharpness > 1.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; sharpness should be in the range 0 - 1. This will be limited to 1.", ERR_MSG_SIZE);
    sharpness = 1.0f;
  }
  if (sizeStar > 100.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; sizeStar should be in the range 0 - 100. This will be limited to 100.", ERR_MSG_SIZE);
    sizeStar = 100.0f;
  }
  if (numRay > 16)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; numRay should be in the range 0 - 16. This will be limited to 16.", ERR_MSG_SIZE);
    numRay = 16;
  }
  if (rotateRay > 360)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; rotateRay should be in the range 0 - 360. This will be limited to 360.", ERR_MSG_SIZE);
    rotateRay = 360;
  }
  if (randomAngle > 1.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; randomAngle should be in the range 0 - 1. This will be limited to 1.", ERR_MSG_SIZE);
    randomAngle = 1.0f;
  }
  if (sprayRay > 1.0f)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; sprayRay should be in the range 0 - 1. This will be limited to 1.", ERR_MSG_SIZE);
    sprayRay = 1.0f;
  }

  if (w1 != w2 || h1 != h2 || w1 <= 0 || h1 <= 0)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; bad input size", ERR_MSG_SIZE);
    return false;
  }
  if (bpp1 != bpp2 || bpp1 != 16)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; ivalid image format; both images must be HDR;", ERR_MSG_SIZE);
    return false;
  }
  if (input == nullptr)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; argument not found: 'in_color' ", ERR_MSG_SIZE);
    return false;
  }
  if (output == nullptr)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; argument not found: 'out_color' ", ERR_MSG_SIZE);
    return false;
  }

  // now run our internal implementation
  //

  ExecutePostProcessHydra1(
    input, output, numThreads, 
    exposure, compress, contrast, 
    saturation, whiteBalance, whitePointColor,
    uniformContrast, normalize, vignette,
    chromAberr, sharpness, sizeStar, 
    numRay, rotateRay, randomAngle, 
    sprayRay, w1, h1);

  return true;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ExecutePostProcessHydra1(
  const float* a_input, float* a_output, unsigned int a_numThreads,
  const float a_exposure, const float a_compress, const float a_contrast,
  const float a_saturation, const float a_whiteBalance, float3 a_whitePointColor,
  const float a_uniformContrast, const float a_normalize, const float a_vignette,
  const float a_chromAberr, const float a_sharpness, const float a_sizeStar,
  const int a_numRay, const int a_rotateRay, const float a_randomAngle, 
  const float a_sprayRay, unsigned int a_width, unsigned int a_height)
{ 
  // Set the desired number of threads
  const unsigned int maxThreads = omp_get_num_procs();
  if (a_numThreads > maxThreads) a_numThreads = maxThreads;
  else if (a_numThreads <= 0)    a_numThreads = maxThreads;
  omp_set_num_threads(a_numThreads);

  float4* image4in = (float4*)a_input;
  float4* image4out = (float4*)a_output;


  // variables and constants
  const int sizeImage = a_width * a_height;
  const float centerImageX = a_width / 2.0f;
  const float centerImageY = a_height / 2.0f;
  const float diagonalImage = Distance(0, 0, (float)a_width, (float)a_height);
  const float radiusImage = diagonalImage / 2.0f;

  const int histogramBin = 10000;

  bool autoWhiteBalance = true;
  if (a_whitePointColor.x > 0.0f || a_whitePointColor.y > 0.0f || a_whitePointColor.z > 0.0f)
  {
    autoWhiteBalance = false;
  }
  float3 summRgb = { 0.0f, 0.0f, 0.0f };
  float3 d65 = { 0.9505f , 1.0f, 1.0888f }; // in XYZ

  float  minRgbSource = 10000.0f;
  float  maxRgbSource = 0.0f;

  float  minAllHistBin = 0.0f;
  float  maxAllHistBin = 1.0f;
  float3 minHistBin(0.0f, 0.0f, 0.0f);
  float3 maxHistBin(histogramBin, histogramBin, histogramBin);

  // arrays
  float* chromAbberR = nullptr;
  float* chromAbberG = nullptr;
  float* lumForSharp = nullptr;
  float3* diffrStars = nullptr;
  float3* histogram = nullptr;


  if (a_chromAberr > 0.0f)
  {
    chromAbberR = new float[sizeImage]; 
    chromAbberG = new float[sizeImage]; 
  }

  if (a_normalize > 0.0f)
    histogram = new float3[histogramBin]; 

  if (a_sharpness > 0.0f)
    lumForSharp = new float[sizeImage];

  if (a_sizeStar > 0.0f)
    diffrStars = new float3[sizeImage];


  //////////////////////////////////////////////////////////////////////////////////////


  // ----- Analization, precalculate and any effects. Loop 1. -----
  //#pragma omp parallel for
  for (int y = 0; y < a_height; ++y)
  {
    for (int x = 0; x < a_width; ++x)
    {
      const int i = y * a_width + x;

      image4out[i] = image4in[i];  

      // ----- Sharpness -----
      if (a_sharpness > 0.0f)
      {
        lumForSharp[i] = (image4out[i].x + image4out[i].y + image4out[i].z) / 3.0f;
        lumForSharp[i] /= (1.0f + lumForSharp[i]);
      }

      // ----- Exposure -----
      if (a_exposure != 1.0f)
      {
        image4out[i].x *= a_exposure;
        image4out[i].y *= a_exposure;
        image4out[i].z *= a_exposure;
      }

      // Minimum, maximum, negative number clamp and NaN.
      MinMaxRgb(image4out, minRgbSource, maxRgbSource, 0.0f, i);

      // Collect all the brightness values by channel.
      if (a_whiteBalance > 0.0f && autoWhiteBalance && image4out[i].x < 1.0f && image4out[i].y < 1.0f && image4out[i].z < 1.0f)
        SummValueOnField(image4out, &summRgb, i);
      

      // ----- Vignette -----
      if (a_vignette > 0.0f)      
        Vignette(image4out, a_width, a_height, a_vignette, diagonalImage, centerImageX, centerImageY, radiusImage, x, y, i);
      

      // ----- Difraction stars -----
      if (a_sizeStar > 0.0f && image4out[i].x > 50.0f || image4out[i].y > 50.0f || image4out[i].z > 50.0f)      
        DiffractionStars(image4out, diffrStars, a_sizeStar, a_numRay, a_rotateRay, a_randomAngle, a_sprayRay, a_width, a_height, sizeImage, radiusImage, x, y, i);
      

      // ----- Chromatic aberration -----
      if (a_chromAberr > 0.0f)
        ChrommAberr(image4out, chromAbberR, chromAbberG, a_width, a_height, sizeImage, a_chromAberr, x, y, i);      
    }
  }


  // ----- Chromatic aberration -----
  if (a_chromAberr > 0.0f)
  {
    #pragma omp parallel for
    for (int i = 0; i < sizeImage; ++i)
    {
      image4out[i].x = chromAbberR[i];
      image4out[i].y = chromAbberG[i];
    }
  }


  // ----- Sharpness -----
  if (a_sharpness > 0.0f)
    Sharp(image4out, lumForSharp, a_sharpness, a_width, a_height);


  // ----- White point for white balance -----
  if (a_whiteBalance > 0.0f)
  {
    if (autoWhiteBalance)
      ComputeWhitePoint(summRgb, &a_whitePointColor, sizeImage);

    ConvertSrgbToXyz(&a_whitePointColor);
    MatrixCat02(&a_whitePointColor);
    MatrixCat02(&d65);

    const float lum = Luminance(a_whitePointColor);

    if (lum > 0.0f)
    {
      a_whitePointColor.x /= lum;
      a_whitePointColor.y /= lum;
      a_whitePointColor.z /= lum;
    }
  }

  
  // ---------- Many filters Loop 2. ----------
  #pragma omp parallel for
  for (int i = 0; i < sizeImage; ++i)
  {
    // ----- Diffraction stars -----
    if (a_sizeStar > 0.0f)
    {
      image4out[i].x += diffrStars[i].x;
      image4out[i].y += diffrStars[i].y;
      image4out[i].z += diffrStars[i].z;
    }

    // ----- White balance -----
    if (a_whiteBalance > 0.0f)
    {
      ConvertSrgbToXyz(&image4out[i]);
      ChromAdaptIcam(&image4out[i], a_whitePointColor, d65, a_whiteBalance);
      ConvertXyzToSrgb(&image4out[i]);

      ClampMinusToZero(&image4out[i]);
    }

    // ----- Saturation  -----
    if (a_saturation != 1.0f)
    {
      const float lum = Luminance(image4out[i]);

      // Return to main array
      image4out[i].x = lum + (image4out[i].x - lum) * a_saturation;
      image4out[i].y = lum + (image4out[i].y - lum) * a_saturation;
      image4out[i].z = lum + (image4out[i].z - lum) * a_saturation;

      ClampMinusToZero(image4out, i);
    }

    // ----- Compress -----
    if (a_compress > 0.0f && maxRgbSource > 1.01f)
    {
      //Global LMS/IPT compress.
      float knee = 10.0f;
      Blend(knee, 2.0f, pow(a_compress, 0.175f)); // lower = softer
      const float antiKnee = 1.0f / knee;

      ConvertSrgbToXyz(&image4out[i]);
      ConvertXyzToLmsPower(&image4out[i], 0.43f);
      ConvertLmsToIpt(&image4out[i]);

      const float compLum = image4out[i].x / pow((1.0f + pow(image4out[i].x, knee)), antiKnee);
      const float multSat = pow(compLum / image4out[i].x, 1.5f);
      image4out[i].x = compLum;
      image4out[i].y *= multSat;
      image4out[i].z *= multSat;

      ConvertIptToLms(&image4out[i]);
      ConvertLmsToXyzPower(&image4out[i]);
      ConvertXyzToSrgb(&image4out[i]);
    }


    // ----- Contrast -----
    if (a_contrast > 1.0f)
    {
      //RGB
      float4 contrRGB = image4out[i];
      contrRGB.x = pow(contrRGB.x, 1.0f / 2.2f);
      contrRGB.y = pow(contrRGB.y, 1.0f / 2.2f);
      contrRGB.z = pow(contrRGB.z, 1.0f / 2.2f);
      ContrastField(contrRGB.x, a_contrast - 1.0f);
      ContrastField(contrRGB.y, a_contrast - 1.0f);
      ContrastField(contrRGB.z, a_contrast - 1.0f);
      contrRGB.x = pow(contrRGB.x, 2.2f);
      contrRGB.y = pow(contrRGB.y, 2.2f);
      contrRGB.z = pow(contrRGB.z, 2.2f);

      //IPT
      ConvertSrgbToXyz(&image4out[i]);
      ConvertXyzToLmsPower(&image4out[i], 0.43f);
      ConvertLmsToIpt(&image4out[i]);

      float lumContr = image4out[i].x;
      ContrastField(lumContr, a_contrast - 1.0f);

      const float multColor = 1.0f + abs(image4out[i].x - lumContr);
      image4out[i].y *= multColor;
      image4out[i].z *= multColor;
      image4out[i].x = lumContr;

      ConvertIptToLms(&image4out[i]);
      ConvertLmsToXyzPower(&image4out[i]);
      ConvertXyzToSrgb(&image4out[i]);

      // Blend RGB and IPT 50x50%
      Blend(image4out[i].x, contrRGB.x, 0.5f);
      Blend(image4out[i].y, contrRGB.y, 0.5f);
      Blend(image4out[i].z, contrRGB.z, 0.5f);
    }
  }
  

  // ----- Uniform contrast. Loop 3 and 4. -----
  if (a_uniformContrast > 0.0f)
  {
    UniformContrastRGBFilter(&image4out[0], sizeImage, histogramBin, a_uniformContrast);
  }

  // ----- Calculate histogram for normalize. Loop 5 and 6. -----
  if (a_normalize > 0.0f)
    NormalizeFilter(&image4out[0], &histogram[0], sizeImage, histogramBin, minHistBin, maxHistBin, minAllHistBin, maxAllHistBin, a_normalize);


  //********************* End filters *********************

  if (a_chromAberr > 0.0f)
  {
    delete[] chromAbberR;
    delete[] chromAbberG;
  }

  if (a_sharpness > 0.0f) 
    delete[] lumForSharp;

  if (a_sizeStar > 0.0f)  
    delete[] diffrStars;
      
  if (a_normalize > 0.0f )
    delete[] histogram;  
}




