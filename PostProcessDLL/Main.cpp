#include "Functions.h"

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
  const float a_exposure, const float a_compress, const float a_contrast, const float a_saturation, const float a_vibrance,
  const float a_whiteBalance, float3 a_whitePointColor, const float a_uniformContrast, const float a_normalize, const float a_vignette,
  const float a_chromAberr, const float a_sharpness, const float a_sizeStar, const int a_numRay, const int a_rotateRay, 
  const float a_randomAngle, const float a_sprayRay, unsigned int a_width, unsigned int a_height); // forward declaraion

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
  float vibrance          = settings.attribute(L"vibrance").as_float();
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
  if (exposure < 0.0F || compress < 0.0F || contrast < 0.0F || saturation < 0.0F || vibrance < 0.0F || whiteBalance < 0.0F ||
    whitePointColor.x < 0.0F || whitePointColor.y < 0.0F || whitePointColor.z < 0.0F || uniformContrast < 0.0F || normalize < 0.0F ||
    vignette < 0.0F || chromAberr < 0.0F || sharpness < 0.0F || sizeStar < 0.0F || numRay < 0.0F || rotateRay < 0.0F || 
    randomAngle < 0.0F || sprayRay < 0.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; Arguments must be greater than zero.", ERR_MSG_SIZE);
    return false;
  }

  if (exposure > 10.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; exposure should be in the range 0 - 10. Default 1. This will be limited to 10.", ERR_MSG_SIZE);
    exposure = 10.0F;
  }
  if (compress > 1.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; compress should be in the range 0 - 1. Default 0. This will be limited to 1.", ERR_MSG_SIZE);
    compress = 1.0F;
  }
  if (contrast < 1.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; contrast should be in the range 1 - 2. Default 1. This will be limited to 1.", ERR_MSG_SIZE);
    contrast = 1.0F;
  }
  else if (contrast > 2.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; contrast should be in the range 1 - 2. Default 1. This will be limited to 2.", ERR_MSG_SIZE);
    contrast = 2.0F;
  }
  if (saturation > 2.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; saturation should be in the range 0 - 2. Default 1. This will be limited to 2.", ERR_MSG_SIZE);
    saturation = 2.0F;
  }
  if (vibrance > 2.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; vibrance should be in the range 0 - 2. Default 1. This will be limited to 2.", ERR_MSG_SIZE);
    vibrance = 2.0F;
  }
  if (whiteBalance > 1.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; whiteBalance should be in the range 0 - 1. Default 0. This will be limited to 1.", ERR_MSG_SIZE);
    whiteBalance = 1.0F;
  }
  if (uniformContrast > 1.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; uniformContrast should be in the range 0 - 1. Default 0. This will be limited to 1.", ERR_MSG_SIZE);
    uniformContrast = 1.0F;
  }
  if (normalize > 1.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; normalize should be in the range 0 - 1. Default 0. This will be limited to 1.", ERR_MSG_SIZE);
    normalize = 1.0F;
  }
  if (vignette > 1.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; vignette should be in the range 0 - 1. Default 0. This will be limited to 1.", ERR_MSG_SIZE);
    vignette = 1.0F;
  }
  if (chromAberr > 1.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; chromAberr should be in the range 0 - 1. This will be limited to 1.", ERR_MSG_SIZE);
    chromAberr = 1.0F;
  }
  if (sharpness > 1.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; sharpness should be in the range 0 - 1. This will be limited to 1.", ERR_MSG_SIZE);
    sharpness = 1.0F;
  }
  if (sizeStar > 100.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; sizeStar should be in the range 0 - 100. This will be limited to 100.", ERR_MSG_SIZE);
    sizeStar = 100.0F;
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
  if (randomAngle > 1.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; randomAngle should be in the range 0 - 1. This will be limited to 1.", ERR_MSG_SIZE);
    randomAngle = 1.0F;
  }
  if (sprayRay > 1.0F)
  {
    wcsncpy_s(m_msg, L"post_process_hydra1; sprayRay should be in the range 0 - 1. This will be limited to 1.", ERR_MSG_SIZE);
    sprayRay = 1.0F;
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

  ExecutePostProcessHydra1( input, output, numThreads, exposure, compress, contrast, saturation, vibrance, whiteBalance, whitePointColor,
    uniformContrast, normalize, vignette, chromAberr, sharpness, sizeStar, numRay, rotateRay, randomAngle, sprayRay, w1, h1);

  return true;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ExecutePostProcessHydra1(
  const float* a_input, float* a_output, unsigned int a_numThreads,
  const float a_exposure, const float a_compress, const float a_contrast, const float a_saturation, const float a_vibrance,
  const float a_whiteBalance, float3 a_whitePointColor, const float a_uniformContrast, const float a_normalize, const float a_vignette,
  const float a_chromAberr, const float a_sharpness, const float a_sizeStar, const int a_numRay, const int a_rotateRay,
  const float a_randomAngle, const float a_sprayRay, unsigned int a_width, unsigned int a_height)
{
  // Set the desired number of threads
  const unsigned int maxThreads = omp_get_num_procs();
  if (a_numThreads > maxThreads) a_numThreads = maxThreads;
  else if (a_numThreads <= 0)    a_numThreads = maxThreads;
  omp_set_num_threads(a_numThreads);

  float4* image4in          = (float4*)a_input;
  float4* image4out         = (float4*)a_output;

  // variables and constants
  const int   sizeImage     = a_width * a_height;
  const float centerImageX  = (float)a_width / 2.0F;
  const float centerImageY  = (float)a_height / 2.0F;
  const float diagonalImage = Distance(0.0F, 0.0F, (float)a_width, (float)a_height);
  const float radiusImage   = diagonalImage / 2.0F;

  const int histogramBin    = 10000;

  bool autoWhiteBalance = true;
  if (a_whitePointColor.x > 0.0F || a_whitePointColor.y > 0.0F || a_whitePointColor.z > 0.0F)
    autoWhiteBalance = false;

  float3 autoWhitePointRgb   = { 0.0F, 0.0F, 0.0F };
  float3 summRgbColorPass    = { 0.0F, 0.0F, 0.0F };
  float3 whitePointColorPass = { 1.0F, 1.0F, 1.0F };
  float3 d65                 = { 0.9505F , 1.0F, 1.0888F }; // in XYZ

  float  minRgb              = 1e6F;
  float  minRgb2             = 1e6F;
  float  maxRgb              = 0.0F;
  float  maxRgb2             = 0.0F;
  float  minAllHistBin       = 0.0F;
  float  maxAllHistBin       = 1.0F;
  float3 minHistBin(0.0F, 0.0F, 0.0F);
  float3 maxHistBin(histogramBin, histogramBin, histogramBin);

  // arrays
  std::vector<float> chromAbberR;
  std::vector<float> chromAbberG;
  std::vector<float> lumForSharp;
  std::vector<float3> diffrStars;
  std::vector<float3> histogram;
  //std::vector<float> blurPassR;
  //std::vector<float> blurPassG;
  //std::vector<float> blurPassB;


  if (a_chromAberr > 0.0F)
  {
    chromAbberR.resize(sizeImage);
    chromAbberG.resize(sizeImage);
  }
  if (a_normalize > 0.0F)
    histogram.resize(histogramBin);
  if (a_sharpness > 0.0F)
    lumForSharp.resize(sizeImage);
  if (a_sizeStar > 0.0F)
    diffrStars.resize(sizeImage);
  //if (a_compress > 0.0F || a_contrast > 1.0F || a_whiteBalance > 0.0F || a_diffuse > 0.0F)
  //{
  //  blurPassR.resize(sizeImage);
  //  blurPassG.resize(sizeImage);
  //  blurPassB.resize(sizeImage);
  //}

  //////////////////////////////////////////////////////////////////////////////////////


  // ----- Analization, precalculate and 3 effects. Loop 1. -----
  //#pragma omp parallel for 
  for (int y = 0; y < a_height; ++y)
  {
    for (int x = 0; x < a_width; ++x)
    {
      const int i = y * a_width + x;
      image4out[i] = image4in[i];
      float3 tempImgRGB = float3(image4out[i].x, image4out[i].y, image4out[i].z);


      // ----- blur pass -----
      //if (a_compress > 0.0F || a_contrast > 1.0F || a_whiteBalance > 0.0F || a_diffuse > 0.0F)
      //{
      //  blurPassR[i] = tempImgRGB.x;
      //  blurPassG[i] = tempImgRGB.y;
      //  blurPassB[i] = tempImgRGB.z;
      //}


      if (a_sharpness > 0.0F)
      {
        lumForSharp[i] = (tempImgRGB.x + tempImgRGB.y + tempImgRGB.z) / 3.0F;
        lumForSharp[i] /= (1.0F + lumForSharp[i]);
      }


      if (a_exposure != 1.0F)
        tempImgRGB *= a_exposure;


      ClampMinusToZero(tempImgRGB);
      MinMaxRgb(tempImgRGB, minRgb, maxRgb);


      if (a_whiteBalance > 0.0F && autoWhiteBalance)
        CalculateAutoWhitePoint(tempImgRGB, autoWhitePointRgb);


      if (a_vignette > 0.0F)
        Vignette(tempImgRGB, a_vignette, centerImageX, centerImageY, radiusImage, x, y);


      if (a_sizeStar > 0.0F && (image4out[i].x > 50.0F || image4out[i].y > 50.0F || image4out[i].z > 50.0F))
        DiffractionStars(tempImgRGB, diffrStars[0], a_sizeStar, a_numRay, a_rotateRay, a_randomAngle, a_sprayRay, a_width, a_height,
          sizeImage, radiusImage, x, y, i);


      if (a_chromAberr > 0.0F)
        ChrommAberr(tempImgRGB, &chromAbberR[0], &chromAbberG[0], a_width, a_height, sizeImage, a_chromAberr, x, y);


      image4out[i].x = tempImgRGB.x;
      image4out[i].y = tempImgRGB.y;
      image4out[i].z = tempImgRGB.z;
    }
  }// end loop 1.


  // ----- Chromatic aberration -----
  if (a_chromAberr > 0.0F)
  {
#pragma omp parallel for 
    for (int a = 0; a < sizeImage; ++a)
    {
      image4out[a].x = chromAbberR[a];
      image4out[a].y = chromAbberG[a];
    }
  }

  // ----- Sharpness -----
  if (a_sharpness > 0.0F)
    Sharp(image4out, &lumForSharp[0], a_sharpness, a_width, a_height);


  // ----- Blur pass-----
  //if (a_compress > 0.0F || a_contrast > 1.0F || a_whiteBalance > 0.0F)
  //{
  //  const int radius = (int)(radiusImage * 0.8F);
  //    Blur(&blurPassR[0], radius, a_width, a_height, sizeImage);
  //    Blur(&blurPassG[0], radius, a_width, a_height, sizeImage);
  //    Blur(&blurPassB[0], radius, a_width, a_height, sizeImage);
  //}

  // ----- White point for white balance -----
  if (a_whiteBalance > 0.0F)
    GetWhitePointForWhiteBalance(autoWhiteBalance, autoWhitePointRgb, a_whitePointColor, sizeImage, d65);


  // ----- Many filters. Loop 2. -----
#pragma omp parallel for
  for (int i = 0; i < sizeImage; ++i)
  {
    float3 tempImgRGB(image4out[i].x, image4out[i].y, image4out[i].z);

    // ----- Diffuse -----
    //if (a_diffuse > 0.0F)
    //{
    //  tempImgRGB.x = blurPassR[i];
    //  tempImgRGB.y = blurPassG[i];
    //  tempImgRGB.z = blurPassB[i];
    //}


    // ----- Diffraction stars -----
    if (a_sizeStar > 0.0F)
      tempImgRGB += diffrStars[i];


    // ----- White balance -----
    if (a_whiteBalance > 0.0F)
      WhiteBalance(tempImgRGB, a_whitePointColor, d65, a_whiteBalance);


    // ----- Saturation  -----
    if (a_saturation != 1.0F)
      Saturation(tempImgRGB, a_saturation);


    ////////// IPT block //////////

    if (a_vibrance != 1.0F || (a_compress > 0.0F && maxRgb > 1.0F) || a_contrast > 1.0F || a_uniformContrast > 0.0F)
    {
      float3 tempImgIPT = tempImgRGB;
      ConvertSrgbToXyz(tempImgIPT);
      ConvertXyzToLmsPower(tempImgIPT, 0.43F);
      ConvertLmsToIpt(tempImgIPT);


      // ----- Compress -----
      if ((a_compress > 0.0F && maxRgb > 1.0F) || a_uniformContrast > 0.0F)
        CompressIPT(tempImgIPT, fmax(a_compress, a_uniformContrast));


      // ----- Contrast -----        
      if (a_contrast > 1.0F)
        ContrastIPT(tempImgIPT, a_contrast);


      // ----- Vibrance -----
      if (a_vibrance != 1.0F)
        VibranceIPT(tempImgIPT, tempImgRGB, a_vibrance);


      ConvertIptToLms(tempImgIPT);
      ConvertLmsToXyzPower(tempImgIPT, 1.0F / 0.43F);
      ConvertXyzToSrgb(tempImgIPT);

      ////////// end IPT block //////////

      tempImgRGB = tempImgIPT;

      // little compress in RGB
      if ((a_compress > 0.0F && maxRgb > 1.0F) || a_contrast > 1.0F || a_uniformContrast > 0.0F)
      {
        CompressWithKnee(tempImgRGB.x, 0.01F);
        CompressWithKnee(tempImgRGB.y, 0.01F);
        CompressWithKnee(tempImgRGB.z, 0.01F);
      }

      ClampMinusToZero(tempImgRGB);
    }

    MinMaxRgb(tempImgRGB, minRgb2, maxRgb2);

    image4out[i].x = tempImgRGB.x;
    image4out[i].y = tempImgRGB.y;
    image4out[i].z = tempImgRGB.z;
  } // end loop 2.


  if (maxRgb2 > 1e-3F) maxRgb = fmin(maxRgb, maxRgb2);


  // ----- Uniform contrast. -----
  if (a_uniformContrast > 0.0F && maxRgb > 0.01F)
    UniformContrastIPTFilter(&image4out[0], sizeImage, histogramBin, a_uniformContrast);


  // ----- Normalize. -----
  if (a_normalize > 0.0F)
    NormalizeFilter(&image4out[0], &histogram[0], sizeImage, histogramBin, minHistBin, maxHistBin, minAllHistBin, maxAllHistBin, a_normalize);
}



