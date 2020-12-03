#pragma once

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


//------------------------------------------------------------------------

void ExecutePostProcessHydra1(
  const float* a_input, float* a_output, unsigned int a_numThreads,
  const float a_exposure, const float a_compress, const float a_contrast, const float a_saturation, const float a_vibrance,
  const float a_whiteBalance, float3 a_whitePointColor, const float a_uniformContrast, const float a_normalize, const float a_vignette,
  const float a_chromAberr, const float a_sharpness, const float a_sizeStar, const int a_numRay, const int a_rotateRay,
  const float a_randomAngle, const float a_sprayRay, unsigned int a_width, unsigned int a_height); 


float Luminance(const float3 a_data);
float Mean(const float3 a_data);
void SpectralToRGB(float &r, float &g, float &b, const float l); // RGB <0,1> <- lambda l <400,700> [nm]
void ConvertRgbToHsv(float4* data);
void ConvertHsvToRgb(float4* data);
void ConvertRgbToXyz65(float4* data);
void ConvertSrgbToXyz(float3& a_data);
void ConvertSrgbToXyzIcam06(float4* data);
void ConvertXyz65ToRgb(float3& a_data);
void ConvertXyzToSrgb(float3& a_data);
void ConvertXyzToSrgbIcam06(float3& a_data);
void ConvertXyzToLab(float3& a_data);
void ConvertXyzToZlab(float4* data, const float3 whitePoint, const float a_whiteBalance);
void ConvertXyzToRlab(float4* data, float power);
void ConvertRlabToXyz(float4* data);
void ConvertLabToXyz(float4* data);
void ConvertZlabToXyz(float4* data, const float3 whitePoint);
void MatrixCat02(float3& a_data);
void InverseMatrixCat02(float3& a_data);
void MatrixHpe(float3& a_data);
void InverseMatrixHpe(float3& a_data);
void ConvertXyzToLmsVonKries(float3& a_data);
void ConvertLmsToXyzVonKries(float3& a_data);
void ConvertXyzToLms(float3& a_data);
void ConvertXyzToLmsPower(float3& a_data, const float a_power);
void ExponentIcam(float3& a_data, const float a_dataYblur);
void ConvertLmsToIpt(float3& a_data);
void ConvertLmsToXyz(float3& a_data);
void ConvertLmsToXyzPower(float3& a_data, const float a_power);
void ConvertIptToLms(float3& a_data);
void ChromAdaptIcam06(float3& a_data, const float3 a_dataW, const float3 a_d65, const float a_D);
void InverseChromAdaptIcam(float3& a_data, const float D);
void ChromAdaptCam02(float3& a_data, const float3 a_dataW, const float a_D);
void InverseChromAdaptCam02(float3& a_data, const float3 a_dataW, const float a_Fl, const float a_D, const float a_c);
void GetWhitePointForWhiteBalance(const bool a_autoWhiteBalance, const float3 a_summRgb, float3& a_whitePointColor, const int a_sizeImage, float3& a_d65);
void WhiteBalance(float3& a_data, const float3 a_whitePointColor, const float3 a_d65, const float a_whiteBalance);
void Saturation(float3& a_data, const float a_saturation);
void VibranceIPT(float3 & a_dataIPT, const float3 a_dataRGB, const float a_saturation);  // saturation of unsaturated colors. 
void CompressIPT(float3 & a_dataIPT, const float a_compress);
void ContrastIPT(float3& a_dataIPT, const float a_contrast);
void ChangeColorWithNewLuminance_IPT(const float & a_newLum, float3& a_dataIPT);
void MoreCompressColor_IPT(float3& a_dataIPT);
void ContrastRGB(float3& a_data, const float a_contrast);


template <typename T>
void Normalize(T& data, const float inMin, const float inMax, const float outMin, const float outMax);
float ContrastField(const float a_data, const float a_amount);
float Distance(const float a, const float b);
float Distance(const float x1, const float y1, const float x2, const float y2);
float Max(const float value1, const float value2, const float value3);
float Min(const float value1, const float value2, const float value3);
float3 Mediana(const float4* inData, const int size, const int a_width, const int a_height, const int x, const int y);
int FloatToInt(const float inData);
void AssignHistogramToData(float & a_data, const float * a_histogram, const int a_histogramBin);
void Bilinear(const float inData, float outData[], float newPosX, float newPosY, const int m_width, const int m_height, const int sizeImage);
void Bilinear3(const float3 inData, float3* outData, float newPosX, float newPosY, const int m_width, const int m_height, const int sizeImage);
// 0 - data1, 1 - data2
void Blend(float& inData1, const float inData2, const float amount);
void Blend(float3& inData1, const float inData2, const float amount);
void Blur(float* a_data, int a_blurRadius, const int a_width, const int a_height, const int a_sizeImage);
void Blur1D(std::vector<float>& a_array, const int a_sizeArray, const int a_blurRadius, std::vector<float>& a_weights);
void CalculateHistogram(const float & a_data, float* a_histogram, const int a_histogramBin, const float a_iterrHistBin);
void ChrommAberr(const float3 inData, float outData1[], float outData2[], const int m_width, const int m_height, const int sizeImage, const float a_chromAberr, const int x, const int y);
void ClampMinusToZero(float3& a_data);
void CompressWithMax(float& data, const float maxRgb);
void CompressWithKnee(float& a_data, const float a_compress);
void ComputeWhitePoint(const float3 a_summRgb, float3& a_whitePoint, const int a_sizeImage);
std::vector<float> CreateGaussKernelWeights1D_HDRImage2(const int blurRadius);
void CummulativeHistogram(float* a_histogram, const int a_histogramBin);
void DrawColumn(float4 data[], const int width, const int height, const int offsetX, const int ImgWidth, const float3 color, const bool transparent);
void FireflyDetects(const float4* image4out, float3* fireflyPixels, const int a_width, const int a_height, const int x, const int y, const int i, const float a_fireflyRemove);
void Mediana(float inData[], const int radius, const int a_width, const int a_height, const int sizeImage);
void MinMaxHistBin(const float* histogram, float& minHistBin, float& maxHistBin, const int sizeImage, const int histogramBin);
void MinMaxRgb(const float3 a_data, float& a_minRgb, float& a_maxRgb);
void NormalizeHistogram(float* a_histogram, const int a_histogramBin);
void OffsetCenter(float& data, const float inMin, const float inMax, const float coef);
void Resize(float* a_data, const int a_sourceWidth, const int a_sourceHeight, const int a_sourceSizeImage, const float a_resize);
void Sharp(float4 image4out[], const float lumForSharp[], const float a_sharpness, const int a_width, const int a_height);
float WeightWhiteBalance(const float a_lum);
void CalculateAutoWhitePoint(const float3 a_data, float3& a_autoWhitePoint);
void UniformContrast(float* a_data, const int a_sizeImage, const int a_histogramBin);
void UniformContrastRgb(float* a_data, const int a_sizeImage, const int a_histogramBin);
void ViewHistorgam(float4 data[], float3* histogram, const int m_width, const int m_height, const int histogramBin);
void Vignette(float3& a_data, const float a_vignette, const float a_centerImageX, const float a_centerImageY, const float a_radiusImage, const int a_x, const int a_y);
void DiffractionStars(const float3 a_data, float3& a_diffrStars, const float a_sizeStar, const int a_numRay, const int a_RotateRay,
  const float a_randomAngle, const float a_sprayRay, const int m_width, const int m_height, const int sizeImage, const float radiusImage, const int x, const int y, const int i);
void UniformContrastRGBFilter(float4* a_data, const int a_sizeImage, const int a_histogramBin, const float a_amount);
void UniformContrastIPTFilter(float4* a_data, const int a_sizeImage, const int a_histogramBin, const float a_amount);
void NormalizeFilter(float4* a_data, float3* a_histogram, const int sizeImage, const int a_histogramBin, float3& minHistBin, 
  float3& maxHistBin, float& minAllHistBin, float& maxAllHistBin, const float amount);
void ViewHistorgamFilter(float4* a_data, std::vector<float3>& a_histogram, const int a_width, const int a_height, const bool a_viewHistogram);
void CalculateStatisticsImage(const float4* a_data, float& a_minRgb, float& a_maxRgb, std::vector<float3>& a_histogram, const int a_width, const int a_height);
