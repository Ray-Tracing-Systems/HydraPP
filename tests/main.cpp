namespace PP_TESTS
{
  bool test301_resample();
  bool test302_median();
  bool test303_median_in_place();
  bool test304_obsolete_tone_mapping();
  bool test305_fbi_from_render();

  bool test306_post_process_hydra1_exposure05();
  bool test307_post_process_hydra1_exposure2();
  bool test308_post_process_hydra1_compress();
  bool test309_post_process_hydra1_contrast();
  bool test310_post_process_hydra1_desaturation();
  bool test311_post_process_hydra1_saturation();
  bool test312_post_process_hydra1_whiteBalance();
  bool test312_2_post_process_hydra1_whitePointColor();
  bool test313_post_process_hydra1_uniformContrast();
  bool test314_post_process_hydra1_normalize();
  bool test315_post_process_hydra1_vignette();
  bool test316_post_process_hydra1_chromAberr();
  bool test317_post_process_hydra1_sharpness();
  bool test318_post_process_hydra1_ECCSWUNSVC();
  bool test319_post_process_hydra1_diffStars();

  bool test320_blur();
  bool test321_median_mostly_bad_pixels();

};

void run_all_ipp_tests(int a_start) 
{
  TestFunc tests[] = {                 
    &PP_TESTS::test301_resample, 
    &PP_TESTS::test302_median,
    &PP_TESTS::test303_median_in_place,
    &PP_TESTS::test304_obsolete_tone_mapping,
    &PP_TESTS::test305_fbi_from_render,
    &PP_TESTS::test306_post_process_hydra1_exposure05,
    &PP_TESTS::test307_post_process_hydra1_exposure2,
    &PP_TESTS::test308_post_process_hydra1_compress,
    &PP_TESTS::test309_post_process_hydra1_contrast,
    &PP_TESTS::test310_post_process_hydra1_desaturation,
    &PP_TESTS::test311_post_process_hydra1_saturation,
    &PP_TESTS::test312_post_process_hydra1_whiteBalance,
    &PP_TESTS::test312_2_post_process_hydra1_whitePointColor,
    &PP_TESTS::test313_post_process_hydra1_uniformContrast,
    &PP_TESTS::test314_post_process_hydra1_normalize,
    &PP_TESTS::test315_post_process_hydra1_vignette,
    &PP_TESTS::test316_post_process_hydra1_chromAberr,
    &PP_TESTS::test317_post_process_hydra1_sharpness,
    &PP_TESTS::test318_post_process_hydra1_ECCSWUNSVC,
    &PP_TESTS::test319_post_process_hydra1_diffStars
  };

  std::ofstream fout("z_test_post_process.txt");

  const int testNum = sizeof(tests) / sizeof(TestFunc);

  for (int i = a_start; i < testNum; i++)
  {
    bool res = tests[i]();
    if (res)
    {
      std::cout          << "pp_test_" << std::setfill('0') << std::setw(3) << 300 + i << "\tPASSED!\t\n";
      fout << std::fixed << "pp_test_" << std::setfill('0') << std::setw(3) << 300 + i << "\tPASSED!\t\n";
    }
    else if (g_testWasIgnored)
    {
      std::cout          << "pp_test_" << std::setfill('0') << std::setw(3) << 300 + i << "\tSKIPPED!\t\n";
      fout << std::fixed << "pp_test_" << std::setfill('0') << std::setw(3) << 300 + i << "\tSKIPPED!\t\n";
    }
    else
    {
      std::cout          << "pp_test_" << std::setfill('0') << std::setw(3) << 300 + i << "\tFAILED! :-: MSE = " << g_MSEOutput << std::endl;
      fout << std::fixed << "pp_test_" << std::setfill('0') << std::setw(3) << 300 + i << "\tFAILED! :-: MSE = " << g_MSEOutput << std::endl;
    }

    fout.flush();

    g_testWasIgnored = false;
  }

  fout.close();
}
