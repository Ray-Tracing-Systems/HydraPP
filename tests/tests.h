#ifndef HYDRA_PP_TESTS_H
#define HYDRA_PP_TESTS_H

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

#endif //HYDRA_PP_TESTS_H
