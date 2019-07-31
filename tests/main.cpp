#include <fstream>
#include <iomanip>
#include <iostream>
#include <csignal>
#include "tests.h"

#include "../../HydraAPI/hydra_api/HydraAPI.h"
#include "../hydra_pp/HydraPostProcessAPI.h"

#ifndef WIN32
  #include <zconf.h>
#else
  #include <windows.h>
#endif

extern float g_MSEOutput;

void run_all_ipp_tests(int a_start = 0)
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

  _hrInitPostProcess();

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

void InfoCallBack(const wchar_t* message, const wchar_t* callerPlace, HR_SEVERITY_LEVEL a_level)
{
  if (a_level >= HR_SEVERITY_WARNING)
  {
    if (a_level == HR_SEVERITY_WARNING)
      std::wcerr << L"WARNING: " << callerPlace << L": " << message; // << std::endl;
    else
      std::wcerr << L"ERROR  : " << callerPlace << L": " << message; // << std::endl;
  }
}

void destroy()
{
  std::cout << "call destroy() --> hrSceneLibraryClose()" << std::endl;
  hrSceneLibraryClose();
}

#ifdef WIN32
BOOL WINAPI HandlerExit(_In_ DWORD fdwControl)
{
  exit(0);
  return TRUE;
}
#else
bool destroyedBySig = false;
void sig_handler(int signo)
{
  if(destroyedBySig)
    return;
  switch(signo)
  {
    case SIGINT : std::cerr << "\nmain_app, SIGINT";      break;
    case SIGABRT: std::cerr << "\nmain_app, SIGABRT";     break;
    case SIGILL : std::cerr << "\nmain_app, SIGINT";      break;
    case SIGTERM: std::cerr << "\nmain_app, SIGILL";      break;
    case SIGSEGV: std::cerr << "\nmain_app, SIGSEGV";     break;
    case SIGFPE : std::cerr << "\nmain_app, SIGFPE";      break;
    default     : std::cerr << "\nmain_app, SIG_UNKNOWN"; break;
      break;
  }
  std::cerr << " --> hrSceneLibraryClose()" << std::endl;
  hrSceneLibraryClose();
  destroyedBySig = true;
}
#endif


int main(int argc, char* argv[])
{
  hrInfoCallback(&InfoCallBack);
  hrErrorCallerPlace(L"main");  // for debug needs only

  std::string workingDir = "../../tests";
  if(argc > 1)
    workingDir = std::string(argv[1]);

  atexit(&destroy);                           // if application will terminated you have to call hrSceneLibraryClose to free all connections with hydra.exe
#if defined WIN32
  SetConsoleCtrlHandler(&HandlerExit, TRUE);  // if some one kill console :)
  wchar_t NPath[512];
  GetCurrentDirectoryW(512, NPath);
#ifdef NEED_DIR_CHANGE
  SetCurrentDirectoryW(L"../../main");
#endif
  std::wcout << L"[main]: curr_dir = " << NPath << std::endl;
#else

  chdir(workingDir.c_str());
  char cwd[1024];
  if (getcwd(cwd, sizeof(cwd)) != nullptr)
    std::cout << "[main]: curr_dir = " << cwd <<std::endl;
  else
    std::cout << "getcwd() error" <<std::endl;

  {
    struct sigaction sigIntHandler;
    sigIntHandler.sa_handler = sig_handler;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = SA_RESETHAND;
    sigaction(SIGINT,  &sigIntHandler, NULL);
    sigaction(SIGSTOP, &sigIntHandler, NULL);
    sigaction(SIGABRT, &sigIntHandler, NULL);
    sigaction(SIGILL,  &sigIntHandler, NULL);
    sigaction(SIGTERM, &sigIntHandler, NULL);
    sigaction(SIGSEGV, &sigIntHandler, NULL);
    sigaction(SIGFPE,  &sigIntHandler, NULL);
  }
#endif

  std::cout << "sizeof(size_t) = " << sizeof(size_t) <<std::endl;

  try
  {
    run_all_ipp_tests();
  }
  catch (std::runtime_error& e)
  {
    std::cout << "std::runtime_error: " << e.what() << std::endl;
  }
  catch (...)
  {
    std::cout << "unknown exception" << std::endl;
  }

  hrErrorCallerPlace(L"main"); // for debug needs only

  hrSceneLibraryClose();

  return 0;
}
