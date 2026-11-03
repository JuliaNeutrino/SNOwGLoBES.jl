/* Event rates of an AEDL file with C GLoBES and no oscillations: for each
 * experiment and channel, the events before and after the detector response.
 * Build against GLoBES 3.0.11 with GLoBES.jl's test/reference/c-fixes.patch. */
#include <stdio.h>
#include <globes/globes.h>

int main(int argc, char **argv)
{
  glbInit(argv[0]);
  for (int k = 1; k < argc; k++)
    glbInitExperiment(argv[k], &glb_experiment_list[0], &glb_num_of_exps);
  glb_params p = glbAllocParams();
  glbDefineParams(p, 0, 0, 0, 0, 0, 0);
  glbSetDensityParams(p, 1.0, GLB_ALL);
  glbSetOscillationParameters(p);
  glbSetRates();
  for (int e = 0; e < glb_num_of_exps; e++)
    for (int c = 0; c < glbGetNumberOfChannels(e); c++) {
      double *pre = glbGetChannelRatePtr(e, c, GLB_PRE), *post = glbGetChannelRatePtr(e, c, GLB_POST);
      double a = 0, b = 0;
      for (int i = 0; i < glbGetNumberOfSamplingPoints(e); i++) a += pre[i];
      for (int i = 0; i < glbGetNumberOfBins(e); i++) b += post[i];
      printf("%d %d %.10e %.10e\n", e, c, a, b);
    }
  glbFreeParams(p);
  return 0;
}
