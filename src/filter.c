#include <fftw3.h>

float * filter_detector_output(CCD * det,float obj_size){
  fftw_complex * in, * out;
  in = fftw_malloc(sizeof(fftw_complex) * det->nx*det->ny);
  out = fftw_malloc(sizeof(fftw_complex) * det->nx*det->ny);
  for(i =0 ;i<det->nx*det->ny;i++){
    in[i][0] = det->real_output[i];
  }
  fftw_plan_dft_2d(det->nx,det->ny,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(p);
  /* zero out high resolution */
  
}
