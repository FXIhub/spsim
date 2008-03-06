#include <hdf5.h>
#include <spimage.h>

int main(int argc, char ** argv){
  Image * img;
  if(argc < 1){
    fprintf(stderr,"Usage: %s <file>\n",argv[0]);
  }
  img = sp_image_read(argv[1],0);
  real sum = 0;
  for(int i = 0; i<sp_image_size(img);i++){
    sum += sp_cabs(img->image->data[i]);
  }
  printf("Integral - %f\n",sum);
  return 0;
}


