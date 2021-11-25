#include "file_ops.h"
#include "stddef.h"

static double *data_double = NULL;

void set_data(double *data_in){
  data_double = data_in;
}
