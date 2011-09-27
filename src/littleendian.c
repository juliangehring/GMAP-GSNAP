static char rcsid[] = "$Id: littleendian.c 40271 2011-05-28 02:29:18Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "littleendian.h"
#include <unistd.h>

void
Littleendian_write_uint (unsigned int value, int fd) {
  unsigned char buf[4];

  buf[0] = (unsigned char) (value & 0xff);
  buf[1] = (unsigned char) ((value >>= 8) & 0xff);
  buf[2] = (unsigned char) ((value >>= 8) & 0xff);
  buf[3] = (unsigned char) ((value >>= 8) & 0xff);
  write(fd,buf,4);
}

