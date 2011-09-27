static char rcsid[] = "$Id: littleendian.c,v 1.1 2010-01-12 22:51:55 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "littleendian.h"
#include <unistd.h>

void
Littleendian_write_uint (unsigned int value, int fd) {
  char buf[4];

  buf[0] = (char) (value & 0xff);
  buf[1] = (char) ((value >>= 8) & 0xff);
  buf[2] = (char) ((value >>= 8) & 0xff);
  buf[3] = (char) ((value >>= 8) & 0xff);
  write(fd,buf,4);
}

