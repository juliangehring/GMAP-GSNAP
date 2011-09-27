static char rcsid[] = "$Id: bigendian.c,v 1.7 2005/02/08 00:00:51 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bigendian.h"


/************************************************************************
 *   Int
 ************************************************************************/

int
Bigendian_convert_int (int littleendian) {
  int bigendian;
  char byte1, byte2, byte3;

  bigendian = littleendian & 0xff;
  byte1 = (littleendian >>= 8);
  byte2 = (littleendian >>= 8);
  byte3 = (littleendian >>= 8);

  /* bigendian = byte0; */
  bigendian <<= 8;
  bigendian |= (byte1 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte2 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte3 & 0xff);

  return bigendian;
}

size_t
Bigendian_fwrite_int (int value, FILE *fp) {
  char buf[4];

  buf[0] = value & 0xff;
  buf[1] = (value >>= 8) & 0xff;
  buf[2] = (value >>= 8) & 0xff;
  buf[3] = (value >>= 8) & 0xff;
  if (fwrite(buf,sizeof(char),4,fp) == 0) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    return 1;
  }
}

size_t
Bigendian_fwrite_ints (int *array, int n, FILE *fp) {
  size_t nitems = 0;
  char buf[4];
  int value, i;

  for (i = 0; i < n; i++) {
    value = array[i];
    buf[0] = value & 0xff;
    buf[1] = (value >>= 8) & 0xff;
    buf[2] = (value >>= 8) & 0xff;
    buf[3] = (value >>= 8) & 0xff;
    if (fwrite(buf,sizeof(char),4,fp) == 0) {
      /* Should set error indicator for stream and set errno */
      return 0;
    }
  }
  return n;
}

size_t
Bigendian_fread_int (int *value, FILE *fp) {
  char buf[4];

  if (fread(buf,sizeof(char),4,fp) < 4) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    *value = (buf[3] & 0xff);
    *value <<= 8;
    *value |= (buf[2] & 0xff);
    *value <<= 8;
    *value |= (buf[1] & 0xff);
    *value <<= 8;
    *value |= (buf[0] & 0xff);
    return 1;
  }
}

size_t
Bigendian_fread_ints (int *array, int n, FILE *fp) {
  char buf[4];
  int value, i;

  for (i = 0; i < n; i++) {
    if (fread(buf,sizeof(char),4,fp) < 4) {
      /* Should set error indicator for stream and set errno */
      return 0;
    } else {
      value = (buf[3] & 0xff);
      value <<= 8;
      value |= (buf[2] & 0xff);
      value <<= 8;
      value |= (buf[1] & 0xff);
      value <<= 8;
      value |= (buf[0] & 0xff);
      array[i] = value;
    }
  }
  return n;
}


/************************************************************************
 *   Unsigned int
 ************************************************************************/

unsigned int
Bigendian_convert_uint (unsigned int littleendian) {
  unsigned int bigendian;
  char byte1, byte2, byte3;

  bigendian = littleendian & 0xff;
  byte1 = (littleendian >>= 8);
  byte2 = (littleendian >>= 8);
  byte3 = (littleendian >>= 8);

  /* bigendian = byte0; */
  bigendian <<= 8;
  bigendian |= (byte1 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte2 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte3 & 0xff);

  return bigendian;
}


size_t
Bigendian_fwrite_uint (unsigned int value, FILE *fp) {
  char buf[4];

  buf[0] = value & 0xff;
  buf[1] = (value >>= 8) & 0xff;
  buf[2] = (value >>= 8) & 0xff;
  buf[3] = (value >>= 8) & 0xff;
  if (fwrite(buf,sizeof(char),4,fp) == 0) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    return 1;
  }
}


size_t
Bigendian_fwrite_uints (unsigned int *array, int n, FILE *fp) {
  char buf[4];
  unsigned int value;
  int i;
  
  for (i = 0; i < n; i++) {
    value = array[i];
    buf[0] = value & 0xff;
    buf[1] = (value >>= 8) & 0xff;
    buf[2] = (value >>= 8) & 0xff;
    buf[3] = (value >>= 8) & 0xff;
    if (fwrite(buf,sizeof(char),4,fp) == 0) {
      /* Should set error indicator for stream and set errno */
      return 0;
    }
  }
  return n;
}

size_t
Bigendian_fread_uint (unsigned int *value, FILE *fp) {
  char buf[4];

  if (fread(buf,sizeof(char),4,fp) < 4) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    *value = (buf[3] & 0xff);
    *value <<= 8;
    *value |= (buf[2] & 0xff);
    *value <<= 8;
    *value |= (buf[1] & 0xff);
    *value <<= 8;
    *value |= (buf[0] & 0xff);
    return 1;
  }
}

size_t
Bigendian_fread_uints (unsigned int *array, int n, FILE *fp) {
  char buf[4];
  unsigned int value;
  int i;

  for (i = 0; i < n; i++) {
    if (fread(buf,sizeof(char),4,fp) < 4) {
      /* Should set error indicator for stream and set errno */
      return 0;
    } else {
      value = (buf[3] & 0xff);
      value <<= 8;
      value |= (buf[2] & 0xff);
      value <<= 8;
      value |= (buf[1] & 0xff);
      value <<= 8;
      value |= (buf[0] & 0xff);
      array[i] = value;
    }
  }
  return n;
}
