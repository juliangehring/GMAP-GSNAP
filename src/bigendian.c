static char rcsid[] = "$Id: bigendian.c 99737 2013-06-27 19:33:03Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bigendian.h"
#include <unistd.h>		/* For read() */


/*************************************************************************
 *  OUTPUT_BIGENDIAN provided to test bigendian code on a littleendian
 *  machine.  To use, compile all programs with WORDS_BIGENDIAN defined
 *  in config.h and define OUTPUT_BIGENDIAN here.
 ************************************************************************/

/************************************************************************
 *   Int
 ************************************************************************/

#ifdef OUTPUT_BIGENDIAN
int
Bigendian_convert_int (int littleendian) {
  return littleendian;
}
#else
int
Bigendian_convert_int (int littleendian) {
  int bigendian;

  bigendian = littleendian & 0xff; /* 0 */
  bigendian <<= 8;
  bigendian |= ((littleendian >>= 8) & 0xff); /* 1 */
  bigendian <<= 8;
  bigendian |= ((littleendian >>= 8) & 0xff); /* 2 */
  bigendian <<= 8;
  bigendian |= ((littleendian >>= 8) & 0xff); /* 3 */

  return bigendian;
}
#endif


#ifdef OUTPUT_BIGENDIAN
size_t
Bigendian_fwrite_int (int value, FILE *fp) {
  unsigned char buf[4];

  buf[3] = value & 0xff;
  buf[2] = (value >>= 8) & 0xff;
  buf[1] = (value >>= 8) & 0xff;
  buf[0] = (value >>= 8) & 0xff;
  if (fwrite(buf,sizeof(unsigned char),4,fp) == 0) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    return 1;
  }
}
#else
size_t
Bigendian_fwrite_int (int value, FILE *fp) {
  unsigned char buf[4];

  buf[0] = (unsigned char) (value & 0xff);
  buf[1] = (unsigned char) ((value >>= 8) & 0xff);
  buf[2] = (unsigned char) ((value >>= 8) & 0xff);
  buf[3] = (unsigned char) ((value >>= 8) & 0xff);
  if (fwrite(buf,sizeof(unsigned char),4,fp) == 0) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    return 1;
  }
}
#endif

#ifdef OUTPUT_BIGENDIAN
size_t
Bigendian_fwrite_ints (int *array, int n, FILE *fp) {
  unsigned char buf[4];
  int value, i;

  for (i = 0; i < n; i++) {
    value = array[i];
    buf[3] = value & 0xff;
    buf[2] = (value >>= 8) & 0xff;
    buf[1] = (value >>= 8) & 0xff;
    buf[0] = (value >>= 8) & 0xff;
    if (fwrite(buf,sizeof(unsigned char),4,fp) == 0) {
      /* Should set error indicator for stream and set errno */
      return 0;
    }
  }
  return n;
}
#else
size_t
Bigendian_fwrite_ints (int *array, int n, FILE *fp) {
  unsigned char buf[4];
  int value, i;

  for (i = 0; i < n; i++) {
    value = array[i];
    buf[0] = (unsigned char) (value & 0xff);
    buf[1] = (unsigned char) ((value >>= 8) & 0xff);
    buf[2] = (unsigned char) ((value >>= 8) & 0xff);
    buf[3] = (unsigned char) ((value >>= 8) & 0xff);
    if (fwrite(buf,sizeof(unsigned char),4,fp) == 0) {
      /* Should set error indicator for stream and set errno */
      return 0;
    }
  }
  return n;
}
#endif

#ifdef OUTPUT_BIGENDIAN
size_t
Bigendian_fread_int (int *value, FILE *fp) {
  unsigned char buf[4];

  if (fread(buf,sizeof(unsigned char),4,fp) < 4) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    *value = (buf[0] & 0xff);
    *value <<= 8;
    *value |= (buf[1] & 0xff);
    *value <<= 8;
    *value |= (buf[2] & 0xff);
    *value <<= 8;
    *value |= (buf[3] & 0xff);
    return 1;
  }
}
#else
size_t
Bigendian_fread_int (int *value, FILE *fp) {
  unsigned char buf[4];

  if (fread(buf,sizeof(unsigned char),4,fp) < 4) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
#if 0
    fprintf(stderr,"Reading %2X %2X %2X %2X, and using last as most sig\n",buf[0],buf[1],buf[2],buf[3]);
#endif
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
#endif

#ifdef OUTPUT_BIGENDIAN
size_t
Bigendian_fread_ints (int *array, int n, FILE *fp) {
  unsigned char buf[4];
  int value, i;

  for (i = 0; i < n; i++) {
    if (fread(buf,sizeof(unsigned char),4,fp) < 4) {
      /* Should set error indicator for stream and set errno */
      return 0;
    } else {
      value = (buf[0] & 0xff);
      value <<= 8;
      value |= (buf[1] & 0xff);
      value <<= 8;
      value |= (buf[2] & 0xff);
      value <<= 8;
      value |= (buf[3] & 0xff);
      array[i] = value;
    }
  }
  return n;
}
#else
size_t
Bigendian_fread_ints (int *array, int n, FILE *fp) {
  unsigned char buf[4];
  int value, i;

  for (i = 0; i < n; i++) {
    if (fread(buf,sizeof(unsigned char),4,fp) < 4) {
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
#endif


/************************************************************************
 *   Unsigned int
 ************************************************************************/

unsigned int
Bigendian_convert_uint (unsigned int littleendian) {
  unsigned int bigendian;

  bigendian = littleendian & 0xff; /* 0 */
  bigendian <<= 8;
  bigendian |= ((littleendian >>= 8) & 0xff); /* 1 */
  bigendian <<= 8;
  bigendian |= ((littleendian >>= 8) & 0xff); /* 2 */
  bigendian <<= 8;
  bigendian |= ((littleendian >>= 8) & 0xff); /* 3 */

  return bigendian;
}


#ifdef OUTPUT_BIGENDIAN
size_t
Bigendian_fwrite_uint (unsigned int value, FILE *fp) {
  unsigned char buf[4];

  buf[3] = value & 0xff;
  buf[2] = (value >>= 8) & 0xff;
  buf[1] = (value >>= 8) & 0xff;
  buf[0] = (value >>= 8) & 0xff;
  if (fwrite(buf,sizeof(unsigned char),4,fp) == 0) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    return 1;
  }
}
#else
size_t
Bigendian_fwrite_uint (unsigned int value, FILE *fp) {
  unsigned char buf[4];

  buf[0] = (unsigned char) (value & 0xff);
  buf[1] = (unsigned char) ((value >>= 8) & 0xff);
  buf[2] = (unsigned char) ((value >>= 8) & 0xff);
  buf[3] = (unsigned char) ((value >>= 8) & 0xff);
  if (fwrite(buf,sizeof(unsigned char),4,fp) == 0) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    return 1;
  }
}
#endif


#ifdef OUTPUT_BIGENDIAN
void
Bigendian_write_uint (unsigned int value, int fd) {
  unsigned char buf[4];

  buf[3] = value & 0xff;
  buf[2] = (value >>= 8) & 0xff;
  buf[1] = (value >>= 8) & 0xff;
  buf[0] = (value >>= 8) & 0xff;
  write(fd,buf,4);
  return;
}
#else
void
Bigendian_write_uint (unsigned int value, int fd) {
  unsigned char buf[4];

  buf[0] = (unsigned char) (value & 0xff);
  buf[1] = (unsigned char) ((value >>= 8) & 0xff);
  buf[2] = (unsigned char) ((value >>= 8) & 0xff);
  buf[3] = (unsigned char) ((value >>= 8) & 0xff);
  write(fd,buf,4);
}
#endif

#ifdef OUTPUT_BIGENDIAN
size_t
Bigendian_fwrite_uints (unsigned int *array, int n, FILE *fp) {
  unsigned char buf[4];
  unsigned int value;
  int i;
  
  for (i = 0; i < n; i++) {
    value = array[i];
    buf[3] = value & 0xff;
    buf[2] = (value >>= 8) & 0xff;
    buf[1] = (value >>= 8) & 0xff;
    buf[0] = (value >>= 8) & 0xff;
    if (fwrite(buf,sizeof(unsigned char),4,fp) == 0) {
      /* Should set error indicator for stream and set errno */
      return 0;
    }
  }
  return n;
}
#else
size_t
Bigendian_fwrite_uints (unsigned int *array, int n, FILE *fp) {
  unsigned char buf[4];
  unsigned int value;
  int i;
  
  for (i = 0; i < n; i++) {
    value = array[i];
    buf[0] = (unsigned char) (value & 0xff);
    buf[1] = (unsigned char) ((value >>= 8) & 0xff);
    buf[2] = (unsigned char) ((value >>= 8) & 0xff);
    buf[3] = (unsigned char) ((value >>= 8) & 0xff);
    if (fwrite(buf,sizeof(unsigned char),4,fp) == 0) {
      /* Should set error indicator for stream and set errno */
      return 0;
    }
  }
  return n;
}
#endif


#ifdef OUTPUT_BIGENDIAN
size_t
Bigendian_fread_uint (unsigned int *value, FILE *fp) {
  unsigned char buf[4];

  if (fread(buf,sizeof(unsigned char),4,fp) < 4) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    *value = (buf[0] & 0xff);
    *value <<= 8;
    *value |= (buf[1] & 0xff);
    *value <<= 8;
    *value |= (buf[2] & 0xff);
    *value <<= 8;
    *value |= (buf[3] & 0xff);
    return 1;
  }
}
#else
size_t
Bigendian_fread_uint (unsigned int *value, FILE *fp) {
  unsigned char buf[4];

  if (fread(buf,sizeof(unsigned char),4,fp) < 4) {
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
#endif

#ifdef OUTPUT_BIGENDIAN
size_t
Bigendian_fread_uints (unsigned int *array, int n, FILE *fp) {
  unsigned char buf[4];
  unsigned int value;
  int i;

  for (i = 0; i < n; i++) {
    if (fread(buf,sizeof(unsigned char),4,fp) < 4) {
      /* Should set error indicator for stream and set errno */
      return 0;
    } else {
      value = (buf[0] & 0xff);
      value <<= 8;
      value |= (buf[1] & 0xff);
      value <<= 8;
      value |= (buf[2] & 0xff);
      value <<= 8;
      value |= (buf[3] & 0xff);
      array[i] = value;
    }
  }
  return n;
}
#else
size_t
Bigendian_fread_uints (unsigned int *array, int n, FILE *fp) {
  unsigned char buf[4];
  unsigned int value;
  int i;

  for (i = 0; i < n; i++) {
    if (fread(buf,sizeof(unsigned char),4,fp) < 4) {
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
#endif


#ifdef OUTPUT_BIGENDIAN
unsigned int
Bigendian_fileio_read_uint (int fd) {
  unsigned int value = 0U;
  unsigned char buf[4];

  read(fd,buf,4);
  value = (buf[0] & 0xff);
  value <<= 8;
  value |= (buf[1] & 0xff);
  value <<= 8;
  value |= (buf[2] & 0xff);
  value <<= 8;
  value |= (buf[3] & 0xff);
  return value;
}
#else
unsigned int
Bigendian_fileio_read_uint (int fd) {
  unsigned int value = 0U;
  unsigned char buf[4];

  read(fd,buf,4);
  value = (buf[3] & 0xff);
  value <<= 8;
  value |= (buf[2] & 0xff);
  value <<= 8;
  value |= (buf[1] & 0xff);
  value <<= 8;
  value |= (buf[0] & 0xff);
  return value;
}
#endif


/************************************************************************
 *   Long unsigned int
 ************************************************************************/

#ifdef HAVE_64_BIT

UINT8
Bigendian_convert_uint8 (UINT8 littleendian) {
  UINT8 bigendian;
  unsigned char byte1, byte2, byte3, byte4, byte5, byte6, byte7;

  bigendian = littleendian & 0xff;
  byte1 = (littleendian >>= 8);
  byte2 = (littleendian >>= 8);
  byte3 = (littleendian >>= 8);
  byte4 = (littleendian >>= 8);
  byte5 = (littleendian >>= 8);
  byte6 = (littleendian >>= 8);
  byte7 = (littleendian >>= 8);

  /* bigendian = byte0; */
  bigendian <<= 8;
  bigendian |= (byte1 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte2 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte3 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte4 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte5 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte6 & 0xff);
  bigendian <<= 8;
  bigendian |= (byte7 & 0xff);

  return bigendian;
}


#ifdef OUTPUT_BIGENDIAN
void
Bigendian_write_uint8 (UINT8 value, int fd) {
  unsigned char buf[8];

  buf[7] = value & 0xff;
  buf[6] = (value >>= 8) & 0xff;
  buf[5] = (value >>= 8) & 0xff;
  buf[4] = (value >>= 8) & 0xff;
  buf[3] = (value >>= 8) & 0xff;
  buf[2] = (value >>= 8) & 0xff;
  buf[1] = (value >>= 8) & 0xff;
  buf[0] = (value >>= 8) & 0xff;
  write(fd,buf,8);
  return;
}
#else
void
Bigendian_write_uint8 (UINT8 value, int fd) {
  unsigned char buf[8];

  buf[0] = (unsigned char) (value & 0xff);
  buf[1] = (unsigned char) ((value >>= 8) & 0xff);
  buf[2] = (unsigned char) ((value >>= 8) & 0xff);
  buf[3] = (unsigned char) ((value >>= 8) & 0xff);
  buf[4] = (unsigned char) ((value >>= 8) & 0xff);
  buf[5] = (unsigned char) ((value >>= 8) & 0xff);
  buf[6] = (unsigned char) ((value >>= 8) & 0xff);
  buf[7] = (unsigned char) ((value >>= 8) & 0xff);
  write(fd,buf,8);
}
#endif



#ifdef OUTPUT_BIGENDIAN
size_t
Bigendian_fwrite_uint8 (UINT8 value, FILE *fp) {
  unsigned char buf[8];

  buf[7] = value & 0xff;
  buf[6] = (value >>= 8) & 0xff;
  buf[5] = (value >>= 8) & 0xff;
  buf[4] = (value >>= 8) & 0xff;
  buf[3] = (value >>= 8) & 0xff;
  buf[2] = (value >>= 8) & 0xff;
  buf[1] = (value >>= 8) & 0xff;
  buf[0] = (value >>= 8) & 0xff;
  if (fwrite(buf,sizeof(unsigned char),8,fp) == 0) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    return 1;
  }
}
#else
size_t
Bigendian_fwrite_uint8 (UINT8 value, FILE *fp) {
  unsigned char buf[8];

  buf[0] = value & 0xff;
  buf[1] = (value >>= 8) & 0xff;
  buf[2] = (value >>= 8) & 0xff;
  buf[3] = (value >>= 8) & 0xff;
  buf[4] = (value >>= 8) & 0xff;
  buf[5] = (value >>= 8) & 0xff;
  buf[6] = (value >>= 8) & 0xff;
  buf[7] = (value >>= 8) & 0xff;
  if (fwrite(buf,sizeof(unsigned char),8,fp) == 0) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    return 1;
  }
}
#endif

#ifdef OUTPUT_BIGENDIAN
size_t
Bigendian_fwrite_uint8s (UINT8 *array, int n, FILE *fp) {
  unsigned char buf[8];
  UINT8 value;
  int i;
  
  for (i = 0; i < n; i++) {
    value = array[i];
    buf[7] = value & 0xff;
    buf[6] = (value >>= 8) & 0xff;
    buf[5] = (value >>= 8) & 0xff;
    buf[4] = (value >>= 8) & 0xff;
    buf[3] = (value >>= 8) & 0xff;
    buf[2] = (value >>= 8) & 0xff;
    buf[1] = (value >>= 8) & 0xff;
    buf[0] = (value >>= 8) & 0xff;
    if (fwrite(buf,sizeof(unsigned char),8,fp) == 0) {
      /* Should set error indicator for stream and set errno */
      return 0;
    }
  }
  return n;
}
#else
size_t
Bigendian_fwrite_uint8s (UINT8 *array, int n, FILE *fp) {
  unsigned char buf[8];
  UINT8 value;
  int i;
  
  for (i = 0; i < n; i++) {
    value = array[i];
    buf[0] = value & 0xff;
    buf[1] = (value >>= 8) & 0xff;
    buf[2] = (value >>= 8) & 0xff;
    buf[3] = (value >>= 8) & 0xff;
    buf[4] = (value >>= 8) & 0xff;
    buf[5] = (value >>= 8) & 0xff;
    buf[6] = (value >>= 8) & 0xff;
    buf[7] = (value >>= 8) & 0xff;
    if (fwrite(buf,sizeof(unsigned char),8,fp) == 0) {
      /* Should set error indicator for stream and set errno */
      return 0;
    }
  }
  return n;
}
#endif


#ifdef OUTPUT_BIGENDIAN
size_t
Bigendian_fread_uint8 (UINT8 *value, FILE *fp) {
  unsigned char buf[8];

  if (fread(buf,sizeof(unsigned char),8,fp) < 8) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    *value = (buf[0] & 0xff);
    *value <<= 8;
    *value |= (buf[1] & 0xff);
    *value <<= 8;
    *value |= (buf[2] & 0xff);
    *value <<= 8;
    *value |= (buf[3] & 0xff);
    *value <<= 8;
    *value |= (buf[4] & 0xff);
    *value <<= 8;
    *value |= (buf[5] & 0xff);
    *value <<= 8;
    *value |= (buf[6] & 0xff);
    *value <<= 8;
    *value |= (buf[7] & 0xff);
    return 1;
  }
}
#else
size_t
Bigendian_fread_uint8 (UINT8 *value, FILE *fp) {
  unsigned char buf[8];

  if (fread(buf,sizeof(unsigned char),8,fp) < 8) {
    /* Should set error indicator for stream and set errno */
    return 0;
  } else {
    *value = (buf[7] & 0xff);
    *value <<= 8;
    *value = (buf[6] & 0xff);
    *value <<= 8;
    *value = (buf[5] & 0xff);
    *value <<= 8;
    *value = (buf[4] & 0xff);
    *value <<= 8;
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
#endif


#ifdef OUTPUT_BIGENDIAN
size_t
Bigendian_fread_uint8s (UINT8 *array, int n, FILE *fp) {
  unsigned char buf[8];
  UINT8 value;
  int i;

  for (i = 0; i < n; i++) {
    if (fread(buf,sizeof(unsigned char),8,fp) < 8) {
      /* Should set error indicator for stream and set errno */
      return 0;
    } else {
      value = (buf[0] & 0xff);
      value <<= 8;
      value |= (buf[1] & 0xff);
      value <<= 8;
      value |= (buf[2] & 0xff);
      value <<= 8;
      value |= (buf[3] & 0xff);
      value <<= 8;
      value |= (buf[4] & 0xff);
      value <<= 8;
      value |= (buf[5] & 0xff);
      value <<= 8;
      value |= (buf[6] & 0xff);
      value <<= 8;
      value |= (buf[7] & 0xff);

      array[i] = value;
    }
  }
  return n;
}
#else
size_t
Bigendian_fread_uint8s (UINT8 *array, int n, FILE *fp) {
  unsigned char buf[8];
  UINT8 value;
  int i;

  for (i = 0; i < n; i++) {
    if (fread(buf,sizeof(unsigned char),8,fp) < 8) {
      /* Should set error indicator for stream and set errno */
      return 0;
    } else {
      value = (buf[7] & 0xff);
      value <<= 8;
      value = (buf[6] & 0xff);
      value <<= 8;
      value = (buf[5] & 0xff);
      value <<= 8;
      value = (buf[4] & 0xff);
      value <<= 8;
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
#endif


#ifdef OUTPUT_BIGENDIAN
UINT8
Bigendian_fileio_read_uint8 (int fd) {
  UINT8 value = 0LU;
  unsigned char buf[8];

  read(fd,buf,8);
  value = (buf[0] & 0xff);
  value <<= 8;
  value |= (buf[1] & 0xff);
  value <<= 8;
  value |= (buf[2] & 0xff);
  value <<= 8;
  value |= (buf[3] & 0xff);
  value <<= 8;
  value |= (buf[4] & 0xff);
  value <<= 8;
  value |= (buf[5] & 0xff);
  value <<= 8;
  value |= (buf[6] & 0xff);
  value <<= 8;
  value |= (buf[7] & 0xff);
  return value;
}
#else
UINT8
Bigendian_fileio_read_uint8 (int fd) {
  UINT8 value = 0LU;
  unsigned char buf[8];

  read(fd,buf,8);
  value = (buf[7] & 0xff);
  value <<= 8;
  value = (buf[6] & 0xff);
  value <<= 8;
  value = (buf[5] & 0xff);
  value <<= 8;
  value = (buf[4] & 0xff);
  value <<= 8;
  value = (buf[3] & 0xff);
  value <<= 8;
  value |= (buf[2] & 0xff);
  value <<= 8;
  value |= (buf[1] & 0xff);
  value <<= 8;
  value |= (buf[0] & 0xff);
  return value;
}
#endif


#endif /* HAVE_64_BIT */

