/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */

#ifndef _IO_PNG_H
#define _IO_PNG_H

#ifdef __cplusplus
extern "C" {
#endif

#define IO_PNG_VERSION "0.20110825"

#include <stddef.h>

    /* io_png.c */
    char *io_png_info(void);
    unsigned char *io_png_read_u8(const char *fname, size_t *nxp, size_t *nyp,
                                  size_t *ncp);
    unsigned char *io_png_read_u8_rgb(const char *fname, size_t *nxp, size_t *nyp);
    unsigned char *io_png_read_u8_gray(const char *fname, size_t *nxp, size_t *nyp);
    float *io_png_read_f32(const char *fname, size_t *nxp, size_t *nyp,
                           size_t *ncp);
    float *io_png_read_f32_rgb(const char *fname, size_t *nxp, size_t *nyp);
    float *io_png_read_f32_gray(const char *fname, size_t *nxp, size_t *nyp);
    int io_png_write_u8(const char *fname, const unsigned char *data,
                        size_t nx, size_t ny, size_t nc);
    int io_png_write_f32(const char *fname, const float *data,
                         size_t nx, size_t ny, size_t nc);
    int ami_write_png(char name[200], unsigned char *red, unsigned char *green, unsigned char *blue, int width, int height);
    int ami_read_png(char name[200], unsigned char **red, unsigned char **green, unsigned char **blue, int *width, int *height);



#ifdef __cplusplus
}
#endif

#endif /* !_IO_PNG_H */
