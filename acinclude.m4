M4 Source Code

m4_define([DATE],['`date +%Y-%m-%d`'])

m4_include([config/pagesize.m4])
m4_include([config/madvise-flags.m4])
m4_include([config/mmap-flags.m4])
m4_include([config/acx_mmap_fixed.m4])
m4_include([config/acx_mmap_variable.m4])
m4_include([config/acx_pthread.m4])
m4_include([config/builtin.m4])
m4_include([config/struct-stat64.m4])
m4_include([config/expand.m4])
m4_include([config/perl.m4])

m4_include([config/fopen.m4])
m4_include([config/popcnt.m4])
m4_include([config/asm-bsr.m4])
