#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([SPQR], [1.0], [spoblete@dlab.cl])

# Checks for programs.
AC_PROG_CC

AC_LANG([C])

#if test -z "$CFLAGS"
#then
#  CFLAGS=-O3
#fi


AC_CHECK_HEADER( [mpi.h], [
  AC_CHECK_FUNC([MPI_Init],[
    AC_DEFINE([MPI])
  ])])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
