dnl -*- Autoconf -*-

AC_INIT([BSseq_pipeline], [0.0.0])

AC_ARG_VAR(GENOME_DIR, full path to genome directory)
if test x"$GENOME_DIR" == x"" ; then
  AC_MSG_ERROR([Please specify GENOME_DIR, the full path to the genome directory.])
fi

AC_ARG_VAR(GENOME_INDEX, full path to genome index directory)
if test x"$GENOME_INDEX" == x"" ; then
  AC_MSG_ERROR([Please specify GENOME_INDEX, the full path to the genome index directory.])
fi

AC_ARG_VAR(GENOME_CHROM_INFO, full path to hg19.chromInfo.txt file)
if test x"$GENOME_CHROM_INFO" == x"" ; then
  AC_MSG_ERROR([Please specify GENOME_CHROM_INFO, the full path to the hg19.chromInfo.txt file.])
fi

AC_ARG_VAR(ADAPTER, full path to bbmap adapter file)
if test x"$ADAPTER" == x"" ; then
  AC_MSG_ERROR([Please specify ADAPTER, the full path to the bbmap adapter file (e.g. truseq.fa.gz).])
fi

dnl Check for required programmes and store their full path in the
dnl given variables.  The variables are used to substitute
dnl placeholders in the scripts.

dnl The variables only need to be set manually when they cannot be
dnl found via the PATH environment variable.
AC_DEFUN([find_or_override_prog],
[AC_ARG_VAR($1, override location of $2 executable)dnl
AS_IF([test -z "$$1"],dnl
      [AC_PATH_PROG([$1], [$2], no)dnl
      if test x"$$1" == x"no" ; then
        AC_MSG_ERROR([Could not find $2.])
      fi
      ],dnl
      [AC_MSG_NOTICE([Using $$1 as $2 executable.])])])

find_or_override_prog(BEDTOOLS,             bedtools)
find_or_override_prog(BOWTIE,               bowtie2)
find_or_override_prog(SAMTOOLS,             samtools)
find_or_override_prog(BISMARK,              bismark)
find_or_override_prog(DEDUP_BISMARK,        deduplicate_bismark)
find_or_override_prog(FASTQC,               fastqc)
find_or_override_prog(BBDUK,                bbduk2.sh)
find_or_override_prog(TRIM_GALORE,          trim_galore)
find_or_override_prog(CUTADAPT,             cutadapt)
find_or_override_prog(BED_GRAPH_TO_BIG_WIG, bedGraphToBigWig)

dnl Substitute placeholders to generate these target files
AC_CONFIG_FILES([config.json])
AC_OUTPUT
