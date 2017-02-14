;;; BSseq_pipeline
;;; Copyright © 2017 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of BSseq_pipeline.
;;;
;;; BSseq_pipeline is free software under the GPL version 3 or any
;;; later version; see the LICENSE file for details.

;;; Run the following command to enter a development environment for
;;; this pipeline:
;;;
;;;  $ git clone https://github.com/bimsbbioinfo/guix-bimsb guix-bimsb
;;;  $ export GUIX_PACKAGE_PATH=$PWD/guix-bimsb
;;;  $ guix environment -l guix.scm

(use-modules ((guix licenses) #:prefix license:)
             (guix packages)
             (guix download)
             (guix utils)
             (guix build-system gnu)
             (gnu packages)
             (gnu packages autotools)
             (gnu packages bioinformatics)
             (bimsb packages staging))

(define-public bseq-pipeline/devel
  (package
    (name "bseq-pipeline-development")
    (version "0.0.0")
    (source #f)
    (build-system gnu-build-system)
    (arguments
     `(#:phases
       (modify-phases %standard-phases
         (add-before 'configure 'autoconf
           (lambda _
             (zero? (system* "autoreconf" "-vif")))))))
    (inputs
     `(("bedtools" ,bedtools)
       ("bowtie2" ,bowtie)
       ("samtools" ,samtools)
       ("bismark" ,bismark)
       ("trim-galore" ,trim-galore)
       ("cutadapt" ,cutadapt)
       ("bbmap" ,bbmap)
       ("fastqc" ,fastqc)
       ("kentutils" ,kentutils)))
    (native-inputs
     `(("autoconf" ,autoconf)
       ("automake" ,automake)))
    (home-page "https://github.com/BIMSBbioinfo/makeNGSnake.git")
    (synopsis "TODO")
    (description "TODO")
    (license license:gpl3+)))

bseq-pipeline/devel
