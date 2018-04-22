#!/bin/bash                                                                  
# 
# Create all output results 
#

# Useful shell settings:

# abort the script if a command fails
set -e

# abort the script if an unitialized shell variable is used
set -u

# make sure the code is up to date

pushd src
make clean
make
popd

# generate the result pictures

# src/imgpro input/test_D_face01.jpg output/test_D_blur_2.jpg -blur 2.0

# src/imgpro input/globos_de_colores.jpg output/globos_brighntess_1.0.jpg \
#     -brightness 1.0

# src/imgpro input/globos_de_colores.jpg output/globos_brighntess_1.5.jpg \
#     -brightness 1.5

# src/imgpro input/test_D_face01.jpg output/test_D_01_blur_in_harris.jpg -harris \
#     2.0 100

# src/imgpro input/test_D_face02.jpg output/test_D_02_harris.jpg -harris \
#     2.0 100    

# src/imgpro input/test_D_face01.jpg output/test_D_trans_100.jpg -matchTranslation \
#     input/test_D_face02.jpg 100 0.1

# src/imgpro input/test_D_face01.jpg output/test_D_face01.jpg -videoStabilization 285 ait_test_in/pictures ait_test_stabilized/
src/imgpro input/test_D_face01.jpg output/test_D_face01.jpg -videoStabilization 222 testSequence/ testSequenceStabilized/