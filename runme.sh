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

# src/imgpro input/globos_de_colores.jpg output/globos_brighntess_0.5.jpg \
#     -brightness 0.5

# src/imgpro input/globos_de_colores.jpg output/globos_brighntess_1.0.jpg \
#     -brightness 1.0

# src/imgpro input/globos_de_colores.jpg output/globos_brighntess_1.5.jpg \
#     -brightness 1.5

# src/imgpro input/test_D_face01.jpg output/test_D_01_harris.jpg -harris \
#     2.0 100

# src/imgpro input/test_D_face02.jpg output/test_D_02_harris.jpg -harris \
#     2.0 100    

src/imgpro input/test_D_face01.jpg output/test_D_trans.jpg -matchTranslation \
    input/test_D_face02.jpg   

# src/imgpro -videoStabilization 222 testSequence/