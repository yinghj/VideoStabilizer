# Video Stabilizer - Computer Vision for Digital Post-Production Final Project
### AIT-Budapest Program, Fall 2016

## FILE STRUCTURE

- Main implementation in src/R2Image.cpp.
    Functionalities implemented include: Gaussian Blur, Sharpen, Sobel Filters, Harris Corner Detector, Translation Tracking Marker, Homography Blender, Gaussian Smooth for Vector Sequences, and Video Stabilizer which utilizes and combines the previous functions. 

- Other files:
    src/ - Directory with source code
    Makefile - Unix/Mac makefile for building the project with "make". 
    imagepro.[vcproj/sln/suo] - Project file for Visual Studio 2005 on Windows
    imgpro.cpp - Main program, parses the command line arguments, and calls the appropriate image functions
    R2Image.[cpp/h] - Image class with processing functions (this is the only file that you need to edit)
    R2Pixel.[cpp/h] - Pixel class 
    R2/ - A library of useful 2D geometric primitives
    jpeg/ - A library for reading/writing JPEG files
    input/ - Contains example input images. 
    output/ - Es empty to start -- it will contain the images produced by your program (see below)
    runme.bat - a script (for Windows) that you will fill in to demonstrate execution of your program
    runme.sh - same as runme.bat, but for Mac OS X

## COMPILATION

If you are developing on a Windows machine and have Visual Studio installed, use the provided project solution file (assn1.sln) in the src/ directory to build the program. If you are developing on a Mac or Linux machine, cd into the src/ directory and type "make". In either case, an executable called imgpro (or imgpro.exe) will be created in
the src/ directory.

## RUNNING VIDEO STABILIZER

Modify the runme bash with appropriate commands. 

For Mac or Linux Users, 
Use ffmepg to turn a video into a sequence of image frames, then
add 
  "src/imgpro <foo - some random input image name> <bar - some random out image name> -videoStabilization <int: number of frames> <folder: input images> <folder: output images>" to runme.sh, then
use ffmepg to turn the image sequence back into a video.

