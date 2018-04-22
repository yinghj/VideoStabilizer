// Source file for image class

// Include files
#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <ctime>
#include <numeric>   
#include <stdlib.h>
#include <float.h>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////


R2Image::
R2Image(void)
  : pixels(NULL),
    npixels(0),
    width(0),
    height(0)
{
}

R2Image::
R2Image(const char *filename)
  : pixels(NULL),
    npixels(0),
    width(0),
    height(0)
{
  // Read image
  Read(filename);
}

R2Image::
R2Image(int width, int height)
  : pixels(NULL),
    npixels(width * height),
    width(width),
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}

R2Image::
R2Image(int width, int height, const R2Pixel *p)
  : pixels(NULL),
    npixels(width * height),
    width(width),
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = p[i];
}

R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width),
    height(image.height)

{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = image.pixels[i];
}

R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}

void R2Image::
svdTest(void)
{
	// fit a 2D conic to five points
	R2Point p1(1.2,3.5);
	R2Point p2(2.1,2.2);
	R2Point p3(0.2,1.6);
	R2Point p4(0.0,0.5);
	R2Point p5(-0.2,4.2);

	// build the 5x6 matrix of equations
	double** linEquations = dmatrix(1,5,1,6);

	linEquations[1][1] = p1[0]*p1[0];
	linEquations[1][2] = p1[0]*p1[1];
	linEquations[1][3] = p1[1]*p1[1];
	linEquations[1][4] = p1[0];
	linEquations[1][5] = p1[1];
	linEquations[1][6] = 1.0;

	linEquations[2][1] = p2[0]*p2[0];
	linEquations[2][2] = p2[0]*p2[1];
	linEquations[2][3] = p2[1]*p2[1];
	linEquations[2][4] = p2[0];
	linEquations[2][5] = p2[1];
	linEquations[2][6] = 1.0;

	linEquations[3][1] = p3[0]*p3[0];
	linEquations[3][2] = p3[0]*p3[1];
	linEquations[3][3] = p3[1]*p3[1];
	linEquations[3][4] = p3[0];
	linEquations[3][5] = p3[1];
	linEquations[3][6] = 1.0;

	linEquations[4][1] = p4[0]*p4[0];
	linEquations[4][2] = p4[0]*p4[1];
	linEquations[4][3] = p4[1]*p4[1];
	linEquations[4][4] = p4[0];
	linEquations[4][5] = p4[1];
	linEquations[4][6] = 1.0;

	linEquations[5][1] = p5[0]*p5[0];
	linEquations[5][2] = p5[0]*p5[1];
	linEquations[5][3] = p5[1]*p5[1];
	linEquations[5][4] = p5[0];
	linEquations[5][5] = p5[1];
	linEquations[5][6] = 1.0;

	printf("\n Fitting a conic to five points:\n");
	printf("Point #1: %f,%f\n",p1[0],p1[1]);
	printf("Point #2: %f,%f\n",p2[0],p2[1]);
	printf("Point #3: %f,%f\n",p3[0],p3[1]);
	printf("Point #4: %f,%f\n",p4[0],p4[1]);
	printf("Point #5: %f,%f\n",p5[0],p5[1]);

	// compute the SVD
	double singularValues[7]; // 1..6
	double** nullspaceMatrix = dmatrix(1,6,1,6);
	svdcmp(linEquations, 5, 6, singularValues, nullspaceMatrix);

	// get the result
	printf("\n Singular values: %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6]);

	// find the smallest singular value:
	int smallestIndex = 1;
	for(int i=2;i<7;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

	// solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
	printf("Conic coefficients: %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex]);

	// make sure the solution is correct:
	printf("Equation #1 result: %f\n",	p1[0]*p1[0]*nullspaceMatrix[1][smallestIndex] +
										p1[0]*p1[1]*nullspaceMatrix[2][smallestIndex] +
										p1[1]*p1[1]*nullspaceMatrix[3][smallestIndex] +
										p1[0]*nullspaceMatrix[4][smallestIndex] +
										p1[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #2 result: %f\n",	p2[0]*p2[0]*nullspaceMatrix[1][smallestIndex] +
										p2[0]*p2[1]*nullspaceMatrix[2][smallestIndex] +
										p2[1]*p2[1]*nullspaceMatrix[3][smallestIndex] +
										p2[0]*nullspaceMatrix[4][smallestIndex] +
										p2[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #3 result: %f\n",	p3[0]*p3[0]*nullspaceMatrix[1][smallestIndex] +
										p3[0]*p3[1]*nullspaceMatrix[2][smallestIndex] +
										p3[1]*p3[1]*nullspaceMatrix[3][smallestIndex] +
										p3[0]*nullspaceMatrix[4][smallestIndex] +
										p3[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #4 result: %f\n",	p4[0]*p4[0]*nullspaceMatrix[1][smallestIndex] +
										p4[0]*p4[1]*nullspaceMatrix[2][smallestIndex] +
										p4[1]*p4[1]*nullspaceMatrix[3][smallestIndex] +
										p4[0]*nullspaceMatrix[4][smallestIndex] +
										p4[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #5 result: %f\n",	p5[0]*p5[0]*nullspaceMatrix[1][smallestIndex] +
										p5[0]*p5[1]*nullspaceMatrix[2][smallestIndex] +
										p5[1]*p5[1]*nullspaceMatrix[3][smallestIndex] +
										p5[0]*nullspaceMatrix[4][smallestIndex] +
										p5[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	R2Point test_point(0.34,-2.8);

	printf("A point off the conic: %f\n",	test_point[0]*test_point[0]*nullspaceMatrix[1][smallestIndex] +
											test_point[0]*test_point[1]*nullspaceMatrix[2][smallestIndex] +
											test_point[1]*test_point[1]*nullspaceMatrix[3][smallestIndex] +
											test_point[0]*nullspaceMatrix[4][smallestIndex] +
											test_point[1]*nullspaceMatrix[5][smallestIndex] +
											nullspaceMatrix[6][smallestIndex]);

	return;
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
  // Brighten the image by multiplying each pixel component by the factor.
  // This is implemented for you as an example of how to access and set pixels
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) *= factor;
      Pixel(i,j).Clamp();
    }
  }
}

void R2Image::
SobelX(void)
{
	// Apply the Sobel oprator to the image in X direction

	// Define a SobelX kernel.
	float weight[3][3] = { { -1,  0,  1 },
						 { -2,  0,  2 },
						 { -1,  0,  1 } };

	R2Image tempImg(width, height);

	// Perform Sobel operation and combine with the original image
	for (int i = 1; i < width - 1; i++) {
		for (int j = 1; j < height - 1; j++) {
			R2Pixel accumulator = R2null_pixel;
			for (int y = -1; y < 2; y++) {
				for (int x = -1; x < 2; x++) {
					accumulator += Pixel(i + x, j + y) * weight[y+1][x+1];
				}
			}
			tempImg.Pixel(i, j) = accumulator;
		}
	}

	*this = tempImg;
}

void R2Image::
SobelY(void)
{
	// Apply the Sobel oprator to the image in Y direction

	// Define a SobelY kernel.
	int weight[3][3] = { { -1,  -2,  -1 },
						 { 0,  0,  0 },
						 { 1,  2,  1 } };

	R2Image tempImg(width, height);

	// Perform Sobel operation and combine with the original image
	for (int i = 1; i < width - 1; i++) {
		for (int j = 1; j < height - 1; j++) {
			tempImg.Pixel(i, j) = R2null_pixel;
			for (int y = -1; y < 2; y++) {
				for (int x = -1; x < 2; x++) {
					tempImg.Pixel(i, j) += Pixel(i + x, j + y) * weight[y + 1][x + 1];
				}
			}
		}
	}

	*this = tempImg;
}

void R2Image::
LoG(void)
{
  // Apply the LoG oprator to the image

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "LoG() not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{
  // Changes the saturation of an image
  // Find a formula that changes the saturation without affecting the image brightness

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "ChangeSaturation(%g) not implemented\n", factor);
}


// Linear filtering ////////////////////////////////////////////////
void R2Image::
Blur(double sigma)
{
  // Gaussian blur of the image. Separable solution is preferred
	static const double PI = 3.1415;
	int sigma_int = sigma;
	int KERNEL_SIZE = sigma_int * 6 + 1;
	double *kernel = new double[KERNEL_SIZE];
	double sum = 0.0;
	for (int i = 0; i < KERNEL_SIZE; i++) {
		int x = i - 3 * sigma_int;
		kernel[i] = 1 / sqrt(2 * PI * sigma * sigma) * exp((0 - x * x) / (2 * sigma * sigma));
		sum += kernel[i];
	}
	// Normalize the kernel.
	for (int i = 0; i < KERNEL_SIZE; i++) {
		kernel[i] = kernel[i] / sum;
	}

	// Apply filter x-direction.
	R2Image tempImgX(width, height);
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			for (int x = 0 - 3 * sigma_int; x <= 3 * sigma_int; x++) {
				int new_X = i + x;
				if (new_X < 0) {
					new_X = 0;
				} else if (new_X >= width) {
					new_X = width - 1;
				}
				tempImgX.Pixel(i, j) += Pixel(new_X, j) * kernel[x + 3 * sigma_int];
			}
		}
	}

	// Apply filter y-direction.
	R2Image tempImgY(width, height);
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			for (int y = 0 - 3 * sigma_int; y <= 3 * sigma_int; y++) {
				int new_Y = j + y;
				if (new_Y < 0) {
					new_Y = 0;
				}
				else if (new_Y >= height) {
					new_Y = height - 1;
				}
				tempImgY.Pixel(i, j) += tempImgX.Pixel(i, new_Y) * kernel[y + 3 * sigma_int];
			}
		}
	}

	*this = tempImgY;
}

// Helper functions

bool custom_sort_func(const HarrisPixel &left, const HarrisPixel &right) {
    return left.harrisScore > right.harrisScore;
}

float ssd(R2Pixel a, R2Pixel b) {
	return ((a[0] - b[0]) * (a[0] - b[0])
		+ (a[1] - b[1]) * (a[1] - b[1]) 
		+ (a[2] - b[2]) * (a[2] - b[2])
		+ (a[3] - b[3]) * (a[3] - b[3]));
}

void R2Image::
line(int x0, int x1, int y0, int y1, float r, float g, float b) {
	if (x0 > x1) {
		int x = y1;
		y1 = y0;
		y0 = x;

		x = x1;
		x1 = x0;
		x0 = x;
	}
	int deltax = x1 - x0;
	int deltay = y1 - y0;
	float error = 0;
	float deltaerr = 0.0;
	if (deltax!=0) deltaerr = fabs(float(float(deltay) / deltax));    // Assume deltax != 0 (line is not vertical),
	// note that this division needs to be done in a way that preserves the fractional part
	int y = y0;
	for (int x= x0; x <= x1; x++) {
		Pixel(x,y).Reset(r,g,b,1.0);
		error = error + deltaerr;
		if (error >= 0.5) {
			if (deltay > 0) y = y + 1;
			else y = y - 1;

			error = error - 1.0;
		}
	}
	if (x0 > 3 && x0 < width - 3 && y0 > 3 && y0 < height - 3) {
		for(int x = x0-3; x <= x0+3; x++) {
			for(int y = y0-3; y<=y0+3; y++) {
				Pixel(x,y).Reset(r,g,b,1.0);
			}
		}
	}
	if (x1 > 3 && x1 < width - 3 && y1 > 3 && y1 < height - 3) {
		for(int x = x1-3; x <= x1+3; x++) {
			for(int y = y1-3; y<=y1+3; y++) {
				Pixel(x,y).Reset(r,b,g,1.0);
			}
		}
	}
}

// Bound normalizes the edges of a filtered image that cannot be filtered by
// the kernel.
int bound(int value, int max) {
	if (value < 0) {
		return 0;
	}
	if (value >= max) {
		return max - 1;
	}
	return value;
}

HarrisPixel R2Image::
Search(R2Image originalImage, R2Image otherImage, HarrisPixel featurePixel, float search_win_perc) {
	// Search matching window within close neighborhood of the original window
	int searchWidth = search_win_perc * otherImage.width;
	int searchHeight = search_win_perc * otherImage.height;
	HarrisPixel matchingHarrisPixel;

	// Find the matching window (minimize sum squared distance)
	const int WINDOW_SIZE = 5;

	float minSSD = FLT_MAX;

	int searchStart_x = featurePixel.posx - searchWidth;
	int searchEnd_x = bound(featurePixel.posx + searchWidth, width);

	int searchStart_y = featurePixel.posy - searchHeight;
	int searchEnd_y = bound(featurePixel.posy + searchHeight, height);

	for (int i = searchStart_x; i <= searchEnd_x; i++) {
		for (int j = searchStart_y; j <= searchEnd_y; j++) {
			// Calculate average pixel-wise SSD in the window.
			float newSSD = 0.0;
			int pix_count = 0;

			for (int x = 0 - WINDOW_SIZE; x <= WINDOW_SIZE; x++) {
				for (int y = 0 - WINDOW_SIZE; y <= WINDOW_SIZE; y++) {
					// Watch for out-of-bounds situation.
					newSSD += ssd(
						originalImage.Pixel(
							bound(featurePixel.posx + x, width), 
							bound(featurePixel.posy + y, height)),
						otherImage.Pixel(
							bound(i + x, width), 
							bound(j + y, width)));
					pix_count += 1;
				}
			}

			newSSD = newSSD / pix_count;
			
			// Check if minSSD.
			if (newSSD < minSSD) {
				matchingHarrisPixel.posx = i;
				matchingHarrisPixel.posy = j;
				minSSD = newSSD;
			}
		}
	}

	return matchingHarrisPixel;
}

std::vector<HarrisPixel> R2Image::
Harris(double sigma, int num_top_features)
{
	// Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
	// Output should be 50% grey at flat regions, white at corners and black/dark near edges

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			double grey = (Pixel(i, j)[0] + Pixel(i, j)[1] + Pixel(i, j)[2]) / 3.0;
			Pixel(i, j).SetRed(grey);
			Pixel(i, j).SetBlue(grey);
			Pixel(i, j).SetGreen(grey);
		}
	}

	R2Image tempSobelX(*this);
	tempSobelX.SobelX();

	R2Image tempSobelY(*this);
	tempSobelY.SobelY();

	R2Image tempSobelXY(width, height);

	// Calculate Autocorrection Matrix [A].
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			// I_x * I_y
			tempSobelXY.Pixel(i, j)[0] = tempSobelX.Pixel(i, j)[0] * tempSobelY.Pixel(i, j)[0];
			// I_x ^ 2
			tempSobelX.Pixel(i, j)[0] = tempSobelX.Pixel(i, j)[0] * tempSobelX.Pixel(i, j)[0];
			// I_y ^ 2
			tempSobelY.Pixel(i, j)[0] = tempSobelY.Pixel(i, j)[0] * tempSobelY.Pixel(i, j)[0];
		}
	}


	tempSobelX.Blur(sigma);
	tempSobelY.Blur(sigma);
	tempSobelXY.Blur(sigma);

	R2Image harris(width, height);

	// Calculate Harris Score for each pixel.
	printf("Calculating Harris score for each pixel...\r");
	fflush(stdout);
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			// Harris Score = Det(A) - alpha * (Tr(A)^2)
			harris.Pixel(i, j)[0] = tempSobelX.Pixel(i, j)[0] * tempSobelY.Pixel(i, j)[0]
								  - tempSobelXY.Pixel(i, j)[0] * tempSobelXY.Pixel(i, j)[0]
								  - 0.04 * ((tempSobelX.Pixel(i, j)[0] + tempSobelY.Pixel(i, j)[0])
										 *  (tempSobelX.Pixel(i, j)[0] + tempSobelY.Pixel(i, j)[0]));
			harris.Pixel(i, j)[1] = harris.Pixel(i, j)[0];
			harris.Pixel(i, j)[2] = harris.Pixel(i, j)[0];
		}
	}

	// Choose local maximums as top Harris features.
	std::vector<HarrisPixel> topFeatures;

	for (int x = 0; x < num_top_features; x++) {
		HarrisPixel emptyHarrisPixel;
		emptyHarrisPixel.posx = -100;
		emptyHarrisPixel.posy = -100;
		emptyHarrisPixel.harrisScore = -100.0;
		topFeatures.push_back(emptyHarrisPixel);
	}

	// Only pick pixels that have high enough harris score and are far apart enough from other feature pixels.
	const int DIST_THRESHOLD = 16*16;
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {

			int minHarrisScore = topFeatures.back().harrisScore;

			HarrisPixel currentHarrisPixel;
			currentHarrisPixel.posx = i;
			currentHarrisPixel.posy = j;
			currentHarrisPixel.harrisScore = harris.Pixel(i, j)[0];

			bool topScore = 1;
			if (minHarrisScore > currentHarrisPixel.harrisScore) {
				topScore = 0;
			}

			bool apart = 1;
			if (topScore) {
				for (int x = 0; x < topFeatures.size(); x++) {
					// Check distance.
					if (((topFeatures[x].posx - i) * (topFeatures[x].posx - i) +
						(topFeatures[x].posy - j) * (topFeatures[x].posy - j) <= DIST_THRESHOLD)) {
						apart = 0;
						break;
					}
				}
			}

			if (topScore * apart) {
				topFeatures.pop_back();
				topFeatures.push_back(currentHarrisPixel);
				std::sort(std::begin(topFeatures), std::end(topFeatures), custom_sort_func);
			}
		}
	}

	return topFeatures;
}

void R2Image::
MarkHarrisFeatures(double sigma, int num_top_features) {
	R2Image harrisImage(*this);

	printf("Finding top %d Harris features...\r", num_top_features);
	fflush(stdout);

	std::vector<HarrisPixel> topFeatures = harrisImage.Harris(sigma, num_top_features);
	
	printf("Marking each Harris feature...\r");
	fflush(stdout);
	for (int i = 0; i < topFeatures.size(); i++) {
		int x = topFeatures[i].posx;
		int y = topFeatures[i].posy;
		// Mark Harris Features with small squares at feature pixels.
		// Define square color
		float r = 0.0;
		float g = 1.0;
		float b = 0.0;

		// Draw square marks
		for (int j = -2; j <= 2; j++) {
			// Watch for out-of-bound cases
			int sqr_x = bound(x + j, width);
			for (int k = -2; k <= 2; k++) {
				int sqr_y = bound(y + k, height);
				harrisImage.Pixel(sqr_x, sqr_y).Reset(r,g,b,1.0);
			}
		}
	}
	*this = harrisImage;
}


void R2Image::
Sharpen()
{
	// Sharpen an image using a linear filter. Use a kernel of your choosing.

	// Define a sharpening kernel.
	int weight[3][3] = { { 0,  -1,  0 },
						 { -1,  4,  -1 },
						 { 0,  -1,  0 } };

	R2Image tempImg(*this);

	// Perform sharpening and combine with the original image
	for (int i = 1; i < width - 1; i++) {
		for (int j = 1; j < height - 1; j++) {
			for (int y = -1; y < 2; y++) {
				for (int x = -1; x < 2; x++) {
					tempImg.Pixel(i, j) += Pixel(i + x, j + y) * weight[y + 1][x + 1];
				}
			}
			tempImg.Pixel(i, j).Clamp();
		}
	}

	*this = tempImg;
}

void R2Image::
matchOtherImageTranslated(R2Image* otherImage, float harris_sigma, int num_features_to_match, float search_win_perc)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching translation (pixel precision is OK), and blend the translated "otherImage"
	// into this image with a 50% opacity.

	R2Image markedImage(*otherImage);

	R2Image harrisImage(*this);
	std::vector<HarrisPixel> topFeatures = harrisImage.Harris(harris_sigma, num_features_to_match);

	std::vector<HarrisPixel> matchingFeatures;
	for (int x = 0; x < topFeatures.size(); x++) {

		printf("Searching for matching Harris feature No. %d...\r", x+1);
		fflush(stdout);

		HarrisPixel matchingPixel = Search(*this, *otherImage, topFeatures[x], search_win_perc);
		matchingFeatures.push_back(matchingPixel);
	}

	printf("\n");
	fflush(stdout);
	
	// RANSAC
	const int N = 10;
	// using sqr of length
	const float THRESHOLD = 9;
	int max_inliner = 0;
	int track_vector_x = 0;
	int track_vector_y = 0;

	for (int n = 0; n < N; n++) {
		int num_inliner = 0;

		int rand_index = rand() % num_features_to_match;

		HarrisPixel randPixel = matchingFeatures[rand_index];
		
		int u1 = randPixel.posx - topFeatures[rand_index].posx;
		int u2 = randPixel.posy - topFeatures[rand_index].posy;

		for (int x = 0; x < topFeatures.size(); x++) {
			int m = matchingFeatures[x].posx;
			int n = matchingFeatures[x].posy;
			int p = topFeatures[x].posx;
			int q = topFeatures[x].posy;
			
			int v1 = m - p;
			int v2 = n - q;
			
			if (((u1 - v1) * (u1 - v1) + (u2 - v2) * (u2 - v2)) <= THRESHOLD) {
				num_inliner++;
			}
		}
		
		if (num_inliner > max_inliner) {
			max_inliner = num_inliner;
			track_vector_x = u1;
			track_vector_y = u2;
		}
	}

	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			markedImage.Pixel(i, j) = markedImage.Pixel(i, j) * 0.6 + Pixel(i, j) * 0.4;
		}
	}

	for (int x = 0; x < topFeatures.size(); x++) {
		int m = matchingFeatures[x].posx;
		int n = matchingFeatures[x].posy;
		int p = topFeatures[x].posx;
		int q = topFeatures[x].posy;

		int v1 = m - p;
		int v2 = n - q;

		if (((track_vector_x - v1) * (track_vector_x - v1) + (track_vector_y - v2) * (track_vector_y - v2)) <= THRESHOLD) {
			markedImage.line(m, p, n, q, 0.0, 1.0, 0.0);
		}
		else {
			markedImage.line(m, p, n, q, 1.0, 0.0, 0.0);
		}
	}
	
	*this = markedImage;
}

std::vector<double> findHomographyMatrix(vector<int> & inputX, vector<int> & inputY, vector<int> & outputX, vector<int> & outputY) {
	std::vector<double> h;
	// construct A to compute Ah = 0;
	int num_points = inputX.size();
	double **A = new double*[num_points * 2 + 1];
	for (int i = 0; i < num_points * 2 + 1; i++) {
		A[i] = new double[10];
	}

	for (int i = 0; i < num_points; i ++) {
		int ax_index = i * 2 + 1;
		int ay_index = ax_index + 1;
		A[ax_index][1] = 0;
		A[ax_index][2] = 0;
		A[ax_index][3] = 0;
		A[ax_index][4] = 0 - inputX[i];
		A[ax_index][5] = 0 - inputY[i];
		A[ax_index][6] = 0 - 1;
		A[ax_index][7] = outputY[i] * inputX[i];
		A[ax_index][8] = outputY[i] * inputY[i];
		A[ax_index][9] = outputY[i];
		A[ay_index][1] = inputX[i];
		A[ay_index][2] = inputY[i];
		A[ay_index][3] = 1;
		A[ay_index][4] = 0;
		A[ay_index][5] = 0;
		A[ay_index][6] = 0;
		A[ay_index][7] = 0 - outputX[i] * inputX[i];
		A[ay_index][8] = 0 - outputX[i] * inputY[i];
		A[ay_index][9] = 0 - outputX[i];
	}

	double *w = new double[9 + 1];
	double **v = dmatrix(1, 9, 1, 9);

	svdcmp(A, num_points * 2, 9, w, v);
	
	int nullspaceCol = -1;
	double epsilon = 0.0000001;
	for (int i = 1; i < 9 + 1; i++) {
		if (sqrt(w[i] * w[i]) < epsilon){
			nullspaceCol = i;
			break;
		}
	}

	for (int row = 1; row < 9 + 1; row++) {
		if (sqrt(v[row][nullspaceCol] * v[row][nullspaceCol]) < epsilon) {
			h.push_back(0.0);
		}
		else {
			h.push_back(v[row][nullspaceCol]);
		}
	}

	free(v);
	free(w);

	return h;
}


void compute_centroid(double *c1x, double *c1y, double *c2x, double *c2y, double *s1, double *s2,
	vector<int> & inputX, vector<int> & inputY, vector<int> & outputX, vector<int> & outputY) {

	// first compute for the first image
	int n = inputX.size();
	int i;
	double xsum = 0;
	double ysum = 0;

	for (i = 0; i<n; i++) {
		xsum += inputX[i];
		ysum += inputY[i];
	}

	*c1x = xsum / n;
	*c1y = ysum / n;

	xsum = 0;
	ysum = 0;

	for (i = 0; i<n; i++) {
		xsum += outputX[i];
		ysum += outputY[i];
	}

	*c2x = xsum / n;
	*c2y = ysum / n;

	double dist1 = 0;
	double dist2 = 0;

	for (i = 0; i<n; i++) {
		dist1 += (inputX[i] - *c1x)*(inputX[i] - *c1x) +
			(inputY[i] - *c1y)*(inputY[i] - *c1y);
		dist2 += (outputX[i] - *c2x)*(outputX[i] - *c2x) +
			(outputY[i] - *c2y)*(outputY[i] - *c2y);
	}

	*s1 = sqrt(2.0*n / dist1);
	*s2 = sqrt(2.0*n / dist2);
}

void allocate_T(double **T, double cx, double cy, double s) {
	T[0][0] = s;
	T[1][0] = 0;
	T[2][0] = 0;

	T[0][1] = 0;
	T[1][1] = s;
	T[2][1] = 0;

	T[0][2] = -cx*s;
	T[1][2] = -cy*s;
	T[2][2] = 1;
}

void inverse(double **m, double **minv) {
	// computes the inverse of a matrix m
	double det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
		m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
		m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

	double invdet = 1 / det;

	minv[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
	minv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
	minv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
	minv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
	minv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
	minv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
	minv[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
	minv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
	minv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;
}

// Assume the matrices are of the right dimension
void matrix_multiply(double **first, double **second, double **result){
	int x, y, i;
	x = y = i = 0;

	for (y = 0; y<3; y++) // result->col
		for (x = 0; x<3; x++) // result->row
			for (i = 0; i<3; i++)         // which is second->row == first->column
				result[x][y] += first[x][i] * second[i][y];
}


void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching homography, and blend the transformed "otherImage" into this image with a 50% opacity.
	const int NUM_FEATURES_MATCH = 150;

	R2Image harrisImage(*this);
	std::vector<HarrisPixel> topFeatures = harrisImage.Harris(2.0, NUM_FEATURES_MATCH);

	std::vector<HarrisPixel> matchingFeatures;
	for (int x = 0; x < topFeatures.size(); x++) {
		HarrisPixel matchingPixel = Search(*this, *otherImage, topFeatures[x], 0.05 /** search_win_perc **/);
		matchingFeatures.push_back(matchingPixel);
	}

	R2Image blendImage(*this);

	const int N = 100;
	// using sqr of length
	const float THRESHOLD = 16;
	int max_inliner = 0;
	 double **homography_matrix = new double*[3];
	 for (int i = 0; i < 3; i++) {
		 homography_matrix[i] = new double[3];
	 }

	for (int n = 0; n < N; n++) {
		int num_inliner = 0;

		// Randomly pick 4 points
		int rand_index_1 = rand() % NUM_FEATURES_MATCH;
		int rand_index_2 = rand() % NUM_FEATURES_MATCH;
		while (rand_index_2 == rand_index_1) {
			rand_index_2 = rand() % NUM_FEATURES_MATCH;
		}
		int rand_index_3 = rand() % NUM_FEATURES_MATCH;
		while (rand_index_3 == rand_index_1 || rand_index_3 == rand_index_2) {
			rand_index_3 = rand() % NUM_FEATURES_MATCH;
		}
		int rand_index_4 = rand() % NUM_FEATURES_MATCH;
		while (rand_index_4 == rand_index_1 || rand_index_4 == rand_index_2 || rand_index_4 == rand_index_3) {
			rand_index_4 = rand() % NUM_FEATURES_MATCH;
		}

		std::vector<int> inputX;
		std::vector<int> inputY;
		std::vector<int> outputX;
		std::vector<int> outputY;

		inputX.push_back(topFeatures[rand_index_1].posx);
		inputX.push_back(topFeatures[rand_index_2].posx);
		inputX.push_back(topFeatures[rand_index_3].posx);
		inputX.push_back(topFeatures[rand_index_4].posx);

		outputX.push_back(matchingFeatures[rand_index_1].posx);
		outputX.push_back(matchingFeatures[rand_index_2].posx);
		outputX.push_back(matchingFeatures[rand_index_3].posx);
		outputX.push_back(matchingFeatures[rand_index_4].posx);

		inputY.push_back(topFeatures[rand_index_1].posy);
		inputY.push_back(topFeatures[rand_index_2].posy);
		inputY.push_back(topFeatures[rand_index_3].posy);
		inputY.push_back(topFeatures[rand_index_4].posy);

		outputY.push_back(matchingFeatures[rand_index_1].posy);
		outputY.push_back(matchingFeatures[rand_index_2].posy);
		outputY.push_back(matchingFeatures[rand_index_3].posy);
		outputY.push_back(matchingFeatures[rand_index_4].posy);

		// // Normalize dlt
		// double c1x, c1y, c2x, c2y, s1, s2;
		// compute_centroid(&c1x, &c1y, &c2x, &c2y, &s1, &s2, inputX, inputY, outputX, outputY);
		// double **T1 = new double*[3];
		// for (int i = 0; i < 3; i++) {
		// 	T1[i] = new double[3];
		// }
		// double **T2 = new double*[3];
		// for (int i = 0; i < 3; i++) {
		// 	T2[i] = new double[3];
		// }

		// allocate_T(T1, c1x, c1y, s1);
		// allocate_T(T2, c2x, c2y, s2);

		// std::vector<int> norm_inputX;
		// std::vector<int> norm_inputY;
		// std::vector<int> norm_outputX;
		// std::vector<int> norm_outputY;
		// for (int i = 0; i < inputX.size(); i++) {
		// 	norm_inputX.push_back(s1*(inputX[i] - c1x));
		// 	norm_inputY.push_back(s1*(inputY[i] - c1y));
		// 	norm_outputX.push_back(s2*(outputX[i] - c2x));
		// 	norm_outputY.push_back(s2*(outputY[i] - c2y));

		// }

		std::vector<double> norm_h = findHomographyMatrix(inputX, inputY, outputX, outputY);

		// // Denormalize h.
		// double **T2_inv = new double*[3];
		// for (int i = 0; i < 3; i++) {
		// 	T2_inv[i] = new double[3];
		// }
		// inverse(T2, T2_inv);

		// double **norm_h_matrix = new double*[3];
		// for (int i = 0; i < 3; i++) {
		// 	norm_h_matrix[i] = new double[3];
		// 	for (int j = 0; j < 3; j++) {
		// 		norm_h_matrix[i][j] = norm_h[i * 3 + j];
		// 	}
		// }

		// double **temp = new double*[3];
		// for (int i = 0; i < 3; i++) {
		// 	temp[i] = new double[3];
		// }

		double **h = new double*[3];
		for (int i = 0; i < 3; i++) {
			h[i] = new double[3];
		}
		// matrix_multiply(T2_inv, norm_h_matrix, temp);
		// matrix_multiply(temp, T1, h);

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				h[i][j] = norm_h[i * 3 + j];
			}
		}

		for (int x = 0; x < topFeatures.size(); x++) {
			int m = matchingFeatures[x].posx;
			int n = matchingFeatures[x].posy;
			int p = topFeatures[x].posx;
			int q = topFeatures[x].posy;

			int v1 = m - p;
			int v2 = n - q;

			int u1 = ((h[0][0] * p + h[0][1] * q + h[0][2] * 1) / (h[2][0] * p + h[2][1] * q + h[2][2] * 1)) - p;
			int u2 = ((h[1][0] * p + h[1][1] * q + h[1][2] * 1) / (h[2][0] * p + h[2][1] * q + h[2][2] * 1)) - q;

			if (((u1 - v1) * (u1 - v1) + (u2 - v2) * (u2 - v2)) <= THRESHOLD) {
				num_inliner++;
			}
		}

		if (num_inliner > max_inliner) {
			max_inliner = num_inliner;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					homography_matrix[i][j] = h[i][j];
				}
			}
		}
	}


	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			std::cout << '#' << homography_matrix[i][j];
			std::cout << "\n";
		}
	}

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			int new_x = (homography_matrix[0][0] * i + homography_matrix[0][1] * j + homography_matrix[0][2]) / (homography_matrix[2][0] * i + homography_matrix[2][1] * j + homography_matrix[2][2]);
			int new_y = (homography_matrix[1][0] * i + homography_matrix[1][1] * j + homography_matrix[1][2]) / (homography_matrix[2][0] * i + homography_matrix[2][1] * j + homography_matrix[2][2]);
			if (new_x < width && new_x > -1 && new_y > -1 && new_y < height) {
				blendImage.Pixel(i, j) = Pixel(i, j) * 0.5 + otherImage->Pixel(new_x, new_y) * 0.5;
			}
		}
	}
/*
	for (int x = 0; x < topFeatures.size(); x++) {
		int m = matchingFeatures[x].posx;
		int n = matchingFeatures[x].posy;
		int p = topFeatures[x].posx;
		int q = topFeatures[x].posy;

		int v1 = m - p;
		int v2 = n - q;

		int u1 = ((homography_matrix[0] * p + homography_matrix[1] * q + homography_matrix[2] * 1)
			/ (homography_matrix[6] * p + homography_matrix[7] * q + homography_matrix[8] * 1)) - p;
		int u2 = ((homography_matrix[3] * p + homography_matrix[4] * q + homography_matrix[5] * 1)
			/ (homography_matrix[6] * p + homography_matrix[7] * q + homography_matrix[8] * 1)) - q;

		if (((u1 - v1) * (u1 - v1) + (u2 - v2) * (u2 - v2)) <= THRESHOLD) {
			markedImage.line(m, p, n, q, 0.0, 1.0, 0.0);
		}
		else {
			markedImage.line(m, p, n, q, 1.0, 0.0, 0.0);
		}
	}*/

	*this = blendImage;

	return;
}

double averageSmooth(std::vector<int> vector) {
	double sum = 0;
	for (int i = 0; i < vector.size(); i++) {
		sum += vector[i];
	}
	sum = sum / vector.size();
	return sum;
}


std::vector<double> gaussianSmooth(std::vector<int> vector, int sigma) {
	// Gaussian smooth
	static const double PI = 3.1415;
	int KERNEL_SIZE = 6 * sigma + 1;
	std::vector<double> kernel;
	double sum = 0.0;
	for (int i = 0; i < KERNEL_SIZE; i++) {
		int x = i - 3 * sigma;
		double kernel_val = 1 / sqrt(2 * PI * sigma * sigma) * exp((0 - x * x) / (2 * sigma * sigma));
		sum += kernel_val;
		kernel.push_back(kernel_val);
	}
	// Normalize the kernel.
	for (int i = 0; i < KERNEL_SIZE; i++) {
		kernel[i] = kernel[i] / sum;
	}
	std::vector<double> tempVector;
	tempVector.reserve(vector.size());

	for (int i = 0; i < vector.size(); i++) {
		tempVector.push_back(0);
		for (int x = 0 - 3 * sigma; x <= 3 * sigma; x++) {
			tempVector[i] += vector[bound(i + x, vector.size())] * kernel[x + 3 * sigma];
		}
	}

	return tempVector;
}

void R2Image::
videoStabilization(int frame_num, char* input_folder, char* output_folder)
{	
	const char* jpg = ".jpg";
	const char* zero = "0";
	const int NUM_FEATURES_TRACKED = 100;
	const float HARRIS_SIGMA = 2.0;
	const int FILENAME_LENGTH = 5;

	// Read frames into a vector
	std::vector<R2Image*> frames;
	for (int i = 1; i <= frame_num; i++) {
		char index[100];
		snprintf(index, 100, "%d", i);
		char filename[100];
		strcpy(filename, input_folder); // copy string one into the result.
		for (int j = 0; j < FILENAME_LENGTH - strlen(index); j++) {
			strcat(filename, zero);
		}
		strcat(filename, index); // append string two to the result.
		strcat(filename, jpg);
		printf("Reading %d images...\r", i);
		R2Image* other_image = new R2Image(filename);
		frames.push_back(other_image);
	}
	printf("\n");

	// Video stabilization

	R2Image harrisImage(*frames[0]);

	std::vector<HarrisPixel> feats;
	
	std::vector<int> x_trans_seq_total;
	std::vector<int> y_trans_seq_total;

	// Store each continuous feature tracking sequence to further smoothing.
	// Reset when a new set of Harris features are found. 
	// std::vector<int> x_trans_seq_tmp;
	// std::vector<int> y_trans_seq_tmp;

	// std::vector<int> smoothX_total;
	// std::vector<int> smoothY_total;

	// // Store each continuous smoothing subsequence.
	// // Reset when a new set of Harris features are found. 
	// std::vector<int> smoothX_tmp;
	// std::vector<int> smoothY_tmp;

	// Track feats onto consecutive frames
	for (int f = 0; f < frames.size()-1; f++) {
		if (feats.size() < 10) {
			// Restart a tracking sequence if there exist too few features in the current frame to track.

			// Record previous sub tracking sequence.

			// // Smooth the previous tracking sequence.
			// printf("Smoothing motion sequence...\n");				
			// std::vector<double> smoothX_tmp = gaussianSmooth(x_trans_seq_tmp, 10);
			// std::vector<double> smoothY_tmp = gaussianSmooth(y_trans_seq_tmp, 10);
			// Record the subsequence.
			// smoothX_total.reserve(smoothX_total.size() + smoothX_tmp.size());
			// smoothX_total.insert(smoothX_total.end(), smoothX_tmp.begin(), smoothX_tmp.end());

			// smoothY_total.reserve(smoothY_total.size() + smoothY_tmp.size());
			// smoothY_total.insert(smoothY_total.end(), smoothY_tmp.begin(), smoothY_tmp.end());

			// x_trans_seq_total.reserve(x_trans_seq_total.size() + x_trans_seq_tmp.size());
			// x_trans_seq_total.insert(x_trans_seq_total.end(), x_trans_seq_tmp.begin(), x_trans_seq_tmp.end());

			// y_trans_seq_total.reserve(y_trans_seq_total.size() + y_trans_seq_tmp.size());
			// y_trans_seq_total.insert(y_trans_seq_total.end(), y_trans_seq_tmp.begin(), y_trans_seq_tmp.end());
			// // Reset tmp sequence.
			// // smoothX_tmp.clear();
			// // smoothY_tmp.clear();
			// x_trans_seq_tmp.clear();
			// y_trans_seq_tmp.clear();

			//double x_avg = averageSmooth(x_trans_seq_tmp);
			//double y_avg = averageSmooth(y_trans_seq_tmp);

			// Start a new tracking sequence

			// Initialize the start of each continuous tracking sequence at the end of previous sequence
			// (or <0,0> if the first frame).
			if (x_trans_seq_total.size() == 0) {
				// x_trans_seq_tmp.push_back(0);
				// y_trans_seq_tmp.push_back(0);
				x_trans_seq_total.push_back(0);
				y_trans_seq_total.push_back(0);
			} else {
				// x_trans_seq_tmp.push_back(x_trans_seq_total.back());
				// y_trans_seq_tmp.push_back(y_trans_seq_total.back());
				x_trans_seq_total.push_back(x_trans_seq_total.back());
				y_trans_seq_total.push_back(y_trans_seq_total.back());
			}
			
			R2Image new_harrisImage(*frames[f]);
			printf("Finding %d features on the %dth frame...\r",NUM_FEATURES_TRACKED, f+1);
			feats = new_harrisImage.Harris(HARRIS_SIGMA, NUM_FEATURES_TRACKED);
			printf("%d features found                        \n\n", NUM_FEATURES_TRACKED);
		
		} else {
			// Else track features from the previous frame on current frame.

			// Feature tracking on the current frame.
			std::vector<HarrisPixel> current_feats;
			current_feats.reserve(feats.size());

			printf("Tracking %lu features...\r", feats.size());
			for (int feat_index = 0; feat_index < feats.size(); feat_index++) {
				HarrisPixel matchPix = frames[f]->Search(*frames[f], *frames[f + 1], feats[feat_index], 0.1 /** search_win_perc **/);
				current_feats.push_back(matchPix);
			}
			printf("%lu features matched.			\r", current_feats.size());

			//RANSAC elimination
			std::vector<HarrisPixel> feats_accepted;
			feats_accepted.reserve(current_feats.size());

			// Try N random vectors as reference.
			const int N = 50;
			// using sqr of length
			const float THRESHOLD = 9;
			int max_inliner = 0;
			int track_vector_x = 0;
			int track_vector_y = 0;

			for (int n = 0; n < N; n++) {
				int num_inliner = 0;
				int rand_index = rand() % current_feats.size();

				std::vector<int> temp_accept;
				temp_accept.reserve(current_feats.size());

				int u1 = current_feats[rand_index].posx - feats[rand_index].posx;
				int u2 = current_feats[rand_index].posy - feats[rand_index].posy;

				for (int x = 0; x < feats.size(); x++) {
					int m = current_feats[x].posx;
					int n = current_feats[x].posy;
					int p = feats[x].posx;
					int q = feats[x].posy;

					int v1 = m - p;
					int v2 = n - q;

					if (((u1 - v1) * (u1 - v1) + (u2 - v2) * (u2 - v2)) <= THRESHOLD) {
						num_inliner++;
						temp_accept.push_back(x);
					}
				}

				if (num_inliner > max_inliner) {
					max_inliner = num_inliner;
					track_vector_x = u1;
					track_vector_y = u2;
					feats_accepted.clear();
					for (int tmp_ind= 0; tmp_ind < temp_accept.size(); tmp_ind++) {
  						feats_accepted.push_back(current_feats[temp_accept[tmp_ind]]);
					}
				}
			}

			feats = feats_accepted;
			printf("%lu features accepted. \n", feats.size());

			// Get translation vector
			printf("Calculating %d translation vector: ", f+1);
			printf("%d, ", track_vector_x);
			printf("%d\n\n", track_vector_y);
			// x_trans_seq_tmp.push_back(track_vector_x + x_trans_seq_tmp[x_trans_seq_tmp.size() - 1]);
			// y_trans_seq_tmp.push_back(track_vector_y + y_trans_seq_tmp[y_trans_seq_tmp.size() - 1]);
			x_trans_seq_total.push_back(track_vector_x + x_trans_seq_total.back());
			y_trans_seq_total.push_back(track_vector_y + y_trans_seq_total.back());
		}
	}


	printf("Smoothing motion sequence...\n");				
	std::vector<double> smoothX_total = gaussianSmooth(x_trans_seq_total, 10);
	std::vector<double> smoothY_total = gaussianSmooth(y_trans_seq_total, 10);

	for (int f_ind = 0; f_ind < frames.size(); f_ind++) {
		if (f_ind % 50 == 0) {
			printf("Stabilizing frame %d...\r", f_ind);
		}
		R2Image stabilized_frame(frames[f_ind]->width, frames[f_ind]->height);

		double x_stab = x_trans_seq_total[f_ind] - smoothX_total[f_ind];
		double y_stab = y_trans_seq_total[f_ind] - smoothY_total[f_ind];

		// Apply the homography matrix to every pixel in each frame.
		for (int x = 0; x < frames[f_ind]->width; x++) {
			for (int y = 0; y < frames[f_ind]->height; y++) {
				int new_x = round(x + x_stab);
				int new_y = round(y + y_stab);
				if (new_x >= 0 && new_x < frames[f_ind]->width && new_y >= 0 && new_y < frames[f_ind]->height) {
					stabilized_frame.Pixel(x, y) = frames[f_ind]->Pixel(new_x, new_y);
				}
				else {
					stabilized_frame.Pixel(x, y) = R2Pixel(0, 0, 0, 1);
				}
			}
		}


		// Write the output jpg series.
		if (f_ind % 50 == 0) {
			printf("Writing out frame %d...\r", f_ind);
		}
		char f_char[100];
		snprintf(f_char, 100, "%d", f_ind + 1);
		char file_out[100];
		strcpy(file_out, output_folder);
		for (int j = 0; j < FILENAME_LENGTH - strlen(f_char); j++) {
			strcat(file_out, zero);
		}
		strcat(file_out, f_char);
		strcat(file_out, jpg);

		stabilized_frame.Write(file_out);
	}
	printf("\n");
	printf("Successfully stabilized %d frames.\n", frame_num);
	return;
}

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp);
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s \n", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);

  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);

  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);

  // Check info header
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }

  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }

  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
	// Read pixel values
	int red, green, blue;
	if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
	  fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
	  fclose(fp);
	  return 0;
	}

	// Assign pixel values
	double r = (double) red / max_value;
	double g = (double) green / max_value;
	double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }

    // Close file
    fclose(fp);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
  extern "C" {
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}




int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width; 	/* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 95, TRUE);
  jpeg_start_compress(&cinfo, TRUE);

  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}





