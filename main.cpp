// Identities degradation by motion blur in PPM image files
//
// The program cuts a 512x512 piece from the images,
// Fourier transforms it and finds line structures in
// transformed that correspond to motion blur.
//
// Sources: 
//
// Author: Nikolaus Mayer
///////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <math.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include <CTensor.h>
#include <CFilter.h>
using namespace std;

#ifndef PI
#define PI 3.1415926535
#endif

int main(int argc, char **args) {
    int mode = 0;

    if (argc < 2)
    {
        cout << "Identification of images degraded by motion blur" << endl;
        cout << "Usage: ./motionblur {findlines, sortout} (if sortout: image1.ppm [image2.ppm image3.ppm ... imageN.ppm])" << endl;
        return 1;
    }
    if (argc < 3)
    {
        cout << "No images given!" << endl;
        return 1;
    }

    if (strcmp(args[1], "findlines") == 0)
        mode = 1;
    else if (strcmp(args[1], "sortout") == 0)
        mode = 2;
    else
    {
        cerr << "Error: First argument must be one of {findlines, sortout}." << endl;
        return 1;
    }



    /// "preprocessing": Canny filtering images to find strong edges
    if (mode == 1)
    {        
        vector<char*> filenames;
        for (int i = 2; i < argc; i++)
            filenames.push_back(args[i]);
        
        for (vector<char*>::iterator iter = filenames.begin(); iter != filenames.end(); ++iter)
        {
            string s = (*iter);
            char *charbuf = new char[255];
            cout << "File: " << s << endl;

            CTensor<double> in_img;
            in_img.readFromPPM(s.c_str());
            int width = in_img.xSize();
            int height = in_img.ySize();

            in_img.downsample(width, height);
            CMatrix<double> in_layerx, in_layery, tmp, mag;
            /// all color layers are equally blurred (?)
            in_layerx = in_img.getMatrix(0);
            in_layery = in_img.getMatrix(0);
            tmp.setSize(width, height);
            mag.setSize(width, height);
            
            //~ in_layer.downsample(512, 512);

            //~ in_img.downsample(512, 512);
            //~ in_img.fft();
            //~ in_img.writeToPPM("fffuuuuu.ppm");
            //~ in_img.normalize(0.0, 255.0);
            //~ in_layer = in_img.getMatrix(0);
            //~ in_layer.writeToPGM("ffu.pgm");
            //~ return 0;


            //~ NFilter::filter(in_layerx, CSmooth<double>::CSmooth(2, 1), 1);
            //~ NFilter::filter(in_layery, CSmooth<double>::CSmooth(2, 1), 1);
            
            /// (debug) write 
            //~ sprintf(charbuf, "resized_blackandwhite_%s.pgm", s.c_str());
            //~ in_layerx.writeToPGM(charbuf);
            
            NFilter::filter(in_layerx, CDerivative<double>::CDerivative(3), 1);
            NFilter::filter(in_layery, 1, CDerivative<double>::CDerivative(3));

            
            for (int y = 0; y < height; y++)
                for (int x = 0; x < width; x++)
                {
                    /// gradient direction
                    if (in_layery(x, y) == 0 && in_layerx(x, y) == 0)
                        tmp(x, y) = 0.0;
                    else
                        tmp(x, y) = atan2(in_layery(x, y), in_layerx(x, y));
                    /// gradient magnitude
                    mag(x, y) = sqrt(in_layerx(x, y)*in_layerx(x, y) + in_layery(x, y)*in_layery(x, y));
                }

            /// (debug) write gradient magnitude image
            mag.normalize(0., 255.);
            //~ sprintf(charbuf, "gradient_magnitude_%s.pgm", s.c_str());
            //~ mag.writeToPGM(charbuf);

            /// non-maximum suppression
            for (int y = 0; y < height; y++)
                for (int x = 0; x < width; x++)
                {
                    if (x == 0 || y == 0 || x == width-1 || y == height-1)
                        { tmp(x, y) = 0; }
                    else if (tmp(x, y) >= 0.875*PI || tmp(x, y) < -0.875*PI || (tmp(x, y) >= -0.125*PI && tmp(x, y) < 0.125*PI))
                    {
                        if (mag(x, y) <= mag(x-1, y) || mag(x, y) <= mag(x+1, y))
                            mag(x, y) = 0;
                    }
                    else if ((tmp(x, y) >= -0.875*PI && tmp(x, y) < -0.625*PI) || (tmp(x, y) >= 0.125*PI && tmp(x, y) < 0.375*PI))
                    {
                        if (mag(x, y) <= mag(x-1, y-1) || mag(x, y) <= mag(x+1, y+1))
                            mag(x, y) = 0;
                    }
                    else if ((tmp(x, y) >= -0.625*PI && tmp(x, y) < -0.375*PI) || (tmp(x, y) >= 0.375*PI && tmp(x, y) < 0.625*PI))
                    {
                        if (mag(x, y) <= mag(x, y+1) || mag(x, y) <= mag(x, y-1))
                            mag(x, y) = 0;
                    }
                    else if ((tmp(x, y) >= -0.375*PI && tmp(x, y) < -0.125*PI) || (tmp(x, y) >= 0.625*PI && tmp(x, y) < 0.875*PI))
                    {
                        if (mag(x, y) <= mag(x+1, y-1) || mag(x, y) <= mag(x-1, y+1))
                            mag(x, y) = 0;
                    }
                    else
                        cout << "something's amiss! x=" << x << ", y=" << y << ", tmp(x, y)=" << tmp(x, y) << endl;
                    
                }


            /// thresholding
            mag.clip(70, 255);

            /// (debug) write lines image
            mag.normalize(0., 255.);
            sprintf(charbuf, "%s_Canny.pgm", s.c_str());
            mag.writeToPGM(charbuf);

            delete[] charbuf;
        }
    }

    /// blur estimation
    else if (mode == 2)
    {
        vector<float> scores;
        vector<char*> filenames;
        vector<char*> ok_files;
        char *charbuf = new char[255];

        /// try to read existing scores file
        ifstream infile("scores.txt");
        if (infile.fail())
        {
            cerr << "Could not read file scores.txt!" << endl;
        }
        else
        {
            int s;
            infile >> s;
            for (int i = 0; i < s; ++i)
            {
                float sc;
                char* f = new char[255];
                infile >> f;
                filenames.push_back(f);
                infile >> sc;
                scores.push_back((float)sc);
            }
            infile.close();
        }

        /// if the scores file could not be found/read, create scores and write the file
        if (scores.size() == 0)
        {
            for (int i = 2; i < argc; i++)
                filenames.push_back(args[i]);

            /// compute scores
            for (vector<char*>::iterator iter = filenames.begin(); iter != filenames.end(); ++iter)
            {
                sprintf(charbuf, "./Canny/%s_Canny.pgm", *iter);
                
                cout << "Reading " << charbuf << " for " << *iter << endl;

                CMatrix<float> img;
                img.readFromPGM(charbuf);

                float score = 0.;
                //~ for (int x = 0; x < img.xSize(); ++x)
                    //~ for (int y = 0; y < img.ySize(); ++y)
                        //~ score += img(x, y);
                for (int x = 0.25*img.xSize(); x < 0.75*img.xSize(); ++x)
                    for (int y = 0.25*img.xSize(); y < 0.75*img.ySize(); ++y)
                        score += img(x, y);
                        
                scores.push_back(score);

            }

            /// debug: write scores to file
            ofstream outfile ("scores.txt");
            outfile << scores.size() << endl;
            for (unsigned int i = 0; i < scores.size(); ++i)
                /// (long) casting to avoid float's scientific number notation (1e+06)
                outfile << filenames[i] << "\t" << (long)scores[i] << endl;
            outfile.close();
        }
            

        /// rate images
        int neighborhood_size = 10;
        for (unsigned int i = 0; i < scores.size(); ++i)
        {
            float mean_score = 0.;
            
            for (int shift = -neighborhood_size; shift < neighborhood_size; shift++)
            {
                if (shift == 0 || i+shift < 0 || i+shift >= scores.size())
                    continue;
                else
                    mean_score += scores[i+shift]/(2*neighborhood_size);
            }

            if (scores[i] < 0.85 * mean_score)
                cout << "Image " << filenames[i] << " seems to be blurry." << endl;
            else
                ok_files.push_back(filenames[i]);
        }

        /// save good file names
        ofstream outfile ("ok.txt");
        outfile << ok_files.size() << endl;
        for (unsigned int i = 0; i < ok_files.size(); ++i)
            outfile << ok_files[i] << endl;
        outfile.close();            
        cout << filenames.size() - ok_files.size() << " of " << filenames.size() << " images have been dismissed. OK images saved to 'ok.txt'." << endl;

        delete[] charbuf;
    }

    
    //~ in_layerx.writeToPGM("x.pgm");
    //~ in_layer.normalize(-255.0, 255.0);
    //~ in_layer.writeToPGM("deriv.pgm");
    
    //~ double *in_data = in_layer.data();
    //~ int width  = in_layer.xSize();
    //~ int height = in_layer.ySize();
    //~ 
    //~ fftw_plan p;
    //~ fftw_complex *in;
    //~ fftw_complex *out;
    //~ 
    //~ in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width*height);
    //~ out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width*height);
    //~ 
    //~ p = fftw_plan_dft_2d(width, height, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    //~ 
    //~ for (int i = 0; i < width*height; i++)
        //~ in[i][0] = in_data[i];
    //~ 
    //~ fftw_execute(p);
    //~ 
    //~ CMatrix<double> out_img;
    //~ out_img.setSize(width, height);
    //~ for (int y = 0; y < height; y++)
        //~ for (int x = 0; x < width; x++)
        //~ {
            //~ double r = out[y*width+x][0] / (width*height);
            //~ double i = out[y*width+x][1] / (width*height);
            //~ 
            //~ if (x < width/2 && y < height/2)
                //~ out_img(width/2 + x, height/2 + y) = log(sqrt(r*r+i*i));
            //~ else if (x >= width/2 && y < height/2)
                //~ out_img(x - width/2, height/2 + y) = log(sqrt(r*r+i*i));
            //~ else if (x < width/2 && y >= height/2)
                //~ out_img(width/2 + x, y - height/2) = log(sqrt(r*r+i*i));
            //~ else if (x >= width/2 && y >= height/2)
                //~ out_img(x - width/2, y - height/2) = log(sqrt(r*r+i*i));
            //~ else
                //~ cout << "Error: No place for x=" << x << ", y=" << y << endl;
        //~ }
    //~ 
    //~ out_img.normalize(1.0, 255.0);
    //~ out_img.writeToPGM("fft.pgm");
    //~ 
    //~ fftw_destroy_plan(p);
    //~ fftw_free(in);
    //~ fftw_free(out);
    //~ 
    //~ return 0;

}

