/*** image.h  Javi Agenjo (javi.agenjo@gmail.com) 2013
	This file defines the class Image that allows to manipulate images.
	It defines all the need operators for Color and Image.
	It has a TGA loader and saver.
***/

#ifndef IMAGE_H
#define IMAGE_H

#include <string.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

float clamp(float x, float a, float b) { return x < a ? a : (x > b ? b : x); }

#pragma warning(disable : C4996)

class Color
{
public:
	union {
		struct { float r; float g; float b; };
		float v[3];
	};

	Color() { r = g = b = 0.0; }
	Color(float r, float g, float b) { this->r = r; this->g = g; this->b = b; }
	void set(float r, float g, float b) { this->r = r; this->g = g; this->b = b; }
    Color abs(){ r = fabs(r); g = fabs(g); b = fabs(b); return *this;}
};

Color operator * (const Color& c,float v) { return Color(c.r*v, c.g*v, c.b*v); }
void operator *= (Color& c,float v) { c = c * v; }
Color operator / (const Color& c,float v) { return Color(c.r/v, c.g/v, c.b/v); }
void operator /= (Color& c,float v) { c = c / v; }
Color operator * (float v, const Color& c) { return Color(c.r*v, c.g*v, c.b*v); }
Color operator + (const Color& a, const Color& b) { return Color(a.r+b.r, a.g+b.g, a.b+b.b); }
void operator += (Color& a,const Color& b) { a = a + b; }
Color operator - (const Color& a, const Color& b) { return Color(a.r-b.r, a.g-b.g, a.b-b.b); }
void operator -= (Color& a,const Color& b) { a = a - b; }
bool operator  < (const Color& a, const int& b) { return a.r+a.g+a.b < b;}

class Image
{
	//a general struct to store all the information about a TGA file
	typedef struct sTGAInfo 
	{
		unsigned int width;
		unsigned int height;
		unsigned int bpp; //bits per pixel
		unsigned char* data; //bytes with the pixel information
	} TGAInfo;

public:
	unsigned int width;
	unsigned int height;
	Color* pixels;

	/* CONSTRUCTORS */
	Image() {
		width = 0; height = 0;
		pixels = NULL;
	}

	Image(unsigned int width, unsigned int height)
	{
		this->width = width;
		this->height = height;
		pixels = new Color[width*height];
		memset(pixels, 0, width * height * sizeof(Color));
	}

	//copy constructor
	Image(const Image& c) {
		if(pixels) delete pixels;
		pixels = NULL;

		width = c.width;
		height = c.height;
		if(c.pixels)
		{
			pixels = new Color[width*height];
			memcpy(pixels, c.pixels, width*height*sizeof(Color));
		}
	}

	//assign operator
	Image& operator = (const Image& c)
	{
		if(pixels) delete pixels;
		pixels = NULL;

		width = c.width;
		height = c.height;
		if(c.pixels)
		{
			pixels = new Color[width*height*sizeof(Color)];
			memcpy(pixels, c.pixels, width*height*sizeof(Color));
		}
		return *this;
	}
    
    void norm(){
        float max = 0;
        for(int i = 0; i < width; ++i){
            for(int j = 0; j < height; ++j){
                if(getPixel(i,j).r > max){
                    max = getPixel(i, j).r;
                }
            }
        }
        
        for(int i = 0; i < width; ++i){
            for(int j = 0; j < height; ++j){
                setPixel(getPixel(i,j)/max, i, j);
            }
        }
    }

	~Image()
	{
		if(pixels) delete pixels;
	}
    
    /*
     * Gira la imagen en el sentido de las agujas del reloj
     */
    void roll(){
        Image aux  = Image(height,width);
        
        for(int i = 0; i < width; ++i){
            for(int j = 0; j < height; ++j){
                aux.setPixel(getPixel(i,j), height-1-j, i);
            }
        }
        
        resize(aux.width,aux.height);
        
        for(int i = 0; i < aux.width; ++i){
            for(int j = 0; j < aux.height; ++j){
                setPixel(aux.getPixel(i,j), i, j);
            }
        }
    }
    
	//get the pixel at position x,y
	Color getPixel(unsigned int x, unsigned int y) const
	{        
		return pixels[y * width + x];
	}

	//set the pixel at position x,y with value C
	void setPixel(const Color& c, unsigned int x, unsigned int y)
	{
		pixels[ y * width + x ] = c;
	}

	//change image size (the old one will remain in the top-left corner)
	void resize(unsigned int width, unsigned int height)
	{
		Color* new_pixels = new Color[width*height];
		unsigned int min_width = this->width > width ? width : this->width;
		unsigned int min_height = this->height > height ? height : this->height;

		for(unsigned int x = 0; x < min_width; ++x)
			for(unsigned int y = 0; y < min_height; ++y)
				new_pixels[ y * width + x ] = getPixel(x,y);

		delete pixels;
		this->width = width;
		this->height = height;
		pixels = new_pixels;
	}

	//flip the image top-down
	void flipY()
	{
		for(unsigned int x = 0; x < width; ++x)
			for(unsigned int y = 0; y < height * 0.5; ++y)
			{
				Color temp = getPixel(x, height - y - 1);
				setPixel( getPixel(x,y), x, height - y - 1);
				setPixel( temp, x, y );
			}
	}

	//flip the image left-right
	void flipX()
	{
		for(unsigned int x = 0; x < width * 0.5; ++x)
			for(unsigned int y = 0; y < height; ++y)
			{
				Color temp = getPixel(width - x - 1, y);
				setPixel( getPixel(x,y), width - x - 1, y);
				setPixel( temp, x, y );
			}
	}

	//fill the image with the color C
	void fill(const Color& c)
	{
		for(unsigned int pos = 0; pos < width*height; ++pos)
			pixels[pos] = c;
	}

	//returns a new image with the area from (startx,starty) of size width,height
	Image getArea(unsigned int start_x, unsigned int start_y, unsigned int width, unsigned int height)
	{
		Image result(width, height);
		for(unsigned int x = 0; x < width; ++x)
			for(unsigned int y = 0; y < height; ++x)
			{
				if( (x + start_x) < this->width && (y + start_y) < this->height) 
					result.setPixel( getPixel(x + start_x,y + start_y), x, y);
			}
		return result;
	}

	#ifndef IGNORE_LAMBDAS

	//applies an algorithm to every pixel in an image
	// you can use lambda sintax:   img.forEachPixel( [](Color c) { return c*2; });
	// or callback sintax:   img.forEachPixel( mycallback ); //the callback has to be Color mycallback(Color c) { ... }
	template <typename F>
	Image& forEachPixel( F callback )
	{
		for(unsigned int pos = 0; pos < width*height; ++pos)
			pixels[pos] = callback(pixels[pos]);
		return *this;
	}

	#endif

	//Loads an image from a TGA file
	bool loadTGA(const char* filename)
	{
		unsigned char TGAheader[12] = {0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		unsigned char TGAcompare[12];
		unsigned char header[6];
		unsigned int bytesPerPixel;
		unsigned int imageSize;

		FILE * file = fopen(filename, "rb");
   		if ( file == NULL || fread(TGAcompare, 1, sizeof(TGAcompare), file) != sizeof(TGAcompare) ||
			memcmp(TGAheader, TGAcompare, sizeof(TGAheader)) != 0 ||
			fread(header, 1, sizeof(header), file) != sizeof(header))
		{
			std::cerr << "File not found: " << filename << std::endl;
			if (file == NULL)
				return NULL;
			else
			{
				fclose(file);
				return NULL;
			}
		}

		TGAInfo* tgainfo = new TGAInfo;
    
		tgainfo->width = header[1] * 256 + header[0];
		tgainfo->height = header[3] * 256 + header[2];
    
		if (tgainfo->width <= 0 || tgainfo->height <= 0 || (header[4] != 24 && header[4] != 32))
		{
			fclose(file);
			delete tgainfo;
			return NULL;
		}
    
		tgainfo->bpp = header[4];
		bytesPerPixel = tgainfo->bpp / 8;
		imageSize = tgainfo->width * tgainfo->height * bytesPerPixel;
    
		tgainfo->data = new unsigned char[imageSize];
    
		if (tgainfo->data == NULL || fread(tgainfo->data, 1, imageSize, file) != imageSize)
		{
			if (tgainfo->data != NULL)
				delete tgainfo->data;
            
			fclose(file);
			delete tgainfo;
			return false;
		}

		fclose(file);

		//save info in image
		if(pixels)
			delete pixels;

		width = tgainfo->width;
		height = tgainfo->height;
		pixels = new Color[width*height];

		//convert to float all pixels
		for(unsigned int y = 0; y < height; ++y)
			for(unsigned int x = 0; x < width; ++x)
			{
				unsigned int pos = y * width * bytesPerPixel + x * bytesPerPixel;
				this->setPixel(Color( tgainfo->data[pos+2] / 255.0f, tgainfo->data[pos+1] / 255.0f, tgainfo->data[pos] / 255.0f),x,height - y - 1);
			}

		delete tgainfo->data;
		delete tgainfo;

		return true;
	}

	// Saves the image to a TGA file
	bool saveTGA(const char* filename)
	{
		unsigned char TGAheader[12] = {0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0};

		FILE *file = fopen(filename, "wb");
		if ( file == NULL )
		{
			fclose(file);
			return false;
		}

		unsigned short header_short[3];
		header_short[0] = width;
		header_short[1] = height;
		unsigned char* header = (unsigned char*)header_short;
		header[4] = 24;
		header[5] = 0;

		//tgainfo->width = header[1] * 256 + header[0];
		//tgainfo->height = header[3] * 256 + header[2];

		fwrite(TGAheader, 1, sizeof(TGAheader), file);
		fwrite(header, 1, 6, file);

		//convert pixels to unsigned char
		unsigned char* bytes = new unsigned char[width*height*3];
		for(unsigned int y = 0; y < height; ++y)
			for(unsigned int x = 0; x < width; ++x)
			{
				Color c = pixels[(height-y-1)*width+x];
				unsigned int pos = (y*width+x)*3;
				bytes[pos+2] = (unsigned char)clamp(c.r * 255,0,255);
				bytes[pos+1] = (unsigned char)clamp(c.g * 255,0,255);
				bytes[pos] = (unsigned char)clamp(c.b * 255,0,255);
			}

		fwrite(bytes, 1, width*height*3, file);
		fclose(file);
		return true;
	}
};

Image getGradient(Image *img);
Image getCumulative(Image *img);

typedef struct Pixel{
    int i;
    int j;
}Pixel;

class OSemmer{public:
    Image image_;
    Image gradient_;
    Image cumulative_;
    
    typedef struct seam{
        std::vector<Pixel> pixels;
        int energy;
    }seam;
    
    std::vector<seam> seams;
    

    OSemmer();
    Image*  seammer(Image *img, std::string dir, int howMany);
    void  getSeam();
    void printStem();
    void eraseSeam(seam min);
    Image flip(Image img);
};

OSemmer::OSemmer(){}

void OSemmer::eraseSeam(seam min){
    //Ara eliminem el seam
    Image aux = Image(image_.width, image_.height-1);
    Image auxg= Image(cumulative_.width, cumulative_.height-1);
    
    for(int i = 0; i < image_.width; ++i){
        for(int j = 0; j < image_.height; ++j){
            if(j < min.pixels.at(i).j){                       //Si no hem arribat al seam
                aux. setPixel(image_.getPixel(i,j), i, j);
                auxg.setPixel(cumulative_.getPixel(i,j),i,j);
            }
            if(j >= min.pixels.at(i).j && j <image_.height-1){//Si som al seam o ens l'hem passat
                aux. setPixel(image_.getPixel(i,j+1),i,j);
                auxg.setPixel(cumulative_.getPixel(i,j+1),i,j);
            }
        }
    }
    image_    = aux;
    cumulative_ = auxg;
}


Image* OSemmer::seammer(Image *img, std::string dir, int howMany){
    image_ = *img;
    //Image tester = Image();
    
    if(dir.compare("HORIZONTAL") == 0){ //Para que el acumulativo esté bien orientado
        image_.roll();
        image_.roll();
        image_.roll();
    }
    
    //gradient_    = getGradient(&image_);
    cumulative_ = getCumulative(&image_);

    image_.roll();
    cumulative_.roll();

    //tester  = image_;
    std::string dirb = "";
    seam min;
    
    for(int t = 0; t < howMany; ++t){               //Por cada píxel que queramos quitar
        if(t%10==0)
            std::cout << "Reducidos " << t << " de " << howMany << " píxeles." << std::endl;
        
        getSeam();                                  //Obtenemos los seams
        
        //Y buscamos el de menor energía
        min = seams.at(1);
        for(int i = 0; i < seams.size(); ++i){
            if(seams.at(i).energy < min.energy){
                min  = seams.at(i);
            }
        }
        seams.clear();
        
        //Guardamos el seam en tester
        //for(int i = 0; i<min.pixels.size(); ++i){
        //    tester.setPixel(Color(0,255,0), min.pixels.at(i).i, min.pixels.at(i).j+t);
        //}
        
        //Y lo eliminamos
        eraseSeam(min);
        image_.roll();
        image_.roll();
        image_.roll();
        cumulative_ = getCumulative(&image_);
        image_.roll();
        cumulative_.roll();
    }
    
    if(dir.compare("VERTICAL") == 0 || !dir.compare("HORIZONTAL") == 0){
        image_.roll();
        image_.roll();
        image_.roll();
    }
    return &image_;
}

/*
 * Este método rellena nuestro vector de seams
 */
void OSemmer::getSeam(){
    Pixel p; p.i = 0; p.j = 0;
    Pixel min;
    
    int t    = 0;
    int tmax = 0;
    int jb   = 0;
    
    seam saux;
    
    for(int j = 0; j < cumulative_.height-1; ++j){
        jb = j;
        
        p.i = 0; p.j = j;
        saux.pixels.push_back(p);
        saux.energy = cumulative_.getPixel(p.i, p.j).r;
        
        for(int i = 0; i < cumulative_.width-1; ++i){
            if(jb == 0) t = 0;
            if(jb > 0 && jb < cumulative_.height-1) t = -1;
            if(jb >= cumulative_.height-1) tmax = 0; else tmax = 1;
            
            min.i = i+1; min.j = jb;
            for(t = t; t <= tmax; t++){
                if(cumulative_.getPixel(i+1, jb+t).r < cumulative_.getPixel(min.i, min.j).r){
                    min.i = i+1;
                    min.j = jb+t;
                }
            }
            
            jb = min.j;
            p.i = min.i;
            p.j = min.j;
            
            //image_.setPixel(Color(0,0,255), p.i, p.j); //Esta línea dibuja el seam en la imagen, modifica la imagen original así que solo debe activarse para testing
            
            saux.pixels.push_back(p);
            saux.energy += cumulative_.getPixel(p.i,p.j).r;
            
        }
        seams.push_back(saux);
        saux.pixels = std::vector<Pixel>();
        saux.energy = 0;
    }
    
    
}

Image getGradient(Image *img){
    Image imgb = Image(img->width,img->height);
    
    Color c = Color();
    int aux;
    int t = 0, tmax = 0, l = 0, lmax = 0;
    
    //Recorremos la imagen
    for(int i = 0; i < img->width; ++i){
        for(int j = 0; j < img->height; ++j){
            t    = (i == 0)            ? 0 : -1;  //Estamos en el borde superior
            tmax = (i == img->width-1) ? 0 :  1;  //Estamos en el borde inferior
            l    = (j == 0)            ? 0 : -1;  //Estamos en el borde izquierdo
            lmax = (j == img->height-1)? 0 :  1;  //Estamos en el borde derecho
            c = Color(0,0,0);
            
            //Recorremos los píxeles vecinos, los cálculos anteriores evitan que nos salgamos de la imagen en los bordes
            //No podemos modificar el valor de 't' y 'l' así que usamos m y n para iterar
            for(int m = t; m <= tmax; ++m){
                for(int n = l; n <= lmax; ++n){
                    if(!(m==0 && n == 0))         //Nos saltamos el pixel i,j, no lo sumamos
                        c += (img->getPixel(i,j) - img->getPixel(i+m, j+n)).abs();
                }
            }
            
            //Hacemos que los tres píxeles pasen a tener el mismo valor y asignamos
            aux = c.r+c.g+c.b;
            c.r = aux;
            c.g = aux;
            c.b = aux;
            
            imgb.setPixel(c,i,j);
        }
    }
    
    return imgb;
}

Image getCumulative(Image *img){
    Image gradient    = Image(img->width, img->height);
    Image cumulative  = Image(img->width, img->height);
    
    gradient = getGradient(img);
    
    //Primer fem la primera fila
    for(int i = 0; i < gradient.width; ++i){
        cumulative.setPixel(gradient.getPixel(i,0),i,0);
    }
    
    //I ara la resta
    float a,b,c;
    for(int j = 1; j < cumulative.height; ++j){
        for(int i = 0; i < cumulative.width; ++i){
            if(i>0)
                a = cumulative.getPixel(i-1,j-1).r;
            else a = 100000;
            if(i>0)
                b = cumulative.getPixel(i,j-1).r;
            else b = 100000;
            if(i>0)
                c = cumulative.getPixel(i+1,j-1).r;
            else c = 100000;
            
            if(a<b){
                if(a<c)  //minim = a
                    cumulative.setPixel(gradient.getPixel(i,j) + cumulative.getPixel(i-1,j-1), i, j);
                else     //minim = c
                    cumulative.setPixel(gradient.getPixel(i,j) + cumulative.getPixel(i+1,j-1), i, j);
            }
            else{
                if(b<c)  //minim = b
                    cumulative.setPixel(gradient.getPixel(i,j) + cumulative.getPixel(i,j-1), i, j);
                else     //minim = c
                    cumulative.setPixel(gradient.getPixel(i,j) + cumulative.getPixel(i+1,j-1), i, j);
            }
        }
    }
    
    //Fem que els valors de color estiguin entre 0 i 1
    //Si fem això el resultat és molt menys precís, és millor deixar-ho com està encara que no sigui visible (l'acumulatiu té valors > 1)
    //acumulative.norm();
    
    return cumulative;
}

#ifndef IGNORE_LAMBDAS

//you can apply and algorithm for two images and store the result in the first one
//forEachPixel( img, img2, [](Color a, Color b) { return a + b; } );
template <typename F>
void forEachPixel(Image& img, const Image& img2, F f) {
	for(unsigned int pos = 0; pos < img.width * img.height; ++pos)
		img.pixels[pos] = f( img.pixels[pos], img2.pixels[pos] );
}

#endif

#endif