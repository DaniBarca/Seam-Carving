#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "image.hpp"

Color myfunc(Color c) { return c * 10; }

int main(int argc,char **argv)
{
	std::cout << "Image processing\n********************" << std::endl;
	if(argc < 3)
	{
		std::cerr << "Parameter missing, usage:\n  app image.tga image2.tga" << std::endl;
		return 1;
	}
    
	const char* input_image = argv[1];
	const char* output_image = argv[2];
    
	Image img;
	if( img.loadTGA(input_image) == false )
	{
		std::cerr << "Error: image not found: " << input_image << std::endl;
		return 1;
	}
    
    Image *img2 = new Image(img.width, img.height);
    
	std::cout << "Image info:\n + Width: " << img.width << "\n + Height: " << img.height << std::endl;
    
    std::cout << "Introduzca número de píxeles a reducir: " << std::endl;
    std::string numPix;
    getline(std::cin, numPix);
    
    std::cout << "Introduzca \"VERTICAL\" o \"HORIZONTAL\" según como quiera reducir la imagen: " << std::endl;
    std::string dir;
    getline(std::cin, dir);
    
    OSemmer st = OSemmer();
    
    img2 = st.seammer(&img, dir, atoi(numPix.c_str()));

    //Por si queremos ver la gradiente y/o el acumulado
    //st.cumulative_.norm();
    st.gradient_.norm();
    st.gradient_.roll();
    st.gradient_.roll();
    st.gradient_.roll();
    
    st.gradient_.saveTGA("/Users/danibarca/Desktop/gradient.tga");
    st.cumulative_.saveTGA("/Users/danibarca/Desktop/cumulative.tga");

    img2->saveTGA(output_image);

	std::cout << "Done" << std::endl;
	return 0;
}