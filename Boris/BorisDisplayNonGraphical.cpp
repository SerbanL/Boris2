#include "stdafx.h"
#include "BorisDisplayNonGraphical.h"

#if GRAPHICS == 0

//from image file extract a raw bitmap (BYTE array with 4 BYTE-sized entries as B-G-R-A for each pixel) to specified pixel size (so image file is rescaled to specified size)
void BorisGraphics::GetBitmapFromImage(std::string fileName, std::vector<unsigned char>& bitmap, INT2 n_plane)
{
	//load texture from file
	sf::Texture texture;
	if (!texture.loadFromFile(fileName))
	{
		std::cout << "Could not load mask file" << std::endl;
		return;
	}

	//copy to cpu memory as rgba data
	sf::Image image = texture.copyToImage();

	//original image dimensions
	INT2 n_orig = INT2(image.getSize().x, image.getSize().y);

	//map original bitmap to required n_plane pixel sizes

	//make sure bitmap has the right size
	bitmap.resize(n_plane.dim() * 4);

	for (int i = 0; i < n_plane.x; i++) {
		for (int j = 0; j < n_plane.y; j++) {

			int io = floor(((double)i / n_plane.x) * n_orig.x);
			int jo = floor(((double)j / n_plane.y) * n_orig.y);

			//get bgra data
			sf::Color color = image.getPixel(io, jo);

			bitmap[(i + j * n_plane.x) * 4] = (unsigned char)color.b;
			bitmap[(i + j * n_plane.x) * 4 + 1] = (unsigned char)color.g;
			bitmap[(i + j * n_plane.x) * 4 + 2] = (unsigned char)color.r;
			bitmap[(i + j * n_plane.x) * 4 + 3] = (unsigned char)color.a;
		}
	}
}

void BorisDisplay::DisplayFormattedConsoleMessage(std::string text)
{
	//first trim all formmatting specifiers in block form (i.e. with extra text info between start and end text formatting specifiers)
	for (int idx = 0; idx < TF::textblockSpecifiers.size(); idx++) {

		text = trimblock(text, TF::textblockSpecifiers[idx].first, TF::textblockSpecifiers[idx].second);
	}

	//now trim all other formatting specifiers
	for (int idx = 0; idx < TF::simpleFormattingSpecifiers.size(); idx++) {

		text = trim(text, TF::simpleFormattingSpecifiers[idx]);
	}

	//finally display the unformatted text
	std::cout << text << std::endl;
}

#endif