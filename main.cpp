#include <SFML/Graphics.hpp>
#include <iostream>
#include <string>
#include <memory>
#include <cmath>
#include "jpeg.h"
#include "image.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <jpg file>\n";
        return 0;
    }
    
    int width, height, zoom;
    std::unique_ptr<sf::Uint8> pixels;
    {
      sf::VideoMode desktop = sf::VideoMode::getDesktopMode();
      // exclude 69 pixels for titlebar, taskbar, border
      // 29 + 36 + 4
      desktop.height -= 69;
      image_t image(argv[1], desktop.width, desktop.height);
      float scale = fmax(image.aWidth  / float(desktop.width ),
                         image.aHeight / float(desktop.height));
      if (scale <= 1) {
        width  = image.aWidth;
        height = image.aHeight;
      } else {
        width  = image.aWidth / scale;
        height = image.aHeight/ scale;
      }
      pixels = std::unique_ptr<sf::Uint8>(
        rescale(image.pixels,
                image.aWidth,
                image.aHeight,
                width,
                height
      ));
      zoom = width / float(image.width) * 100;
    }
    
    std::string title = std::string(argv[1]) + " - ";
    title += std::to_string(zoom) + "% -- JPEG Viewer";
    sf::RenderWindow window(
        sf::VideoMode(width, height),
        title,
        sf::Style::Titlebar | sf::Style::Close
    );
    sf::Texture texture;
    texture.create(width, height);
    texture.update(pixels.get());
    sf::Sprite sprite;
    sprite.setTexture(texture);
    
    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear();
        window.draw(sprite);
        window.display();
    }
}
