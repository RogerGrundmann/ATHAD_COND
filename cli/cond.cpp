#include <iostream>
#include <stdlib.h>

#include "cAtmosphereModel.h"

int main(int argc, char **argv){

    cAtmosphereModel model;

    if(argc != 2){
        std::cout << std::endl << "ATHAD — Atmosphere of the Earth in the Hadean Eon" << std::endl;
        std::cout << std::endl;
        std::cout << "Invalid Command Line Parameter" << std::endl;
        std::cout << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "\t" << "./cond <<XML configuration file path>>" << std::endl;
        std::cout << "\t" << "For example: ./cond config_cond.xml" << std::endl;
        std::cout << std::endl;
        exit(1);
    }

    try{
        model.LoadConfig(argv[1]);                                      // loading config_cond.xml located at /ATHAD/cli and /ATHAD/python
        model.Run();                                                    // located in cAtmosphereModel.cpp as Run()
    }

    catch(const std::exception &exc){
        std::cerr << exc.what() << std::endl;
    }
}
