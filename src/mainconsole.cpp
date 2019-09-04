#define OUTPUT
#define OUTPUT2
//#define DETERMINIST
#include <mecacell/mecacell.h>
#include "core/GingiScenario.hpp"
#include <omp.h>



int main(int argc, char *argv[]) {
    //omp_set_num_threads(4);
    long valMS = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch()).count();
    std::stringstream dirname;
    dirname << "../Simulations/" << valMS;
    mkdir(dirname.str().c_str(), 0777);
    using scenario_t = GingiScenario;
    for (int AS = 8; AS>=3; AS--) {
        for (int AI = 8; AI>=3; AI--) {
            for (int DS = 5; DS>=2; DS--) {
                #pragma omp parallel for
                for (int DI = 5; DI>=2; DI--) {
                    stringstream path;
                    //std::ofstream myfile;
                    //std::cout << "DI=" << DI << std::endl;
                    path << dirname.str() << "/data" << AS <<"_"<< AI <<"_"<< DS <<"_"<< DI << ".csv"; // nouveau nom
                    //myfile.open(path, std::ofstream::app);
                    scenario_t scenar;
                    scenar.init(pow(10.,AS*-1.),pow(10.,AI*-1.),pow(10.,DS*-1.),pow(10.,DI*-1.),path.str());
                    int i = 0;
                    while (!scenar.stop()) {
                        scenar.loop();
                        i++;
                    }
                    scenar.myfile.close();
                }
            }
        }
    }

        

}
