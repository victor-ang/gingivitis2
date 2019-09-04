#define OUTPUT
#define OUTPUT2
//#define DETERMINIST
#include <mecacell/mecacell.h>
#include "core/GingiScenario.hpp"

int main(int argc, char *argv[]) {


    using scenario_t = GingiScenario;

        scenario_t scenar;
		
		scenar.init();		

		int i = 0;
        while (!scenar.stop()) {
            scenar.loop();
			i++;
        }

}
