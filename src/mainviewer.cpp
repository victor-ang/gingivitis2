#include <mecacell/mecacell.h>
#include <mecacellviewer/viewer.h>
#include <mecacellviewer/menu/button.hpp>
#include "core/GingiScenario.hpp"
#include "core/GingiCell.hpp"

GingiScenario sc;
bool showAliveStroma;
bool xCut;

// todo list
// - Bug deplacement apoptose => DONE
// - Recrutement circulant => DONE but not optimal until I have the diffusion grid
// - Pro-résolutif
// - Transmission de la nécrose via l'inflammation => DONE

int main(int argc, char *argv[]) {

    //using scenario_t = GingiScenario;
	using Cell = GingiCell<Bodies3D>;

    //scenario_t sc;
	sc.init();

	MecacellViewer::Viewer<GingiScenario> v(sc);
	
	showAliveStroma = true;
	xCut = false;

	v.addButton("hideStromaButton", "Controls", "hide/show stroma", std::function<void(void*, void*)>([](void* a, void* b) {
		showAliveStroma  = ! showAliveStroma;
		auto& w = sc.getWorld();
		for (auto *c : w.cells) {
			c->setVisible(	(showAliveStroma || !(c->type == Cell::Stroma && c->state == Cell::Alive)) && 
							(!xCut || c->getPosition().z() < 0));
		}
	}));

	v.addButton("zCutButton", "Controls", "Z Cut", std::function<void(void*, void*)>([](void* a, void* b) {
		xCut = ! xCut;
		auto& w = sc.getWorld();
		for (auto *c : w.cells) {
			c->setVisible(	(showAliveStroma || !(c->type == Cell::Stroma && c->state == Cell::Alive)) &&
							(!xCut || c->getPosition().z() <= 0));	
		}
	}));

	v.addButton("OneNecroButton", "Controls", "1 Necrosis", std::function<void(void*, void*)>([](void* a, void* b) {
	    auto& w = sc.getWorld();
		std::uniform_int_distribution<unsigned int> dist(0, w.cells.size());
        w.cells[dist(MecaCell::Config::globalRand())]->state = Cell::Necrosis;
	}));

	v.addButton("TenNecroButton", "Controls", "10 Necrosis", std::function<void(void*, void*)>([](void* a, void* b) {
	    auto& w = sc.getWorld();
		std::uniform_int_distribution<unsigned int> dist(0, w.cells.size());
        for (int i=0; i<10; i++) {
			w.cells[dist(MecaCell::Config::globalRand())]->state = Cell::Necrosis;
		}
	}));

	v.addButton("HundredNecroButton", "Controls", "100 Necrosis", std::function<void(void*, void*)>([](void* a, void* b) {
	    auto& w = sc.getWorld();
		std::uniform_int_distribution<unsigned int> dist(0, w.cells.size());
        for (int i=0; i<100; i++) {
			w.cells[dist(MecaCell::Config::globalRand())]->state = Cell::Necrosis;
		}
	}));

	v.addButton("clusterNecroButton", "Controls", "cluster necrosis", std::function<void(void*, void*)>([](void* a, void* b) {
	    auto& w = sc.getWorld();
		std::uniform_int_distribution<unsigned int> dist(0, w.cells.size());
		Cell *c = w.cells[dist(MecaCell::Config::globalRand())];
		c->state = Cell::Necrosis;
		for (auto* nc : c->getConnectedCells()) {
			for (auto *nnc : nc->getConnectedCells()) {
				nnc->state = Cell::Necrosis;
			}
			
		}
	}));

	return v.exec();
}
