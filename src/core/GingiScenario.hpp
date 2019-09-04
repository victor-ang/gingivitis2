#ifndef _GINGISCENARIO_HPP_
#define _GINGISCENARIO_HPP_
#include "Bodies3D.hpp"
#include "GingiCell.hpp"
#include <fstream>
#include <mecacell/mecacell.h>
#include <sstream>
#include <string>
#include <chrono>
#include <sys/stat.h>

class GingiScenario {

public:
  using Cell = GingiCell<Bodies3D>;
  using World = MecaCell::World<Cell>;

  World w;
  std::ofstream myfile;
  //std::string pathOutput;
  
  

  nlohmann::json config;

  GingiScenario() {}

  static double BOX_HALF_SIZE;

  inline World &getWorld() { return w; }

  void init(double apopStroma, double apopImmune, double divStroma, double divImmune, std::string path) {
    //std::cout << apopStroma << "_" << apopImmune << "_"<< divStroma << "_"<< divImmune << "_"<< _pathOutput << std::endl;
    myfile.open(path, std::ofstream::app);

    config = loadJsonConfig("../j.json");
    Molecule infla(400.0, 0.0, 1.0, 0.1);
    Molecule reso(400.0, 0.0, 1.0, 0.1);
    Molecule eatme(600.0, 0.0, 1.0, 0.1);

    w.cellPlugin.diffusionPlugin.addMolecule(infla);
    w.cellPlugin.diffusionPlugin.addMolecule(reso);
    w.cellPlugin.diffusionPlugin.addMolecule(eatme);

    

    w.setDt(1);
    // Add cells
    //		for (int i = -5; i <= 5; i++) {
    //			for (int j=-5; j <= 5; j++) {
    //				for (int k=-5; k <= 5; k++) {
    //					Cell* c = new Cell(MecaCell::Vec(i *
    //MecaCell::Config::DEFAULT_CELL_RADIUS * 1.75f, 													 j *
    //MecaCell::Config::DEFAULT_CELL_RADIUS * 1.75f, 													 k *
    //MecaCell::Config::DEFAULT_CELL_RADIUS * 1.75f)); 					w.addCell(c);
    //				}
    //			}
    //		}
    // std::normal_distribution<double>
    // nDist(-5.0f*MecaCell::Config::DEFAULT_CELL_RADIUS * 1.75f, 5.0f *
    // MecaCell::Config::DEFAULT_CELL_RADIUS * 1.75f);
    std::uniform_real_distribution<double> nDist(-BOX_HALF_SIZE, BOX_HALF_SIZE);

    // Immune cells
    int nCells = config["nbCells"]; // Number of cells at the beginning
    float ratio = config["ratioImmuneStromal"]; // Ratio Immune/Stromale
    int k = 0;
    Cell *c;
    while (k < nCells * ratio) {
      k++;
      MecaCell::Vec pos(nDist(MecaCell::Config::globalRand()),
                        nDist(MecaCell::Config::globalRand()),
                        nDist(MecaCell::Config::globalRand()));
      c = new Cell(pos, w.cellPlugin.diffusionPlugin.getGrid());
      c->init(Cell::Immune, Cell::Resident, &config,apopImmune,divImmune);
      w.addCell(c);
    }
    while (k < nCells) {
      k++;
      MecaCell::Vec pos(nDist(MecaCell::Config::globalRand()),
                        nDist(MecaCell::Config::globalRand()),
                        nDist(MecaCell::Config::globalRand()));
      c = new Cell(pos, w.cellPlugin.diffusionPlugin.getGrid());
      c->init(Cell::Stroma, Cell::None, &config,apopStroma,divStroma);
      w.addCell(c);
    }

    // Number of ImmuneMaker cells : normal distribution

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(config["mean"],config["variance"]);
    int nbImmuneMaker = d(gen);

    for (int k = 0; k < nbImmuneMaker; k++) {
      MecaCell::Vec pos(nDist(MecaCell::Config::globalRand()),
                        nDist(MecaCell::Config::globalRand()),
                        nDist(MecaCell::Config::globalRand()));
      c = new Cell(pos, w.cellPlugin.diffusionPlugin.getGrid());
      c->init(Cell::ImmuneMaker, Cell::None, &config,0,0);
      w.addCell(c);
    }

    w.addNewCells();
  }

  void loop() {
    

    if (!stop()) {

      w.update();

      if (w.getNbUpdates() == 1) {
        myfile << "nbIter"
              << ","
              << "nbCells"
              << ","
              << "nbAliveStroma"
              << ","
              << "nbApopStroma"
              << ","
              << "nbNecroStroma"
              << ","
              << "nbAliveImmune"
              << ","
              << "nbAliveImmuneResident"
              << ","
              << "nbAliveImmuneCirculatory"
              << ","
              << "nbApopImmune"
              << ","
              << "nbNecroImmune"
              << ","
              << "nbImmuneMaker"
              << ","
              << "avgIngla" << std::endl;
      }

      if (w.getNbUpdates() % 1 == 0) {

        

        // Création de fichiers (1 par step) qui contiennent les infos des cellules

        /* std::ofstream myfile;
        std::string fileName = "../CSV/csvStep"; // file name
        std::string extension = ".csv"; // extension
        std::stringstream fichier; // fichier à ouvrir
        fichier << fileName << w.getNbUpdates() << extension; // nouveau nom
        myfile.open(fichier.str().c_str(),std::ofstream::app); */

        int nbAliveStroma = 0;
        int nbApopStroma = 0;
        int nbAliveImmune = 0;
        int nbApopImmune = 0;
        int nbNecroStroma = 0;
        int nbAliveImmuneResident = 0;
        int nbAliveImmuneCirculatory = 0;
        int nbNecroImmune = 0;
        int nbImmuneMaker = 0;
        double avgInfla = 0;
        for (Cell *c : w.cells) {

          // Ecriture dans le fichier csv

          // myfile << c->type << "," << c->health << "," << c->infla << "," <<
          // c->reso << "," << c->getPosition() << std::endl;

          if (c->type == Cell::Immune) {
            if (c->state == Cell::Alive) {
              nbAliveImmune++;
              if (c->immuneType == Cell::Resident)
                nbAliveImmuneResident++;
              else if (c->immuneType == Cell::Circulatory)
                nbAliveImmuneCirculatory++;
            } else if (c->state == Cell::Apoptosis) {
              nbApopImmune++;
            } else if (c->state == Cell::Necrosis) {
              nbNecroImmune++;
            }
          } else if (c->type == Cell::Stroma) {
            if (c->state == Cell::Alive)
              nbAliveStroma++;
            else if (c->state == Cell::Apoptosis)
              nbApopStroma++;
            else if (c->state == Cell::Necrosis)
              nbNecroStroma++;
          } else if (c->type == Cell::ImmuneMaker) {
            nbImmuneMaker++;
          }
          avgInfla += c->getBody().getQuantities()[SIGNAL::INFLAMMATORY];
        }
        avgInfla /= w.cells.size();

        cerr << w.getNbUpdates() << "\t" << w.cells.size() << "\t"
             << nbAliveStroma << "\t" << nbApopStroma << "\t" << nbNecroStroma
             << "\t" << nbAliveImmune << "\t" << nbAliveImmuneResident << "\t"
             << nbAliveImmuneCirculatory << "\t" << nbApopImmune << "\t"
             << nbNecroImmune << "\t" << nbImmuneMaker << "\t" << avgInfla
             << endl;
             
        myfile << w.getNbUpdates() << "," << w.cells.size() << "," << nbAliveStroma << ","
           << nbApopStroma << "," << nbNecroStroma << "," << nbAliveImmune
           << "," << nbAliveImmuneResident << "," << nbAliveImmuneCirculatory
           << "," << nbApopImmune << "," << nbNecroImmune << ","
           << nbImmuneMaker << "," << avgInfla << std::endl;


        // myfile.close(); // Fermeture du csv contenant toutes les cellules
      }




      // if (w.getNbUpdates() == 10) {
      //   auto& w = getWorld();
      //   std::uniform_int_distribution<unsigned int> dist(0, w.cells.size());
      //   Cell *c = w.cells[dist(MecaCell::Config::globalRand())];
      //   c->state = Cell::Necrosis;
      //   for (auto* nc : c->getConnectedCells()) {
      //     for (auto *nnc : nc->getConnectedCells()) {
      //       nnc->state = Cell::Necrosis;
      //     }
      //   }
      // }
    }
  }

  bool stop() { return w.getNbUpdates() >= 2000; } // return false pour ne pas s'arreter

  // Import json file
  nlohmann::json loadJsonConfig(std::string fileName) {
    nlohmann::json config;
    std::ifstream i(fileName, std::ifstream::in);
    i >> config;
    i.close();
    return (config);
  }

  // Function to export values step by step (1row-> average of cells) to csv
  // void writeToCSV(std::ofstream myfile, int nbIter, int nbCells, int nbAliveStroma, int nbApopStroma,
  //                 int nbNecroStroma, int nbAliveImmune,
  //                 int nbAliveImmuneResident, int nbAliveImmuneCirculatory,
  //                 int nbApopImmune, int nbNecroImmune, int nbImmuneMaker,
  //                 float avgInfla) {
  //   //std::ofstream file = myfile;

  //   std::ofstream myfile;
    
  //   myfile << nbIter << "," << nbCells << "," << nbAliveStroma << ","
  //          << nbApopStroma << "," << nbNecroStroma << "," << nbAliveImmune
  //          << "," << nbAliveImmuneResident << "," << nbAliveImmuneCirculatory
  //          << "," << nbApopImmune << "," << nbNecroImmune << ","
  //          << nbImmuneMaker << "," << avgInfla << std::endl;
    
  // }
};

double GingiScenario::BOX_HALF_SIZE = 250.0;

#endif // _GINGISCENARIO_H
