#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>

#define BOOST_DISABLE_ASSERTS

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

#include "util.h"
#include "version.h"
#include "tithread.h"
#include "telomere.h"
#include "amplicon.h"
#include "ecc.h"
#include "extract.h"
#include "repeat.h"
#include "pct.h"
#include "ecov.h"
#include "convert.h"
#include "dot.h"
#include "ncov.h"
#include "comp.h"
#include "stats.h"

using namespace lorax;

inline void
displayUsage() {
  std::cout << "Usage: lorax <command> <arguments>" << std::endl;
  std::cout << std::endl;
  std::cout << "Linear reference alignments:" << std::endl;
  std::cout << std::endl;
  std::cout << "    tithreads     templated insertion threads" << std::endl;
  std::cout << "    amplicon      amplicon read selection for targeted assembly" << std::endl;
  std::cout << "    pct           percent identity" << std::endl;
  std::cout << "    extract       extract matches and fasta for selected reads" << std::endl;
  std::cout << std::endl;
  std::cout << "Telomeres:" << std::endl;
  std::cout << std::endl;
  std::cout << "    telomere      reads with neo- and intra-telomere repeats" << std::endl;  
  std::cout << "    repeat        telomeric repeat counting" << std::endl;
  std::cout << std::endl;
  std::cout << "Pan-genome graphs:" << std::endl;
  std::cout << std::endl;
  std::cout << "    stats         basic graph statistics" << std::endl; 
  std::cout << "    components    connnected components of a pan-genome graph" << std::endl;
  std::cout << "    convert       convert pan-genome graph alignment to BAM" << std::endl;
  std::cout << "    gfa2dot       convert pan-genome graph to dot (graphviz) format" << std::endl;
  std::cout << "    ncov          node coverage" << std::endl;
  std::cout << "    ecov          edge coverage" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}

int main(int argc, char **argv) {
  if (argc < 2) { 
    printTitle("Lorax");
    displayUsage();
    return 0;
  }
  
  if ((std::string(argv[1]) == "version") || (std::string(argv[1]) == "--version") || (std::string(argv[1]) == "--version-only") || (std::string(argv[1]) == "-v")) {
    std::cout << "Lorax version: v" << loraxVersionNumber << std::endl;
    std::cout << " using Boost: v" << BOOST_VERSION / 100000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100 << std::endl;
    std::cout << " using HTSlib: v" << hts_version() << std::endl;
    return 0;
  }
  else if ((std::string(argv[1]) == "help") || (std::string(argv[1]) == "--help") || (std::string(argv[1]) == "-h") || (std::string(argv[1]) == "-?")) {
    printTitle("Lorax");
    displayUsage();
    return 0;
  }
  else if ((std::string(argv[1]) == "warranty") || (std::string(argv[1]) == "--warranty") || (std::string(argv[1]) == "-w")) {
    displayWarranty();
    return 0;
  }
  else if ((std::string(argv[1]) == "license") || (std::string(argv[1]) == "--license") || (std::string(argv[1]) == "-l")) {
    bsd();
    return 0;
  }
  else if ((std::string(argv[1]) == "tithreads")) {
    return tithread(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "telomere")) {
    return telomere(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "repeat")) {
    return repeat(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "amplicon")) {
    return amplicon(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "eccdna")) {
    return eccDna(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "extract")) {
    return extract(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "pct")) {
    return pct(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "ncov")) {
    return ncov(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "ecov")) {
    return ecov(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "convert")) {
    return convert(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "gfa2dot")) {
    return gfa2dot(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "components")) {
    return components(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "stats")) {
    return stats(argc-1,argv+1);
  } else {
    std::cerr << "Unrecognized command " << std::string(argv[1]) << std::endl;
    return 1;
  }
}

