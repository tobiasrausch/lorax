#ifndef GENO_H
#define GENO_H

#define BOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX

#include <fstream>
#include <iomanip>

#include <boost/dynamic_bitset.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/icl/split_interval_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include "gfa.h"
#include "gaf.h"

namespace lorax
{

  struct GenoConfig {
    bool seqCoords;
    std::string outprefix;
    std::string sampleid;
    boost::filesystem::path gfafile;
    boost::filesystem::path seqfile;
    boost::filesystem::path sample;
  };


  template<typename TConfig>
  inline bool
  genotypeLinks(TConfig const& c, Graph const& g, std::vector<AlignRecord> const& aln) {
    // Sort links
    typedef std::vector<LinkCargo> TLinks;
    TLinks links(g.links.size());
    for(uint32_t i = 0; i < g.links.size(); ++i) links[i] = g.links[i];
    std::sort(links.begin(), links.end(), SortLinks<LinkCargo>());

    // Iterate alignments
    for(uint32_t id = 0; id < aln.size(); ++id) {
      for(uint32_t i = 0; i < aln[id].path.size(); ++i) {
	uint32_t j = i + 1;
	if (j < aln[id].path.size()) {
	  Link lk(aln[id].path[i].first, aln[id].path[j].first, aln[id].path[i].second, aln[id].path[j].second);
	  typename TLinks::iterator iter = std::lower_bound(links.begin(), links.end(), lk, SortLinks<LinkCargo>());
	  bool found = false;
	  for(;((iter != links.end()) && (iter->from == lk.from) && (iter->to == lk.to));++iter) {
	    if ((iter->fromfwd == lk.fromfwd) && (iter->tofwd == lk.tofwd)) {
	      found = true;
	      break;
	    }
	  }
	  if (!found) {
	    Link lkSwap(!aln[id].path[j].first, !aln[id].path[i].first, aln[id].path[j].second, aln[id].path[i].second);
	    iter = std::lower_bound(links.begin(), links.end(), lkSwap, SortLinks<LinkCargo>());
	    for(;((iter != links.end()) && (iter->from == lkSwap.from) && (iter->to == lkSwap.to)); ++iter) {
	      if ((iter->fromfwd == lkSwap.fromfwd) && (iter->tofwd == lkSwap.tofwd)) {
		found = true;
		break;
	      }
	    }
	  }
	  if (!found) {
	    std::cerr << "Inconsistent alignment edge!" << std::endl;
	    return false;
	  }

	  // Increase support
	  ++iter->support;
	  iter->mapq += aln[id].mapq;
	}
      }
    }

    // Median support
    uint32_t medsup = 0;
    {
      std::vector<uint32_t> lsup;
      for(uint32_t i = 0; i < links.size(); ++i) {
	if (links[i].support > 0) lsup.push_back(links[i].support);
      }
      std::sort(lsup.begin(), lsup.end());
      medsup = lsup[lsup.size()/2];
    }
      
    // Output
    std::ofstream sfile;
    std::string filen = c.outprefix + ".tsv";
    sfile.open(filen.c_str());
    for(uint32_t i = 0; i < links.size(); ++i) {
      if (links[i].support > 0) links[i].mapq /= links[i].support;
      sfile << g.chrnames[g.segments[links[i].from].tid];
      if (links[i].fromfwd) {
	sfile << "\t" << (g.segments[links[i].from].pos + g.segments[links[i].from].len);
	sfile << "\t+";
      }
      else {
	sfile << "\t" << g.segments[links[i].from].pos;
	sfile << "\t-";
      }
      sfile << "\t" << g.chrnames[g.segments[links[i].to].tid];
      if (links[i].tofwd) {
	sfile << "\t" << g.segments[links[i].to].pos;
	sfile << "\t+";
      } else {
	sfile << "\t" << (g.segments[links[i].to].pos + g.segments[links[i].to].len);
	sfile << "\t-";
      }
      sfile << "\t" << links[i].support;
      sfile << "\t" << links[i].mapq;
      if (medsup > 0) sfile << "\t" << (double) links[i].support / (double) medsup << std::endl;
      else sfile << "\tNA" << std::endl;
    }
    sfile.close();
    return true;
  }


  template<typename TConfig>
  inline int32_t
  runGeno(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Load pan-genome graph
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Load pan-genome graph" << std::endl;
    Graph g;
    parseGfa(c, g);
    
    // Parse alignments
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse alignments" << std::endl;
    std::vector<AlignRecord> aln;
    parseGaf(c, g, aln);

    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Genotype links" << std::endl;
    if (!genotypeLinks(c, g, aln)) return -1;
    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  
  int geno(int argc, char** argv) {
    GenoConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("graph,g", boost::program_options::value<boost::filesystem::path>(&c.gfafile), "GFA pan-genome graph")
      ("sequences,s", boost::program_options::value<boost::filesystem::path>(&c.seqfile), "stable sequences")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("out"), "output tprefix")
      ("seqcoords,c", "GAF uses sequence coordinates")
      ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.sample), "input file")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("graph"))) {
      std::cout << "Usage:" << std::endl;
      std::cout << "lorax " << argv[0] << " [OPTIONS] -g <pangenome.hg38.gfa.gz> <sample.gaf>" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Sequence coordinates
    if (vm.count("seqcoords")) c.seqCoords = true;
    else c.seqCoords = false;

    // Sequence coordinates
    if (c.seqCoords) {
      // Seqfile has the stable sequences
      if (!vm.count("sequences")) {
	std::cerr << "Using stable sequences requires -s <stable.sequences.fa.gz>" << std::endl;
	return -1;
      }
    } else {
      if (vm.count("sequences")) {
	std::cerr << "Please do not specify stable sequences if -c is not used!" << std::endl;
	return -1;
      }
      // Temporary file for the node sequence information
      boost::uuids::uuid uuid = boost::uuids::random_generator()();
      c.seqfile = c.outprefix + "." + boost::lexical_cast<std::string>(uuid) + ".fa";
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "lorax ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runGeno(c);
  }

}

#endif
