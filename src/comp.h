#ifndef COMP_H
#define COMP_H

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
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include "gfa.h"
#include "gaf.h"

namespace lorax
{
  
  struct CompConfig {
    bool splitGraph;
    std::string prefix;
    boost::filesystem::path outfile;
    boost::filesystem::path seqfile;
    boost::filesystem::path gfafile;
  };

  template<typename TConfig>
  inline void
  splitGraph(TConfig const& c, Graph const& g) {
    // Vertex map
    std::vector<std::string> idSegment(g.smap.size());
    for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;


    for(uint32_t k = 0; k < g.numComp; ++k) {
      // Temporary output file
      std::string filename = c.prefix + ".comp" + boost::lexical_cast<std::string>(k) + ".gfa";
      // Output rGFA
      std::ofstream sfile;
      sfile.open(filename.c_str());

      // Output graph
      faidx_t* fai = fai_load(c.seqfile.string().c_str());
      for(uint32_t i = 0; i < g.segments.size(); ++i) {
	if (g.segments[i].comp == k) {
	  std::string seqid = idSegment[i];
	  int32_t seqlen;
	  char* seq = faidx_fetch_seq(fai, seqid.c_str(), 0, faidx_seq_len(fai, seqid.c_str()), &seqlen);
	  sfile << "S\t" << seqid;
	  sfile << "\t" << seq;
	  sfile << "\tLN:i:" << g.segments[i].len;
	  sfile << "\tSN:Z:" << g.chrnames[g.segments[i].tid];
	  sfile << "\tSO:i:" << g.segments[i].pos;
	  sfile << "\tSR:i:" << g.ranks[g.segments[i].tid];
	  sfile << std::endl;
	  free(seq);
	}
      }
      fai_destroy(fai);

      // Output links
      for(uint32_t i = 0; i < g.links.size(); ++i) {
	if ((g.segments[g.links[i].from].comp == k) && (g.segments[g.links[i].to].comp == k)) {
	  sfile << "L\t" << idSegment[g.links[i].from];
	  if (g.links[i].fromfwd) sfile << "\t+";
	  else sfile << "\t-";
	  sfile << "\t" << idSegment[g.links[i].to];
	  if (g.links[i].tofwd) sfile << "\t+";
	  else sfile << "\t-";
	  sfile << "\t0M" << std::endl;
	}
      }
      sfile.close();
    }
  }
  

  template<typename TConfig>
  inline int32_t
  runComp(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Load pan-genome graph
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Load pan-genome graph" << std::endl;
    Graph g;
    if (c.splitGraph) parseGfa(c, g, true);
    else parseGfa(c, g, false);

    // Connected components
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Compute connected components" << std::endl;
    connectedComponents(g);

    // Vertex map
    std::vector<std::string> idSegment(g.smap.size());
    for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;
    
    // Open output file
    std::streambuf * buf;
    std::ofstream of;
    if(c.outfile.string() != "-") {
      of.open(c.outfile.string().c_str());
      buf = of.rdbuf();
    } else {
      buf = std::cout.rdbuf();
    }
    std::ostream out(buf);
    out << "segment\tlength\tcomponent" << std::endl;
    for(uint32_t i = 0; i < g.segments.size(); ++i) {
      out << idSegment[i] << "\t" << g.segments[i].len << "\t" << g.segments[i].comp << std::endl;
    }
    // Close file
    if(c.outfile.string() != "-") of.close();

    // Split graph
    if (c.splitGraph) {
      std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Split graph" << std::endl;
      splitGraph(c, g);

      // Clean-up
      boost::filesystem::remove(c.seqfile.string());
      boost::filesystem::remove(c.seqfile.string() + ".fai");
    }


    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  
  int components(int argc, char** argv) {
    CompConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("prefix,p", boost::program_options::value<std::string>(&c.prefix), "output prefix to split graph into components")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output file")
      ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("gfafile", boost::program_options::value<boost::filesystem::path>(&c.gfafile), "input file")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("gfafile", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("gfafile"))) {
      std::cerr << "Usage:" << std::endl;
      std::cerr << "lorax " << argv[0] << " [OPTIONS] <pangenome.hg38.gfa.gz>" << std::endl;
      std::cerr << visible_options << "\n";
      return -1;
    }

    // Split graph
    if (vm.count("prefix")) {
      c.splitGraph = true;
      c.seqfile = boost::filesystem::unique_path().replace_extension(".fa");
    }
    else c.splitGraph = false;
    
    // Check outfile
    if (!vm.count("outfile")) c.outfile = "-";
    else {
      if (c.outfile.string() != "-") {
	if (!_outfileValid(c.outfile)) return 1;
      }
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cerr << "lorax ";
    for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
    std::cerr << std::endl;
    
    return runComp(c);
  }

}

#endif
