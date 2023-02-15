#ifndef DOT_H
#define DOT_H

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
  
  struct DotConfig {
    uint32_t radius;
    std::string nodeid;
    boost::filesystem::path outfile;
    boost::filesystem::path seqfile;
    boost::filesystem::path gfafile;
  };


  template<typename TConfig>
  inline int32_t
  runGfa2Dot(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Load pan-genome graph
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Load pan-genome graph" << std::endl;
    Graph g;
    parseGfa(c, g, false);

    // Vertex map
    std::vector<std::string> idSegment(g.smap.size());
    for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;

    // Which nodes should be plotted?
    std::vector<bool> nodeplot;
    if (c.nodeid == "all") nodeplot.resize(g.smap.size(), true);
    else {
      nodeplot.resize(g.smap.size(), false);
      // Get neighborhood
      typedef std::set<uint32_t> TGang;
      TGang gang;
      gang.insert(g.smap[c.nodeid]);
      for(uint32_t r = 0; r<c.radius; ++r) {
	TGang newmembers;
	for(uint32_t i = 0; i < g.links.size(); ++i) {
	  if (gang.find(g.links[i].from) != gang.end()) newmembers.insert(g.links[i].to);
	  if (gang.find(g.links[i].to) != gang.end()) newmembers.insert(g.links[i].from);
	}
	gang.insert(newmembers.begin(), newmembers.end());
      }
      for(typename TGang::iterator it = gang.begin(); it != gang.end(); ++it) nodeplot[*it] = true;
    }

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
    out << "strict digraph {" << std::endl;
    for(uint32_t i = 0; i < g.segments.size(); ++i) {
      if (nodeplot[i]) out << "   " << idSegment[i] << " [label=\"" << idSegment[i] << "\"];" << std::endl;
    }
    for(uint32_t i = 0; i < g.links.size(); ++i) {
      if ((nodeplot[g.links[i].from]) && (nodeplot[g.links[i].to])) {
	out << "   " << idSegment[g.links[i].from] << " -> " << idSegment[g.links[i].to];
	out << " [dir=both,arrowtail=";
	if (g.links[i].fromfwd) out << "inv";
	else out << "normal";
	out << ",arrowhead=";
	if (g.links[i].tofwd) out << "normal";
	else out << "inv";
	out << "];" << std::endl;
      }
    }
    out << "}" << std::endl;
    
    // Close file
    if(c.outfile.string() != "-") of.close();

#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  
  int gfa2dot(int argc, char** argv) {
    DotConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("radius,r", boost::program_options::value<uint32_t>(&c.radius)->default_value(1), "radius around selected node")
      ("nodeid,n", boost::program_options::value<std::string>(&c.nodeid)->default_value("all"), "node to plot (all: all nodes)")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output dot file")
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
    
    return runGfa2Dot(c);
  }

}

#endif
