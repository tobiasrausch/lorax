#ifndef ECOV_H
#define ECOV_H

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

  struct EcovConfig {
    std::string name;
    boost::filesystem::path outfile;
    boost::filesystem::path gfafile;
    boost::filesystem::path seqfile;
    boost::filesystem::path sample;
  };

  template<typename TConfig>
  inline void
  outputEdgeCoverage(TConfig const& c, Graph const& g, std::vector< std::pair<uint32_t, uint32_t> >& sm) {
    // Vertex map
    std::vector<std::string> idSegment(g.smap.size());
    for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;
    
    // Output file
    std::streambuf* buf;
    std::ofstream of;
    if (c.outfile != "-") {
      of.open(c.outfile.string().c_str());
      buf = of.rdbuf();
    } else {
      buf = std::cout.rdbuf();
    }
    std::ostream ofile(buf);

    // Summary statistics
    ofile << "from\tfromOrient\tto\ttoOrient\tsupport\tmapq" << std::endl;

    // Output
    for(uint32_t i = 0; i < g.links.size(); ++i) {
      if (sm[i].first > 0) sm[i].second /= sm[i].first;
      ofile << idSegment[g.links[i].from];
      if (g.links[i].fromfwd) ofile << "\t+";
      else ofile << "\t-";
      ofile << "\t" << idSegment[g.links[i].to];
      if (g.links[i].tofwd) ofile << "\t+";
      else ofile << "\t-";
      ofile << "\t" << sm[i].first;
      ofile << "\t" << sm[i].second;
      ofile << std::endl;
    }

    /*
    for(uint32_t i = 0; i < g.links.size(); ++i) {
      if (sm[i].first > 0) sm[i].second /= sm[i].first;
      ofile << g.chrnames[g.segments[g.links[i].from].tid];
      if (g.links[i].fromfwd) {
	ofile << "\t" << (g.segments[g.links[i].from].pos + g.segments[g.links[i].from].len);
	ofile << "\t+";
      }
      else {
	ofile << "\t" << g.segments[g.links[i].from].pos;
	ofile << "\t-";
      }
      ofile << "\t" << g.chrnames[g.segments[g.links[i].to].tid];
      if (g.links[i].tofwd) {
	ofile << "\t" << g.segments[g.links[i].to].pos;
	ofile << "\t+";
      } else {
	ofile << "\t" << (g.segments[g.links[i].to].pos + g.segments[g.links[i].to].len);
	ofile << "\t-";
      }
      ofile << "\t" << sm[i].first;
      ofile << "\t" << sm[i].second;
      if (medsup > 0) ofile << "\t" << (double) sm[i].first / (double) medsup << std::endl;
      else ofile << "\tNA" << std::endl;
    }
    */
  }

  template<typename TConfig>
  inline void
  coverageLinks(TConfig const& c, Graph const& g, std::vector< std::pair<uint32_t, uint32_t> >& sm) {
    typedef std::vector<Link> TLinks;
    
    // Open GAF
    std::ifstream gafFile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    if (is_gz(c.sample)) {
      gafFile.open(c.sample.string().c_str(), std::ios_base::in | std::ios_base::binary);
      dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
    } else gafFile.open(c.sample.string().c_str(), std::ios_base::in);
    dataIn.push(gafFile);

    // Parse GAF
    std::istream instream(&dataIn);
    bool parseAR = true;
    while (parseAR) {
      AlignRecord ar;
      std::string qname;
      if (parseAlignRecord(instream, g, ar, qname)) {
	for(uint32_t i = 0; i < ar.path.size(); ++i) {
	  uint32_t j = i + 1;
	  if (j < ar.path.size()) {
	    Link lk(ar.path[i].first, ar.path[j].first, ar.path[i].second, ar.path[j].second);
	    typename TLinks::const_iterator iter = std::lower_bound(g.links.begin(), g.links.end(), lk, SortLinks<Link>());
	    bool found = false;
	    for(; ((iter != g.links.end()) && (iter->from == lk.from)); ++iter) {
	      if (iter->to == lk.to) {
		if ((iter->fromfwd == lk.fromfwd) && (iter->tofwd == lk.tofwd)) {
		  found = true;
		  break;
		}
	      }
	    }
	    if (!found) {
	      Link lkSwap(!ar.path[j].first, !ar.path[i].first, ar.path[j].second, ar.path[i].second);
	      iter = std::lower_bound(g.links.begin(), g.links.end(), lkSwap, SortLinks<Link>());
	      for(;((iter != g.links.end()) && (iter->from == lkSwap.from)); ++iter) {
		if (iter->to == lkSwap.to) {
		  if ((iter->fromfwd == lkSwap.fromfwd) && (iter->tofwd == lkSwap.tofwd)) {
		    found = true;
		    break;
		  }
		}
	      }
	    }
	    if (!found) {
	      std::cerr << "Error: Inconsistent alignment edge!" << std::endl;
	    }
	    
	    // Increase support
	    int index = iter - g.links.begin();
	    ++sm[index].first;
	    sm[index].second += ar.mapq;
	  }
	}
      } else parseAR = false;
    }

    // Close file
    dataIn.pop();
    if (is_gz(c.sample)) dataIn.pop();
    gafFile.close();
  }


  template<typename TConfig>
  inline int32_t
  runECov(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Load pan-genome graph
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Load pan-genome graph" << std::endl;
    Graph g;
    parseGfa(c, g, false);

    // Sort links
    std::sort(g.links.begin(), g.links.end(), SortLinks<Link>());
    
    // Coverage links
    std::vector< std::pair<uint32_t, uint32_t> > sm(g.links.size(), std::make_pair(0, 0));
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse GAF alignments" << std::endl;
    coverageLinks(c, g, sm);

    // Output
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Output edge coverage" << std::endl;
    outputEdgeCoverage(c, g, sm);
    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  
  int ecov(int argc, char** argv) {
    EcovConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("graph,g", boost::program_options::value<boost::filesystem::path>(&c.gfafile), "GFA pan-genome graph")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output statistics")
      ("name,n", boost::program_options::value<std::string>(&c.name), "sample name")
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

    // Check outfile
    if (!vm.count("outfile")) c.outfile = "-";
    else {
      if (c.outfile.string() != "-") {
	if (!_outfileValid(c.outfile)) return 1;
      }
    }

    // Sample name
    if (c.name.empty()) {
      std::string delim(".");
      std::vector<std::string> parts;
      boost::split(parts, c.sample.stem().string(), boost::is_any_of(delim));
      c.name = parts[0];
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "lorax ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runECov(c);
  }

}

#endif
