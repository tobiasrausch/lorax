#ifndef NCOV_H
#define NCOV_H

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

  struct NCovConfig {
    bool basecov;
    std::string name;
    boost::filesystem::path outfile;
    boost::filesystem::path seqfile;
    boost::filesystem::path gfafile;
    boost::filesystem::path sample;
  };


  template<typename TConfig, typename TSegmentCoverage>
  inline void
  outputNodeCoverage(TConfig const& c, Graph const& g, TSegmentCoverage& segcov) {
    typedef typename TSegmentCoverage::value_type TBpCoverage;
    typedef typename TBpCoverage::value_type TValue;
    uint32_t maxCoverage = std::numeric_limits<TValue>::max();
    
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

    // Base coverage
    if (c.basecov) ofile << "sample\tsegment\tseglen\tpos\tcoverage" << std::endl;

    // Output coverage
    std::vector<TValue> medvec(g.segments.size(), 0);
    for(uint32_t k = 0; k < g.segments.size(); ++k) {
      if (c.basecov) {
	std::string segname = idSegment[k];
	uint32_t seglen = segcov[k].size();
	for(uint32_t p = 0; p < seglen; ++p) {
	  ofile << c.name << '\t' << segname << '\t' << seglen << '\t' << p+1 << '\t' << segcov[k][p] << std::endl;
	}
      } else {
	// Attention: Changes order of base coverage!!!
	auto m = segcov[k].begin() + segcov[k].size() / 2;
	std::nth_element(segcov[k].begin(), m, segcov[k].end());
	medvec[k] = segcov[k][segcov[k].size() / 2];
      }
    }

    // Summary statistics
    if (!c.basecov) {
      std::vector<TValue> filtered;
      // Drop segments with 0 or maximum coverage
      for(uint32_t i = 0; i < medvec.size(); ++i) {
	if ((medvec[i] > 0) && (medvec[i] < maxCoverage)) filtered.push_back(medvec[i]);
      }
      auto m = filtered.begin() + filtered.size() / 2;
      std::nth_element(filtered.begin(), m, filtered.end());
      TValue median = filtered[filtered.size() / 2];

      // Output normalized coverage
      ofile << "sample\tsegment\tseglen\tcoverage\tcnest" << std::endl;
      for(uint32_t k = 0; k < g.segments.size(); ++k) {
	std::string segname = idSegment[k];
	uint32_t seglen = segcov[k].size();
	double cnest = (double) medvec[k] / (double) median;
	ofile << c.name << '\t' << segname << '\t' << seglen << '\t' << medvec[k] << '\t' << cnest << std::endl;
      }
    }
  }
  
  template<typename TConfig, typename TSegmentCoverage>
  inline void
  parseGafCoverage(TConfig const& c, Graph const& g, TSegmentCoverage& segcov) {
    typedef typename TSegmentCoverage::value_type TBpCoverage;
    typedef typename TBpCoverage::value_type TMaxCoverage;
    uint32_t maxCoverage = std::numeric_limits<TMaxCoverage>::max();
    
    // Vertex map
    //std::vector<std::string> idSegment(g.smap.size());
    //for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;
    
    // Parse GAF
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
      //if (parseAlignRecord(instream, g, ar, qname)) {
      if (parseAlignRecord(instream, g, ar)) {
	uint32_t refstart = 0;
	for(uint32_t i = 0; i < ar.path.size(); ++i) {
	  //std::string seqname = idSegment[ar.path[i].second];
	  uint32_t seqlen = g.segments[ar.path[i].second].len;
	  uint32_t plen = seqlen;
	  if (i == 0) plen -= ar.pstart;
	  if (i + 1 == ar.path.size()) plen = ar.pend - ar.pstart - refstart;
	  uint32_t refend = refstart + plen;
	  //std::cerr << refstart << ',' << refend << ':' << ar.pstart << ',' << ar.pend << ';' << seqlen << ',' << refstart << std::endl;
	  uint32_t rp = 0;
	  uint32_t sp = 0;
	  std::string cigout;
	  for (uint32_t ci = 0; ci < ar.cigarop.size(); ++ci) {
	    if ((ar.cigarop[ci] == BAM_CEQUAL) || (ar.cigarop[ci] == BAM_CDIFF)) {
	      for(uint32_t k = 0; k < ar.cigarlen[ci]; ++k, ++sp, ++rp) {
		if ((rp >= refstart) && (rp < refend)) {
		  cigout += "M";
		}
	      }
	    } else if (ar.cigarop[ci] == BAM_CDEL) {
	      for(uint32_t k = 0; k < ar.cigarlen[ci]; ++k, ++rp) {
		if ((rp >= refstart) && (rp < refend)) cigout += "D";
	      }
	    } else if (ar.cigarop[ci] == BAM_CINS) {
	      for(uint32_t k = 0; k < ar.cigarlen[ci]; ++k, ++sp) {
		if ((rp >= refstart) && (rp < refend)) cigout += "I";
	      }
	    } else {
	      std::cerr << "Warning: Unknown Cigar option " << ar.cigarop[ci] << std::endl;
	    }
	  }
	  if (!ar.path[i].first) std::reverse(cigout.begin(), cigout.end());
	  // Compute coverage
	  rp = 0;
	  for(uint32_t k = 0; k < cigout.size(); ++k) {
	    if (cigout[k] == 'M') {
	      if (segcov[ar.path[i].second][rp] < maxCoverage) ++segcov[ar.path[i].second][rp];
	      ++rp;
	    } else if (cigout[k] == 'D') ++rp;
	  }	      
	  // Next segment
	  //std::cerr << qname << ',' << seqname << ',' << cigout << std::endl;
	  refstart = refend;
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
  runNCov(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Load pan-genome graph
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Load pan-genome graph" << std::endl;
    Graph g;
    parseGfa(c, g, false);

    // Allocate coverage vector
    typedef uint16_t TMaxCoverage;
    typedef std::vector<TMaxCoverage> TBpCoverage;
    typedef std::vector<TBpCoverage> TSegmentCoverage;
    TSegmentCoverage segcov(g.segments.size());
    for(uint32_t i = 0; i < g.segments.size(); ++i) segcov[i].resize(g.segments[i].len, 0);
    
    // Parse alignments
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse GAF alignments" << std::endl;
    parseGafCoverage(c, g, segcov);

    // Output
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Output node coverage" << std::endl;
    outputNodeCoverage(c, g, segcov);
    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  
  int ncov(int argc, char** argv) {
    NCovConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("graph,g", boost::program_options::value<boost::filesystem::path>(&c.gfafile), "GFA pan-genome graph")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output statistics")
      ("name,n", boost::program_options::value<std::string>(&c.name), "sample name")
      ("basecov,b", "basepair coverage")
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
      std::cout << "lorax " << argv[0] << " [OPTIONS] -g <pangenome.gfa.gz> <sample.gaf.gz>" << std::endl;
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

    // Base coverage
    if (vm.count("basecov")) c.basecov = true;
    else c.basecov = false;

    // Sample name
    if (c.name.empty()) {
      std::string delim(".");
      std::vector<std::string> parts;
      boost::split(parts, c.sample.string(), boost::is_any_of(delim));
      c.name = parts[0];
    }
    std::cout << c.name << std::endl;
    exit(-1);
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "lorax ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runNCov(c);
  }

}

#endif