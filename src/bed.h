#ifndef BED_H
#define BED_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/algorithm/string.hpp>

#include <htslib/sam.h>

#include "util.h"

namespace lorax
{

  struct BedEntry {
    uint32_t start;
    uint32_t end;
    uint32_t cid;

    BedEntry(uint32_t const s, uint32_t const e, uint32_t const c) : start(s), end(e), cid(c) {}
  };
    
  
  // Flattens overlapping intervals
  template<typename TConfig, typename TRegionsGenome>
  inline int32_t
  _parseBedIntervals(TConfig const& c, TRegionsGenome& bedRegions) {
    // Open file handles
    samFile* samfile = sam_open(c.tumor.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Typedef
    typedef typename TRegionsGenome::value_type TChrIntervals;
    
    int32_t intervals = 0;
    bedRegions.resize(hdr->n_targets, TChrIntervals());
    std::ifstream chrFile(c.bedfile.string().c_str(), std::ifstream::in);
    if (chrFile.is_open()) {
      while (chrFile.good()) {
	std::string chrFromFile;
	getline(chrFile, chrFromFile);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t,;");
	Tokenizer tokens(chrFromFile, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  std::string chrName = *tokIter++;
	  int32_t tid = bam_name2id(hdr, chrName.c_str());
	  if (tid >= 0) {
	    if (tokIter!=tokens.end()) {
	      int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	      int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	      bedRegions[tid].push_back(BedEntry(start, end, intervals));
	      ++intervals;
	    }
	  }
	}
      }
      chrFile.close();
    }

    // Clean-up
    bam_hdr_destroy(hdr);
    sam_close(samfile);

    return intervals;
  }


}

#endif
