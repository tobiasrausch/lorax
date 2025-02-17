#ifndef GAF_H
#define GAF_H


#include <iostream>
#include <fstream>

#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>

#include "util.h"

namespace lorax
{

  struct AlignRecord {
    int32_t qlen;
    int32_t qstart;
    int32_t qend;
    int32_t plen;
    int32_t pstart;
    int32_t pend;
    int32_t matches;
    int32_t alignlen;
    int32_t mapq;
    
    char strand;
    char hap;
    std::size_t seed;
    std::vector<std::pair<bool, uint32_t> > path; // simplified path (forward, tid)
    std::vector<int8_t> cigarop;
    std::vector<uint32_t> cigarlen;

    AlignRecord() : qstart(0), hap('*'), seed(0) {}
    AlignRecord(int32_t const q, std::size_t const s) : qstart(q), hap('*'), seed(s) {}

    bool operator<(const AlignRecord& s2) const {
      return ((seed < s2.seed) || ((seed == s2.seed) && (qstart < s2.qstart)));
    }
  };
  

  inline void
  parseGafCigar(std::string const& cigar, AlignRecord& ar) {
    uint32_t nstart = 0;
    uint32_t nend = 0;
    for(uint32_t i = 0; i < cigar.size(); ++i) {
      if (isdigit(cigar[i])) ++nend;
      else {
	uint32_t oplen = boost::lexical_cast<uint32_t>(cigar.substr(nstart, nend - nstart));
	ar.cigarlen.push_back(oplen);
	ar.cigarop.push_back(bam_cigar_table[(int) cigar[i]]);
	nstart = i + 1;
	nend = i + 1;
      }
    }
  }

  inline bool
  parseGafPath(std::string const& path, Graph const& g, AlignRecord& ar) {
    if (path.size()) {
      if ((path[0] == '>') || (path[0] == '<')) {
	std::size_t index = 0;
	std::vector<uint32_t> breaks;
	while ((index = path.find("<", index)) != std::string::npos) breaks.push_back(index++);
	index = 0;
	while ((index = path.find(">", index)) != std::string::npos) breaks.push_back(index++);
	std::sort(breaks.begin(), breaks.end());
	for(uint32_t i = 0; i < breaks.size(); ++i) {
	  bool forward = 1;
	  if (path[breaks[i]] == '<') forward = false;
	  std::string segment;
	  if (i + 1 < breaks.size()) segment = path.substr(breaks[i] + 1, breaks[i+1] - breaks[i] - 1);
	  else segment = path.substr(breaks[i] + 1);
	  typename Graph::TSegmentIdMap::const_iterator it = g.smap.find(segment);
	  if (it != g.smap.end()) ar.path.push_back(std::make_pair(forward, it->second));
	  else {
	    std::cerr << "Unknown segment " << segment << std::endl;
	    return false;
	  }
	}
      } else {
	std::cerr << "Unknown path format!" << std::endl;
	return false;
      }
    } else {
      std::cerr << "Empty path!" << std::endl;
      return false;
    }
    return true;
  }

  inline bool
  parseAlignRecord(std::istream& instream, Graph const& g, AlignRecord& ar, std::string& qname) {
    std::string gline;
    if(std::getline(instream, gline)) {      
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep("\t");
      Tokenizer tokens(gline, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter == tokens.end()) { std::cerr << "GAF parsing error!" << std::endl; return false; }	
      qname = *tokIter;
      ar.seed = hash_lr(qname);
      if (++tokIter == tokens.end()) { std::cerr << "GAF parsing error!" << std::endl; return false; }	
      ar.qlen = boost::lexical_cast<int32_t>(*tokIter);
      if (++tokIter == tokens.end()) { std::cerr << "GAF parsing error!" << std::endl; return false; }	
      ar.qstart = boost::lexical_cast<int32_t>(*tokIter);
      if (++tokIter == tokens.end()) { std::cerr << "GAF parsing error!" << std::endl; return false; }	
      ar.qend = boost::lexical_cast<int32_t>(*tokIter);
      if (++tokIter == tokens.end()) { std::cerr << "GAF parsing error!" << std::endl; return false; }
      ar.strand = boost::lexical_cast<char>(*tokIter);
      if (++tokIter == tokens.end()) { std::cerr << "GAF parsing error!" << std::endl; return false; }
      if (!g.empty()) {
	if (!parseGafPath(*tokIter, g, ar)) return false;
      }
      if (++tokIter == tokens.end()) { std::cerr << "GAF parsing error!" << std::endl; return false; }
      ar.plen = boost::lexical_cast<int32_t>(*tokIter);
      if (++tokIter == tokens.end()) { std::cerr << "GAF parsing error!" << std::endl; return false; }
      ar.pstart = boost::lexical_cast<int32_t>(*tokIter);
      if (++tokIter == tokens.end()) { std::cerr << "GAF parsing error!" << std::endl; return false; }
      ar.pend = boost::lexical_cast<int32_t>(*tokIter);
      if (++tokIter == tokens.end()) { std::cerr << "GAF parsing error!" << std::endl; return false; }
      ar.matches = boost::lexical_cast<int32_t>(*tokIter);
      if (++tokIter == tokens.end()) { std::cerr << "GAF parsing error!" << std::endl; return false; }
      ar.alignlen = boost::lexical_cast<int32_t>(*tokIter);
      if (++tokIter == tokens.end()) { std::cerr << "GAF parsing error!" << std::endl; return false; }
      ar.mapq = boost::lexical_cast<int32_t>(*tokIter);
      ++tokIter;
      for(; tokIter != tokens.end(); ++tokIter) {
	// Optional fields
	boost::char_separator<char> kvsep(":");
	Tokenizer tokopt(*tokIter, kvsep);
	Tokenizer::iterator tikv = tokopt.begin();
	if (*tikv == "cg") {
	  ++tikv; ++tikv;
	  parseGafCigar(*tikv, ar);
	}
      }
      return true;
    } else return false;
  }

  inline bool
  parseAlignRecord(std::istream& instream, Graph const& g, AlignRecord& ar) {
    std::string qname;
    return parseAlignRecord(instream, g, ar, qname);
  }
  
  inline bool
  parseAlignRecord(std::istream& instream, AlignRecord& ar) {
    Graph g;
    return parseAlignRecord(instream, g, ar);
  }
  
			    
  template<typename TConfig>
  inline bool
  parseGaf(TConfig const& c, Graph const& g, std::vector<AlignRecord>& aln, uint32_t const maxread, std::set<std::size_t>& seedSet) {
    // Open GAF
    std::ifstream gafFile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    if (is_gz(c.sample)) {
      gafFile.open(c.sample.string().c_str(), std::ios_base::in | std::ios_base::binary);
      dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
    } else gafFile.open(c.sample.string().c_str(), std::ios_base::in);
    dataIn.push(gafFile);

    // Parse GAF
    std::set<std::size_t> newReads;
    std::istream instream(&dataIn);
    bool parseAR = true;
    bool gafdone = true;
    while (parseAR) {
      AlignRecord ar;
      if (parseAlignRecord(instream, g, ar)) {
	if (maxread) {
	  if (seedSet.find(ar.seed) == seedSet.end()) {
	    // New read
	    if (newReads.find(ar.seed) == newReads.end()) {
	      if (newReads.size() >= maxread) {
		gafdone = false;
		continue; // No break here because all records are required for primary/secondary alignments
	      }
	      newReads.insert(ar.seed);
	    }
	  } else continue; // Already processed
	}
	aln.push_back(ar);
	//std::cerr << ar.seed << ',' << ar.qlen << ',' << ar.qstart << ',' << ar.qend << ',' << ar.strand << ',' << ar.plen << ',' << ar.pstart << ',' << ar.pend << ',' << ar.matches << ',' << ar.alignlen << ',' << ar.mapq << std::endl;
      }
      else parseAR = false;
    }
    // Keep track of processed reads
    if (maxread) seedSet.insert(newReads.begin(), newReads.end());

    // Close file
    dataIn.pop();
    if (is_gz(c.sample)) dataIn.pop();
    gafFile.close();

    // Sort alignment records by read
    std::sort(aln.begin(), aln.end());
    
    return gafdone;
  }

  template<typename TConfig>
  inline bool
  parseGaf(TConfig const& c, Graph const& g, std::vector<AlignRecord>& aln) {
    std::set<std::size_t> empty;
    return parseGaf(c, g, aln, 0, empty);
  }

}

#endif
