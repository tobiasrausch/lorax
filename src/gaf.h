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

  struct Path {
    bool forward;
    uint32_t tid;
    int32_t start;
    int32_t end;

    Path() {}
    Path(bool const f, uint32_t const t, uint32_t const s, uint32_t const e) : forward(f), tid(t), start(s), end(e) {}
  };
    

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
    std::size_t seed;
    std::vector<Path> path;  // stable ids (sequence coordinates)
    std::vector<int8_t> cigarop;
    std::vector<uint32_t> cigarlen;

    AlignRecord() : qstart(0), seed(0) {}
    AlignRecord(int32_t const q, std::size_t const s) : qstart(q), seed(s) {}
  };

  template<typename TRecord>
  struct SortAlignRecord : public std::binary_function<TRecord, TRecord, bool> {
    inline bool operator()(TRecord const& s1, TRecord const& s2) const {
      return ((s1.seed < s2.seed) || ((s1.seed == s2.seed) && (s1.qstart < s2.qstart)));
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
	  if (it != g.smap.end()) ar.path.push_back(Path(forward, it->second, 0, 0));
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


  template<typename TChrMap>
  inline bool
  parseGafPathStableId(std::string const& path, TChrMap const& chrmap, AlignRecord& ar) {
    if (path.size()) {
      if ((path[0] == '>') || (path[0] == '<')) {
	//std::cerr << path << std::endl;
	std::size_t index = 0;
	std::vector<uint32_t> breaks;
	while ((index = path.find("<", index)) != std::string::npos) breaks.push_back(index++);
	index = 0;
	while ((index = path.find(">", index)) != std::string::npos) breaks.push_back(index++);
	std::sort(breaks.begin(), breaks.end());
	for(uint32_t i = 0; i < breaks.size(); ++i) {
	  bool forward = true;
	  if (path[breaks[i]] == '<') forward = false;
	  std::string coords;
	  if (i + 1 < breaks.size()) coords = path.substr(breaks[i] + 1, breaks[i+1] - breaks[i] - 1);
	  else coords = path.substr(breaks[i] + 1);
	  index = coords.find(":", 0);
	  if (index != std::string::npos) {
	    std::string chrn = coords.substr(0, index);
	    coords = coords.substr(index+1);
	    index = coords.find("-", 0);
	    if (index != std::string::npos) {
	      int32_t pstart = boost::lexical_cast<int32_t>(coords.substr(0, index));
	      int32_t pend = boost::lexical_cast<int32_t>(coords.substr(index+1));
	      //std::cerr << path[breaks[i]] << "," << chrn << "," << pstart << "," << pend << "," << (int) forward <<  std::endl;
	      typename TChrMap::const_iterator it = chrmap.find(chrn);
	      if (it == chrmap.end()) {
		std::cerr << "Unknown chromosome name " << chrn << std::endl;
		return false;
	      }
	      ar.path.push_back(Path(forward, it->second, pstart, pend));
	    } else {
	      std::cerr << "Path lacks coordinates!" << std::endl;
	      return false;
	    }
	  } else {
	    std::cerr << "Path lacks chromosome name and coordinates!" << std::endl;
	    return false;
	  }
	}
      } else {
	// Only stable id
	typename TChrMap::const_iterator it = chrmap.find(path);
	if (it == chrmap.end()) {
	  std::cerr << "Unknown chromosome name " << path << std::endl;
	  return false;
	}
	ar.path.push_back(Path(true, it->second, 0, ar.plen));
      }
    }
    return true;
  }

  template<typename TConfig>
  inline bool
  parseGaf(TConfig const& c, Graph const& g, std::vector<AlignRecord>& aln) {
    // Chromosome map
    typedef std::map<std::string, uint32_t> TChrMap;
    TChrMap chrmap; // Each chromosome is mapped to a rank and globally unique integer id
    if (c.seqCoords) {
      for(uint32_t i = 0; i < g.chrnames.size(); ++i) chrmap.insert(std::make_pair(g.chrnames[i], i));
    }

    // Parse GAF
    uint32_t id = 0;
    std::ifstream gafFile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    if (is_gz(c.sample)) {
      gafFile.open(c.sample.string().c_str(), std::ios_base::in | std::ios_base::binary);
      dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
    } else gafFile.open(c.sample.string().c_str(), std::ios_base::in);
    dataIn.push(gafFile);
    std::istream instream(&dataIn);
    std::string gline;
    while(std::getline(instream, gline)) {      
      aln.resize(id + 1, AlignRecord());
		   
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep("\t");
      Tokenizer tokens(gline, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter != tokens.end()) {
	std::string qname = *tokIter;
	aln[id].seed = hash_lr(qname);
	++tokIter;
	if (tokIter != tokens.end()) {
	  aln[id].qlen = boost::lexical_cast<int32_t>(*tokIter);
	  ++tokIter;
	  if (tokIter != tokens.end()) {
	    aln[id].qstart = boost::lexical_cast<int32_t>(*tokIter);
	    ++tokIter;
	    if (tokIter != tokens.end()) {
	      aln[id].qend = boost::lexical_cast<int32_t>(*tokIter);
	      ++tokIter;
	      if (tokIter != tokens.end()) {
		aln[id].strand = boost::lexical_cast<char>(*tokIter);
		++tokIter;
		if (tokIter != tokens.end()) {
		  std::string path;
		  if (c.seqCoords) path = *tokIter;
		  else {
		    // Vertex coordinates
		    if (!parseGafPath(*tokIter, g, aln[id])) return false;
		  }
		  ++tokIter;
		  if (tokIter != tokens.end()) {
		    aln[id].plen = boost::lexical_cast<int32_t>(*tokIter);
		    ++tokIter;
		    if (tokIter != tokens.end()) {
		      aln[id].pstart = boost::lexical_cast<int32_t>(*tokIter);
		      ++tokIter;
		      if (tokIter != tokens.end()) {
			aln[id].pend = boost::lexical_cast<int32_t>(*tokIter);
			++tokIter;
			if (c.seqCoords) {
			  // Sequence coordinates
			  if (!parseGafPathStableId(path, chrmap, aln[id])) return false;
			}
			if (tokIter != tokens.end()) {
			  aln[id].matches = boost::lexical_cast<int32_t>(*tokIter);
			  ++tokIter;
			  if (tokIter != tokens.end()) {
			    aln[id].alignlen = boost::lexical_cast<int32_t>(*tokIter);
			    ++tokIter;
			    if (tokIter != tokens.end()) {
			      aln[id].mapq = boost::lexical_cast<int32_t>(*tokIter);
			      ++tokIter;
			      for(; tokIter != tokens.end(); ++tokIter) {
				// Optional fields
				boost::char_separator<char> kvsep(":");
				Tokenizer tokens(*tokIter, kvsep);
				Tokenizer::iterator tikv = tokens.begin();
				if (*tikv == "cg") {
				  ++tikv; ++tikv;
				  parseGafCigar(*tikv, aln[id]);
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      ++id; // Next alignment record
    }
    dataIn.pop();
    if (is_gz(c.sample)) dataIn.pop();
    gafFile.close();

    // Sort alignment records by read
    std::sort(aln.begin(), aln.end(), SortAlignRecord<AlignRecord>());
    
    return true;
  }
  
}

#endif
