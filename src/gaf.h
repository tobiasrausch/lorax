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
    std::size_t seed;
    std::vector<int8_t> cigarop;
    std::vector<uint32_t> cigarlen;
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

  
  inline void
  parseGafPath(std::string const& path, AlignRecord& ar) {
    std::cerr << path << std::endl;
  }
  
  template<typename TConfig>
  inline void
  parseGaf(TConfig& c, std::vector<AlignRecord>& aln) {
    // Parse GAF
    uint32_t id = 0;
    std::ifstream gafFile(c.sample.string().c_str(), std::ios_base::in);
    if(gafFile.is_open()) {
      while (gafFile.good()) {
	std::string gline;
	getline(gafFile, gline);
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
		    parseGafPath(*tokIter, aln[id]);
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
      gafFile.close();
    }
  }


  

}

#endif
