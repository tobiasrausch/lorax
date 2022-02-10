#ifndef TITHREAD_H
#define TITHREAD_H

#include <fstream>
#include <iomanip>

#include <boost/dynamic_bitset.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
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

#include "extract.h"

namespace lorax
{

  struct CallConfig {
    bool extract;
    uint16_t minMapQual;
    uint16_t minClip;
    uint16_t minSplit;
    uint32_t minSegmentSize;
    uint32_t maxSegmentSize;
    uint32_t minChrLen;
    float contam;
    std::string outprefix;
    boost::filesystem::path genome;
    boost::filesystem::path tumor;
    boost::filesystem::path control;
  };

  struct Breakpoint {
    bool left;
    uint32_t pos;
    uint32_t splits;
    float obsexp;

    Breakpoint(bool const l, uint32_t const p, uint32_t const s, float oe) : left(l), pos(p), splits(s), obsexp(oe) {}
  };
 

  template<typename TBreakpoint>
  struct SortBreakpoints : public std::binary_function<TBreakpoint, TBreakpoint, bool>
  {
    inline bool operator()(TBreakpoint const& bp1, TBreakpoint const& bp2) {
      return ((bp1.pos < bp2.pos) || ((bp1.pos == bp2.pos) && (bp1.left)));
    }
  };


  struct Segment {
    uint32_t refIndex;
    uint32_t start;
    uint32_t end;
    uint32_t cid;
    float cn;

    Segment(uint32_t const c, uint32_t const s, uint32_t const e, uint32_t const uid, float const cnval) : refIndex(c), start(s), end(e), cid(uid), cn(cnval) {}
  };


  inline bool
  mappedReads(hts_idx_t* idx, int32_t refIndex, std::string const& str) {
    // Any data?
    std::string suffix("cram");
    if ((str.size() >= suffix.size()) && (str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0)) return true;
    uint64_t mapped = 0;
    uint64_t unmapped = 0;
    hts_idx_get_stat(idx, refIndex, &mapped, &unmapped);
    if (mapped) return true;
    else return false;
  }
  
  template<typename TConfig, typename TVector, typename TChrReadPos>
  inline void
  parseChr(TConfig& c, samFile* samfile, hts_idx_t* idx, bam_hdr_t* hdr, int32_t refIndex, TVector& left, TVector& right, TVector& cov, TChrReadPos& read, bool const trackreads) {
    typedef typename TVector::value_type TValue;
    // Max value
    TValue maxval = std::numeric_limits<TValue>::max();

    // Read alignments
    hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
    bam1_t* rec = bam_init1();
    while (sam_itr_next(samfile, iter, rec) >= 0) {
      if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FUNMAP)) continue;
      if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
      std::size_t seed = hash_string(bam_get_qname(rec));
	
      // SV detection using single-end read
      uint32_t rp = rec->core.pos; // reference pointer
      uint32_t sp = 0; // sequence pointer
      
      // Parse the CIGAR
      uint32_t* cigar = bam_get_cigar(rec);
      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	  for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
	    if (cov[rp] < maxval) ++cov[rp];
	    ++rp;
	    ++sp;
	  }
	} else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	  rp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	  sp += bam_cigar_oplen(cigar[i]);
	} else if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
	  if (bam_cigar_oplen(cigar[i]) >= c.minClip) {
	    if (sp == 0) {
	      if (left[rp] < maxval) left[rp] += 1;
	    } else {
	      if (right[rp] < maxval) right[rp] += 1;
	    }
	    if (trackreads) read.push_back(std::make_pair(seed, rp));
	  }
	  sp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	  rp += bam_cigar_oplen(cigar[i]);
	} else {
	  std::cerr << "Warning: Unknown Cigar operation!" << std::endl;
	}
      }
    }
    bam_destroy1(rec);
    hts_itr_destroy(iter);
  }


  template<typename TBitSet, typename TVector>
  inline void
  covParams(TBitSet const& nrun, TVector const& cov, uint32_t const seedwin, uint32_t& avgcov, uint32_t& sdcov) {
    // Collect coverage values
    std::vector<uint32_t> vcov;
    for(uint32_t i = seedwin; i < nrun.size(); i = i + seedwin) {
      uint32_t lcov = 0;
      bool ncontent = false;
      for(uint32_t k = i - seedwin; k < i; ++k) {
	if (nrun[k]) {
	  ncontent = true;
	  break;
	}
	lcov += cov[k];
      }
      if (!ncontent) vcov.push_back(lcov);
    }
    // Drop lowest and highest 25%
    uint32_t ist = 0;
    uint32_t ien = vcov.size();
    if (ien > 1000) {
      ist = (uint32_t) (0.25 * vcov.size());
      ien = (uint32_t) (0.75 * vcov.size());
      std::sort(vcov.begin(), vcov.end());
    }
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance> > acc;
    for(uint32_t i = ist; i < ien; ++i) acc(vcov[i]);
    sdcov = sqrt(boost::accumulators::variance(acc));
    avgcov = boost::accumulators::mean(acc);
  }

  template<typename TBitSet, typename TVector, typename TValue>
  inline bool
  getcov(TBitSet const& nrun, TVector const& cov, uint32_t const start, uint32_t const end, TValue& lcov) {
    lcov = 0;
    for(uint32_t k = start; k < end; ++k) {
      if (nrun[k]) return false;
      lcov += cov[k];
    }
    return true;
  }
  

  template<typename TReadSegment, typename TEdgeSupport>
  inline void
  computelinks(TReadSegment const& read, TEdgeSupport& es) {
    std::size_t oldid = 0;
    int32_t kleft = 0;
    for(uint32_t k = 0; k < read.size(); ++k) {
      if (read[k].first == oldid) {
	for(int32_t walk = k - 1; walk >= kleft; --walk) {
	  uint32_t id1 = read[k].second;
	  uint32_t id2 = read[walk].second;
	  if (id1 > id2) {
	    uint32_t tmp = id1;
	    id1 = id2;
	    id2 = tmp;
	  }
	  if (es.find(std::make_pair(id1, id2)) == es.end()) es.insert(std::make_pair(std::make_pair(id1, id2), 1));
	  else ++es[std::make_pair(id1, id2)];
	}
      } else {
	// New split-read
	oldid = read[k].first;
	kleft = k;
      }
    }
  }
  
  template<typename TConfig, typename TEdgeSupport, typename TSegments>
  inline void
  segconnect(TConfig const& c, TEdgeSupport& es, TSegments& sgm) {
    // Segments are nodes, split-reads are edges
    for(uint32_t id1 = 0; id1 < sgm.size(); ++id1) {
      for(uint32_t id2 = id1 + 1; id2 < sgm.size(); ++id2) {
	if (es.find(std::make_pair(id1, id2)) != es.end()) {
	  if (es[std::make_pair(id1, id2)] >= c.minSplit) {
	    // Vertices in different components?
	    if (sgm[id1].cid != sgm[id2].cid) {
	      uint32_t oldid = sgm[id2].cid;
	      for(uint32_t i = 0; i < sgm.size(); ++i) {
		if (sgm[i].cid == oldid) sgm[i].cid = sgm[id1].cid;
	      }
	    }
	  }
	}
      }
    }
  }
  

  
  template<typename TConfig>
  inline int32_t
  runTithread(TConfig& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Open file handles
    samFile* samfile = sam_open(c.tumor.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.tumor.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    samFile* cfile = sam_open(c.control.string().c_str(), "r");
    hts_set_fai_filename(cfile, c.genome.string().c_str());
    hts_idx_t* cidx = sam_index_load(cfile, c.control.string().c_str());
    faidx_t* fai = fai_load(c.genome.string().c_str());
    char* seq = NULL;
    
    // Parse genome, process chromosome by chromosome
    typedef std::vector<Segment> TSegments;
    TSegments sgm;
    typedef std::pair<std::size_t, uint32_t> TReadPos;
    typedef std::vector<TReadPos> TChrReadPos;
    TChrReadPos readSeg;
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      // Any data?
      if ((!mappedReads(idx, refIndex, c.tumor.string())) || (!mappedReads(idx, refIndex, c.control.string()))) continue;

      // Large enough chromosome?
      if (hdr->target_len[refIndex] > c.minChrLen) {
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();	  
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Parsing " << hdr->target_name[refIndex] << std::endl;
      } else continue;

      // Tumor
      std::vector<uint16_t> left(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> right(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> cov(hdr->target_len[refIndex], 0);
      TChrReadPos reads;
      parseChr(c, samfile, idx, hdr, refIndex, left, right, cov, reads, true);
      
      // Control
      std::vector<uint16_t> cleft(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> cright(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> ccov(hdr->target_len[refIndex], 0);
      parseChr(c, cfile, cidx, hdr, refIndex, cleft, cright, ccov, reads, false);

      // Load sequence
      int32_t seqlen = -1;
      seq = faidx_fetch_seq(fai, hdr->target_name[refIndex], 0, hdr->target_len[refIndex], &seqlen);
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet nrun(hdr->target_len[refIndex], 0);
      for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	if ((seq[i] == 'n') || (seq[i] == 'N')) nrun[i] = 1;
      }
      if (seq != NULL) free(seq);

      uint32_t seedwin = 2 * c.minSegmentSize;
      std::map<uint32_t, uint32_t> possegmentmap;
      if (2 * seedwin < hdr->target_len[refIndex]) {
	// Get background coverage
	uint32_t sdcov = 0;
	uint32_t avgcov = 0;
	covParams(nrun, cov, seedwin, avgcov, sdcov);
	//std::cerr << "Tumor avg. coverage and SD coverage " << avgcov << "," << sdcov << std::endl;
	uint32_t csdcov = 0;
	uint32_t cavgcov = 0;
	covParams(nrun, ccov, seedwin, cavgcov, csdcov);
	//std::cerr << "Control avg. coverage and SD coverage " << cavgcov << "," << csdcov << std::endl;
	float expratio = (float) (avgcov) / (float) (cavgcov);

	// Identify candidate breakpoints
	typedef std::vector<Breakpoint> TBreakpointVector;
	TBreakpointVector bpvec;
	for(uint32_t i = seedwin; i < hdr->target_len[refIndex] - seedwin; ++i) {
	  // Left soft-clips
	  if (left[i] >= c.minSplit) {
	    uint32_t threshold = (uint32_t) (c.contam * left[i]);
	    if (cleft[i] <= threshold) {
	      uint32_t lcov = 0;
	      uint32_t rcov = 0;
	      if (!getcov(nrun, cov, i - seedwin, i, lcov)) continue;
	      if (!getcov(nrun, cov, i, i+seedwin, rcov)) continue;
	      if ((lcov * 1.5 < rcov) && (rcov > avgcov + 3 * sdcov)) {
		uint32_t controllcov = 0;
		uint32_t controlrcov = 0;
		if (!getcov(nrun, ccov, i - seedwin, i, controllcov)) continue;
		if (!getcov(nrun, ccov, i, i+seedwin, controlrcov)) continue;
		if ((controllcov * 1.5 < controlrcov) || (controlrcov > cavgcov + 3 * csdcov)) continue;
		if (controlrcov > 0) {
		  float obsratio = rcov / controlrcov;
		  if (obsratio / expratio > 1.5) bpvec.push_back(Breakpoint(true, i, left[i], obsratio / expratio));
		}
	      }
	    }
	  }
	  // Right soft-clips
	  if (right[i] >= c.minSplit) {
	    uint32_t threshold = (uint32_t) (c.contam * right[i]);
	    if (cright[i] <= threshold) {
	      uint32_t lcov = 0;
	      uint32_t rcov = 0;
	      if (!getcov(nrun, cov, i - seedwin, i, lcov)) continue;
	      if (!getcov(nrun, cov, i, i+seedwin, rcov)) continue;
	      if ((rcov * 1.5 < lcov) && (lcov > avgcov + 3 * sdcov)) {
		uint32_t controllcov = 0;
		uint32_t controlrcov = 0;
		if (!getcov(nrun, ccov, i - seedwin, i, controllcov)) continue;
		if (!getcov(nrun, ccov, i, i+seedwin, controlrcov)) continue;
		if ((controlrcov * 1.5 < controllcov) || (controllcov > cavgcov + 3 * csdcov)) continue;
		if (controllcov > 0) {
		  float obsratio = lcov / controllcov;
		  if (obsratio / expratio > 1.5) bpvec.push_back(Breakpoint(false, i, right[i], obsratio / expratio));
		}
	      }
	    }
	  }
	}
	
	// Merge left and right breakpoints into candidate regions
	if (bpvec.size()) {
	  std::sort(bpvec.begin(), bpvec.end(), SortBreakpoints<Breakpoint>());
	  uint32_t lastRight = 0;
	  for(uint32_t i = 0; i < bpvec.size() - 1; ++i) {
	    if (i < lastRight) continue;
	    if ((bpvec[i].left) && (!bpvec[i+1].left) && (bpvec[i+1].pos - bpvec[i].pos < c.maxSegmentSize)) {
	      // Split-read switchpoint (extend if possible)
	      uint32_t bestLeft = i;
	      for(int32_t k = i - 1; k >= 0; --k) {
		if (!bpvec[k].left) break;
		if (bpvec[i].pos - bpvec[k].pos > c.maxSegmentSize) break;
		if (bpvec[k].obsexp / bpvec[i].obsexp < 0.5) break;
		bestLeft = k;
	      }
	      uint32_t bestRight = i + 1;
	      for(uint32_t k = i + 2; k < bpvec.size(); ++k) {
		if (bpvec[k].left) break;
		if (bpvec[k].pos - bpvec[i+1].pos > c.maxSegmentSize) break;
		if (bpvec[k].obsexp / bpvec[i+1].obsexp < 0.5) break;
		bestRight = k;
	      }
	      uint32_t segsize = bpvec[bestRight].pos - bpvec[bestLeft].pos;
	      if ((segsize > c.minSegmentSize) && (segsize < c.maxSegmentSize)) {
		// New candidate segment
		lastRight = bestRight;
		uint64_t tmrcov = 0;
		if (getcov(nrun, cov, bpvec[bestLeft].pos, bpvec[bestRight].pos, tmrcov)) {
		  uint64_t ctrcov = 0;
		  if (getcov(nrun, ccov, bpvec[bestLeft].pos, bpvec[bestRight].pos, ctrcov)) {
		    if (ctrcov > 0) {
		      float obsratio = (float) (tmrcov) / (float) (ctrcov);
		      float obsexp = obsratio / expratio;
		      if (obsexp > 1.5) {
			uint32_t lid = sgm.size();
			sgm.push_back(Segment(refIndex, bpvec[bestLeft].pos, bpvec[bestRight].pos, lid, obsexp * 2.0));  // Assumes diploid
			for(uint32_t k = bpvec[bestLeft].pos; k <= bpvec[bestRight].pos; ++k) {
			  // ToDo: Replace with interval tree !!!
			  possegmentmap.insert(std::make_pair(k, lid)); 
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
      // Carry-over all split-reads
      if (!possegmentmap.empty()) {
	for(uint32_t i = 0; i < reads.size(); ++i) {
	  if (possegmentmap.find(reads[i].second) != possegmentmap.end()) {
	    // Keep track of seed and segment
	    readSeg.push_back(std::make_pair(reads[i].first, possegmentmap[reads[i].second]));
	  }
	}
      }
    }

    // Compute links
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();	  
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Computing segment links" << std::endl;
    // Sort split-reads by read ID
    std::sort(readSeg.begin(), readSeg.end());
    // Output split-reads
    //for(uint32_t i = 0; i < readSeg.size(); ++i) std::cerr << readSeg[i].first << '\t' << readSeg[i].second << std::endl;
    
    // Edges
    typedef std::pair<uint32_t, uint32_t> TEdge;
    typedef std::map<TEdge, uint32_t> TEdgeSupport;
    TEdgeSupport es;
    computelinks(readSeg, es);
    
    // Segment connections    
    now = boost::posix_time::second_clock::local_time();	  
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Computing connected components" << std::endl;
    segconnect(c, es, sgm);

    // Filter singletons or clusters where all segments are nearby
    std::vector<bool> confirmed(sgm.size(), false);  // Confirmed by component id (cid)
    for(uint32_t i = 0; i < sgm.size(); ++i) {
      if (confirmed[sgm[i].cid]) continue;
      for(uint32_t j = i + 1; j < sgm.size(); ++j) {
	if (sgm[i].cid != sgm[j].cid) continue;
	// Different chromosomes?
	if (sgm[i].refIndex != sgm[j].refIndex) {
	  confirmed[sgm[i].cid] = true;
	  break;
	}
	// Different pos? (segments are ordered)
	if (sgm[j].start - sgm[i].end > c.maxSegmentSize) {
	  confirmed[sgm[i].cid] = true;
	  break;
	}
      }
    }

    // Compute node degree (without self edges)
    std::vector<uint32_t> degree(sgm.size(), 0);
    for(uint32_t i = 0; i < sgm.size(); ++i) {
      for(uint32_t j = i + 1; j < sgm.size(); ++j) {
	if (es.find(std::make_pair(i, j)) != es.end()) {
	  if (es[std::make_pair(i, j)] >= c.minSplit) {
	    degree[i] += es[std::make_pair(i, j)];
	    degree[j] += es[std::make_pair(i, j)];
	  }
	}
      }
    }
	
    // Output segments
    std::string filename = c.outprefix + ".bed";
    std::ofstream ofile(filename.c_str());
    ofile << "chr\tstart\tend\tnodeid\tselfdegree\tdegree\testcn\tclusterid\tedges" << std::endl;
    for(uint32_t i = 0; i < sgm.size(); ++i) {
      if (confirmed[sgm[i].cid]) {
	ofile << hdr->target_name[sgm[i].refIndex] << '\t' << sgm[i].start << '\t' << sgm[i].end << '\t';
	ofile << i << "[label=\"" << hdr->target_name[sgm[i].refIndex] << ':' << sgm[i].start << '-' << sgm[i].end << "(" << sgm[i].cid << ")" <<  "\"];" << '\t';
	if ((es.find(std::make_pair(i, i)) != es.end()) && (es[std::make_pair(i, i)] >= c.minSplit)) {
	  ofile << es[std::make_pair(i, i)] << '\t';
	} else {
	  ofile << "0\t";
	}
	ofile << degree[i] << '\t' << sgm[i].cn << '\t' << sgm[i].cid << '\t';
	for(uint32_t id2 = i; id2 < sgm.size(); ++id2) {
	  if (es.find(std::make_pair(i, id2)) != es.end()) {
	    if (es[std::make_pair(i, id2)] >= c.minSplit) {
	      ofile << i << "--" << id2 << "[label=\"" << es[std::make_pair(i, id2)] << "\"];";
	    }
	  }
	}
	ofile << std::endl;
      }
    }
    ofile.close();

    // Clean-up
    bam_hdr_destroy(hdr);
    fai_destroy(fai);
    hts_idx_destroy(idx);
    sam_close(samfile);
    hts_idx_destroy(cidx);
    sam_close(cfile);

    // Get valid reads
    std::set<std::size_t> selectedReads;
    std::vector<bool> segconfirmed(sgm.size(), false);  // Confirmed by original id (index)
    for(uint32_t i = 0; i < sgm.size(); ++i) {
      if (confirmed[sgm[i].cid]) segconfirmed[i] = true;
    }
    for(uint32_t i = 0; i < readSeg.size(); ++i) {
      if (segconfirmed[readSeg[i].second]) selectedReads.insert(readSeg[i].first);
    }

    // Extract read-level information
    if (c.extract) {
      now = boost::posix_time::second_clock::local_time();	  
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Fetch alignments for " << selectedReads.size() << " reads" << std::endl;
      ExtractConfig ec;
      ec.genome = c.genome;
      ec.bamfile = c.tumor;
      ec.outfile = boost::filesystem::path(c.outprefix + ".match.gz");
      ec.fafile = boost::filesystem::path(c.outprefix + ".fa.gz");
      fetchAlignments(ec, selectedReads);
    }
    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    return 0;
  }


  
  int tithread(int argc, char** argv) {
    CallConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Options");
    generic.add_options()
      ("help,?", "show help message")
      ("qual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("clip,c", boost::program_options::value<uint16_t>(&c.minClip)->default_value(25), "min. clipping length")
      ("split,s", boost::program_options::value<uint16_t>(&c.minSplit)->default_value(3), "min. split-read support")
      ("chrlen,l", boost::program_options::value<uint32_t>(&c.minChrLen)->default_value(40000000), "min. chromosome length")
      ("minsize,i", boost::program_options::value<uint32_t>(&c.minSegmentSize)->default_value(100), "min. segment size")
      ("maxsize,j", boost::program_options::value<uint32_t>(&c.maxSegmentSize)->default_value(10000), "max. segment size")
      ("contam,n", boost::program_options::value<float>(&c.contam)->default_value(0), "max. fractional tumor-in-normal contamination")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("matched,m", boost::program_options::value<boost::filesystem::path>(&c.control), "matched control BAM")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("out"), "output prefix")
      ("extract,e", "extract matches and sequences")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.tumor), "input file")
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
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome")) || (!vm.count("matched"))) {
      std::cout << "Usage: lorax " << argv[0] << " [OPTIONS] -g <ref.fa> -m <control.bam> <tumor.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Extract sequences
    if (vm.count("extract")) c.extract = true;
    else c.extract = false;

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "lorax ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runTithread(c);
  }

}

#endif
