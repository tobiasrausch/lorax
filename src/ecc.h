#ifndef ECC_H
#define ECC_H

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

namespace lorax
{

  struct EccConfig {
    uint16_t minQual;
    uint32_t minOverlap;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
    boost::filesystem::path tumor;
    boost::filesystem::path control;
  };


  struct ReadMatch {
    int32_t tid;
    int32_t gstart;
    int32_t gend;
    int32_t rstart;
    int32_t rend;
    bool fwd;
    bool begalign;
    std::size_t seed;

    ReadMatch(int32_t const t, int32_t const gs, int32_t const ge, int32_t const rs, int32_t const re, bool f, bool b, std::size_t se) : tid(t), gstart(gs), gend(ge), rstart(rs), rend(re), fwd(f), begalign(b), seed(se) {}
  };


  template<typename TReadMatch>
  struct SortReadMatch : public std::binary_function<TReadMatch, TReadMatch, bool>
  {
    inline bool operator()(TReadMatch const& rm1, TReadMatch const& rm2) {
      return ((rm1.seed < rm2.seed) || ((rm1.seed == rm2.seed) && (rm1.begalign)));
    }
  };

  struct CircleSegment {
    int32_t cid;
    int32_t tid;
    int32_t gstart;
    int32_t gend;
    int32_t rstart;
    int32_t rend;
    bool fwd;
    uint16_t qual;
    std::size_t seed;
    std::string qname;

    CircleSegment(int32_t const cd, int32_t const t, int32_t const gs, int32_t const ge, int32_t const rs, int32_t const re, bool const val, uint16_t const qval, std::size_t s, std::string const& qn) : cid(cd), tid(t), gstart(gs), gend(ge), rstart(rs), rend(re), fwd(val), qual(qval), seed(s), qname(qn) {}
  };

  
  
  inline int32_t
  overlap(std::vector<ReadMatch> const& candidates, int32_t const i, int32_t const j) {
    if (candidates[i].tid != candidates[j].tid) return 0;
    if ((candidates[i].gend < candidates[j].gstart) || (candidates[j].gend < candidates[i].gstart)) return 0;
    std::vector<int32_t> coord(4);
    coord[0] = candidates[i].gstart;
    coord[1] = candidates[i].gend;
    coord[2] = candidates[j].gstart;
    coord[3] = candidates[j].gend;
    std::sort(coord.begin(), coord.end());

    // Contained segments cannot be
    if ((coord[1] == candidates[i].gstart) && (coord[2] == candidates[i].gend)) return 0;
    if ((coord[1] == candidates[j].gstart) && (coord[2] == candidates[j].gend)) return 0;

    // Overlap must be in expected order
    if (candidates[i].fwd) {
      // both forward
      if (candidates[i].begalign) {
	if (candidates[i].gstart < candidates[j].gstart) return 0;
      } else {
	if (candidates[j].gstart < candidates[i].gstart) return	0;
      }
    } else {
      // both reverse
      if (candidates[i].begalign) {
	if (candidates[i].gend > candidates[j].gend) return 0;
      } else {
	if (candidates[j].gend > candidates[i].gend) return 0;
      }
    }

    // Return overlap size
    return coord[2] - coord[1];
  }


  template<typename TConfig>
  inline void
  circleSegments(TConfig const& c, std::string const& filename, std::set<std::size_t> const& candidates, std::vector<CircleSegment>& mp, std::vector<CircleSegment> const& ctrl) {
    // For the tumor run
    std::set<std::size_t> invalidReads;
    bool tumor_run = false;
    if (filename == c.tumor.string()) tumor_run = true;

    // Open file handles
    samFile* samfile = sam_open(filename.c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, filename.c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Parse BAM
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing... " << hdr->target_name[refIndex] << std::endl;

      // Mask control regions
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet mask;
      if (tumor_run) {
	mask.resize(hdr->target_len[refIndex], 0);
	for(uint32_t i = 0; i < ctrl.size(); ++i) {
	  if (ctrl[i].tid == refIndex) {
	    for(int32_t k = ctrl[i].gstart; ((k < ctrl[i].gend) && (k < (int) hdr->target_len[refIndex])); ++k) mask[k] = 1;
	  }
	}
      }
      
      // Iterate alignments 
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY)) continue;
	std::size_t seed = hash_string(bam_get_qname(rec));
	if (candidates.find(seed) != candidates.end()) {
	  
	  // Parse CIGAR
	  uint32_t* cigar = bam_get_cigar(rec);
	  int32_t gp = rec->core.pos; // Genomic position
	  int32_t gpStart = -1; //Match start
	  int32_t gpEnd = -1; //Match end
	  int32_t sp = 0; // Sequence position
	  int32_t seqStart = -1;  // Match start
	  int32_t seqEnd = -1; // Match end
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	      if (seqStart == -1) {
		seqStart = sp;
		gpStart = gp;
	      }
	      gp += bam_cigar_oplen(cigar[i]);
	      sp += bam_cigar_oplen(cigar[i]);
	      seqEnd = sp;
	      gpEnd = gp;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	      if (seqStart == -1) {
		seqStart = sp;
		gpStart = gp;
	      }
	      sp += bam_cigar_oplen(cigar[i]);
	      seqEnd = sp;
	      gpEnd = gp;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	      if (seqStart == -1) {
		seqStart = sp;
		gpStart = gp;
	      }
	      gp += bam_cigar_oplen(cigar[i]);
	      seqEnd = sp;
	      gpEnd = gp;
	    } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	      gp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else {
	      std::cerr << "Warning: Unknown Cigar options!" << std::endl;
	    }
	  }

	  // Alignment direction
	  bool fwd = true;
	  if (rec->core.flag & BAM_FREVERSE) {
	    fwd = false;
	    int32_t seqTmp = seqStart;
	    seqStart = sp - seqEnd;
	    seqEnd = sp - seqTmp;
	  }
	  if (gpStart < gpEnd) {
	    bool insertSegment = true;
	    if (tumor_run) {
	      // For the tumor, check that mappings are absent in the control
	      for(int32_t k = gpStart; ((k < gpEnd) && (k < (int) hdr->target_len[refIndex])); ++k) {
		if (mask[k]) {
		  insertSegment = false;
		  invalidReads.insert(seed);
		  break;
		}
	      }
	    }
	    //std::cerr << bam_get_qname(rec) << ',' << hdr->target_name[rec->core.tid] << ',' << gpStart << std::endl;
	    if (insertSegment) mp.push_back(CircleSegment(mp.size(), rec->core.tid, gpStart, gpEnd, seqStart, seqEnd, fwd, rec->core.qual, seed, bam_get_qname(rec)));
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }

    // Filter reads masked by control
    if (invalidReads.size()) {
    }

    // Debug
    std::cerr << invalidReads.size() << std::endl;
    if (tumor_run) {
      for(uint32_t i = 0; i < mp.size(); ++i) {
	if (invalidReads.find(mp[i].seed) == invalidReads.end()) {
	  std::cerr << hdr->target_name[mp[i].tid] << '\t' << mp[i].gstart << '\t' << mp[i].gend << '\t' << mp[i].qname << '\t' << mp[i].rstart << '\t' << mp[i].rend << '\t' << (int) mp[i].fwd << std::endl;
	}
      }
    }
    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }

  
  template<typename TConfig>
  inline int32_t
  circleCandidates(TConfig const& c, std::string const& filename, std::set<std::size_t>& selreads) {
    std::vector<ReadMatch> candidates;
  
    // Open file handles
    samFile* samfile = sam_open(filename.c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, filename.c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Parse BAM
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing... " << hdr->target_name[refIndex] << std::endl;

      // Iterate alignments 
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FUNMAP)) continue;
	if (rec->core.qual <= c.minQual) continue;

	// Parse CIGAR
	uint32_t* cigar = bam_get_cigar(rec);
	int32_t gp = rec->core.pos; // Genomic position
	int32_t gpStart = -1; //Match start
	int32_t gpEnd = -1; //Match end
	int32_t sp = 0; // Sequence position
	int32_t seqStart = -1;  // Match start
	int32_t seqEnd = -1; // Match end
	uint32_t softClip = 0;
	for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	  if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	    if (seqStart == -1) {
	      seqStart = sp;
	      gpStart = gp;
	    }
	    gp += bam_cigar_oplen(cigar[i]);
	    sp += bam_cigar_oplen(cigar[i]);
	    seqEnd = sp;
	    gpEnd = gp;
	  } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	    if (seqStart == -1) {
	      seqStart = sp;
	      gpStart = gp;
	    }
	    sp += bam_cigar_oplen(cigar[i]);
	    seqEnd = sp;
	    gpEnd = gp;
	  } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	    if (seqStart == -1) {
	      seqStart = sp;
	      gpStart = gp;
	    }
	    gp += bam_cigar_oplen(cigar[i]);
	    seqEnd = sp;
	    gpEnd = gp;
	  } else if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
	    if (bam_cigar_oplen(cigar[i]) >= c.minOverlap) {
	      if (sp == 0) softClip += 1;
	      else softClip += 2;
	    }
	    sp += bam_cigar_oplen(cigar[i]);
	  } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	    gp += bam_cigar_oplen(cigar[i]);
	  } else std::cerr << "Unknown Cigar options" << std::endl;
	}
	bool fwd = true;
	if (rec->core.flag & BAM_FREVERSE) {
	  fwd = false;
	  int32_t seqTmp = seqStart;
	  seqStart = sp - seqEnd;
	  seqEnd = sp - seqTmp;
	}
	// Debug
	//std::cerr << hdr->target_name[refIndex] << '\t' << gpStart << '\t' << gpEnd << '\t' << bam_get_qname(rec) << '\t' << seqStart << '\t' << seqEnd << '\t' << (int) fwd << '\t' << (int) rec->core.qual << std::endl;
	
	// Leading / trailing soft-clip of sufficient length ?
	if ((softClip > 0) && (softClip < 3)) {
	  if ((seqStart < (int32_t) c.minOverlap) && (seqEnd + (int32_t) c.minOverlap < sp)) {
	    // Read begin alignment
	    std::size_t seed = hash_string(bam_get_qname(rec));
	    candidates.push_back(ReadMatch(refIndex, gpStart, gpEnd, seqStart, seqEnd, fwd, true, seed));
	  } else if ((seqEnd + (int32_t) c.minOverlap > sp) && (seqStart > (int32_t) c.minOverlap)) {
	    // Read end alignment
	    std::size_t seed = hash_string(bam_get_qname(rec));
	    candidates.push_back(ReadMatch(refIndex, gpStart, gpEnd, seqStart, seqEnd, fwd, false, seed));
	  } else {
	    // Nop
	    // ToDo
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }

    // Sort matches
    std::sort(candidates.begin(), candidates.end(), SortReadMatch<ReadMatch>());

    // Select circular reads
    for(uint32_t i = 0; i < candidates.size(); ++i) {
      for(uint32_t j = i + 1; j < candidates.size(); ++j) {
	if (candidates[i].seed != candidates[j].seed) break;
	if (candidates[i].fwd == candidates[j].fwd) {
	  if (candidates[i].begalign != candidates[j].begalign) {
	    if (overlap(candidates, i, j) > (int32_t) c.minOverlap) {
	      selreads.insert(candidates[i].seed);
	      break;
	    }
	  }
	}
      }
    }
    
    // Debug
    //for(uint32_t i = 0; i < candidates.size(); ++i) {
    //if (selreads.find(candidates[i].seed) != selreads.end()) {
    //std::cerr << hdr->target_name[candidates[i].tid] << '\t' << candidates[i].gstart << '\t' << candidates[i].gend << '\t' << candidates[i].seed << '\t' << candidates[i].rstart << '\t' << candidates[i].rend << '\t' << (int) candidates[i].fwd << '\t' << (int) candidates[i].begalign << std::endl;
    //}
    //}
    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    return selreads.size();
  }

  template<typename TConfig>
  inline int32_t
  runEccDna(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Tumor reads
    std::set<std::size_t> treads;
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Collecting tumor BAM candidates." << std::endl;
    int32_t tsize = circleCandidates(c, c.tumor.string(), treads);

    // Control reads
    std::set<std::size_t> creads;
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Collecting control BAM candidates." << std::endl;
    int32_t csize = circleCandidates(c, c.control.string(), creads);
    std::cerr << tsize << ',' << csize << std::endl;

    // Control regions
    std::vector<CircleSegment> ctrlsegments;
    std::vector<CircleSegment> tumsegments;
    circleSegments(c, c.control.string(), creads, ctrlsegments, tumsegments);

    // Tumor somatic eccDNA
    circleSegments(c, c.tumor.string(), treads, tumsegments, ctrlsegments);

    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    return 0;
  }


  
  int eccDna(int argc, char** argv) {
    EccConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Options");
    generic.add_options()
      ("help,?", "show help message")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minQual)->default_value(10), "min. mapping quality")
      ("overlap,p", boost::program_options::value<uint32_t>(&c.minOverlap)->default_value(200), "min. overlap")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("matched,m", boost::program_options::value<boost::filesystem::path>(&c.control), "matched control BAM")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.bed.gz"), "gzipped output file")
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

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "lorax ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runEccDna(c);
  }

}

#endif
