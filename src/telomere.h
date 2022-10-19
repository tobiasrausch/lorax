#ifndef TELOMERE_H
#define TELOMERE_H

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


#include "junction.h"

namespace lorax
{

  struct TelomereConfig {
    typedef std::set<std::string> TStringSet;

    uint16_t minMapQual;
    uint32_t minClip;
    uint32_t minSupport;
    uint32_t period;
    uint32_t replen;
    uint32_t delta;
    uint32_t maxTelLen;
    uint32_t minChrLen;

    std::string outprefix;    
    std::string repeatStr;
    TStringSet fwdmotif;
    TStringSet revmotif;
    boost::filesystem::path genome;
    boost::filesystem::path sample;
  };


  struct TelomereRecord {
    typedef boost::dynamic_bitset<> TBitSet;
    std::string qname;
    std::string sequence;
    TBitSet telfwd;
    TBitSet telrev;
  };
  
  struct Mapping {
    int32_t tid;
    int32_t gstart;
    int32_t gend;
    int32_t rstart;
    int32_t rend;
    bool fwd;
    uint16_t qual;
    std::size_t seed;

    Mapping(int32_t const t, int32_t const gs, int32_t const ge, int32_t const rs, int32_t const re, bool const val, uint16_t const qval, std::size_t const sin) : tid(t), gstart(gs), gend(ge), rstart(rs), rend(re), fwd(val), qual(qval), seed(sin) {}
  };


  template<typename TConfig, typename TReadBp, typename TEdges>
  inline void
  links(TConfig const& c, TReadBp const& bp, TEdges& es) {
    typedef typename TReadBp::mapped_type TJunctionVector;
    
    // Segments are nodes, read connections define edges
    for(typename TReadBp::const_iterator itFirst = bp.begin(); itFirst != bp.end(); ++itFirst) {
      typename TReadBp::const_iterator itSecond = itFirst;
      ++itSecond;
      for(; itSecond != bp.end(); ++itSecond) {
	// Any common junction?
	bool run = true;
	for(typename TJunctionVector::const_iterator itI = itFirst->second.begin(); ((run) && (itI != itFirst->second.end())); ++itI) {
	  for(typename TJunctionVector::const_iterator itJ = itSecond->second.begin(); ((run) && (itJ != itSecond->second.end())); ++itJ) {
	    if ((itI->refidx == itJ->refidx) && (itI->scleft == itJ->scleft) && (std::abs(itI->refpos - itJ->refpos) < c.delta)) {
	      es.insert(std::make_pair(itFirst->first, itSecond->first));
	      run = false;
	    }
	  }
	}
      }
    }
  }
  
  template<typename TConfig, typename TReadBp>
  inline void
  concomp(TConfig const& c, TReadBp const& bp, std::map<std::size_t, uint32_t>& comp) {
    typedef std::map<std::size_t, uint32_t> TComponentMap;
    
    // Compute edges
    typedef std::pair<std::size_t, std::size_t> TEdge;
    std::set<TEdge> es;
    links(c, bp, es);

    // Put each node in its unique component
    comp.clear();
    uint32_t idx = 0;
    for(typename TReadBp::const_iterator it = bp.begin(); it != bp.end(); ++it, ++idx) comp.insert(std::make_pair(it->first, idx));
    
    // Merge components based on edges
    for(typename TReadBp::const_iterator itFirst = bp.begin(); itFirst != bp.end(); ++itFirst) {
      typename TReadBp::const_iterator itSecond = itFirst;
      ++itSecond;
      for(; itSecond != bp.end(); ++itSecond) {
	if (es.find(std::make_pair(itFirst->first, itSecond->first)) != es.end()) {
	  // Vertices in different components?
	  if (comp[itFirst->first] != comp[itSecond->first]) {
	    uint32_t oldcomp = comp[itSecond->first];
	    for(typename TComponentMap::iterator cit = comp.begin(); cit != comp.end(); ++cit) {
	      if (cit->second == oldcomp) cit->second = comp[itFirst->first];
	    }
	  }
	}
      }
    }
  }
  
  template<typename TConfig, typename TReadBp>
  inline void
  collectCandidates(TConfig const& c, TReadBp& readBp, std::map<std::size_t, TelomereRecord>& candidates) {
    uint32_t motiflen = c.period * c.replen;
    
    // Open file handles
    samFile* samfile = sam_open(c.sample.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.sample.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Parse BAM
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing... " << hdr->target_name[refIndex] << std::endl;

      // Iterate alignments 
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	std::size_t seed = hash_lr(rec);
	if (readBp.find(seed) != readBp.end()) {

	  // Load sequence
	  TelomereRecord telrec;
	  telrec.qname = std::string(bam_get_qname(rec));
	  telrec.sequence.resize(rec->core.l_qseq, 'N');
	  uint8_t* seqptr = bam_get_seq(rec);
	  for (int32_t i = 0; i < rec->core.l_qseq; ++i) telrec.sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	  if (rec->core.flag & BAM_FREVERSE) reverseComplement(telrec.sequence);  // Junctions are always reported on the forward strand
	  telrec.telfwd.resize(rec->core.l_qseq, 0);
	  telrec.telrev.resize(rec->core.l_qseq, 0);
	
	  // Mark telomere hits
	  for(std::set<std::string>::const_iterator it = c.fwdmotif.begin(); it != c.fwdmotif.end(); ++it) {
	    std::size_t index = 0;
	    while ((index = telrec.sequence.find(*it, index)) != std::string::npos) telrec.telfwd[index++] = 1;
	  }
	  for(std::set<std::string>::const_iterator it = c.revmotif.begin(); it != c.revmotif.end(); ++it) {
	    std::size_t index = 0;
	    while ((index = telrec.sequence.find(*it, index)) != std::string::npos) telrec.telrev[index++] = 1;
	  }

	  // Candidate telomere read
	  if (telrec.telfwd.count() + telrec.telrev.count() >= motiflen) {
	    // Debug
	    //std::cerr << seed << '\t' << bam_get_qname(rec) << '\t' << telrec.telfwd.count() << '\t' << telrec.telrev.count() << std::endl;
	    //for(uint32_t k = 0; k < readBp[seed].size(); ++k) std::cerr << hdr->target_name[readBp[seed][k].refidx] << '\t' << readBp[seed][k].refpos << '\t' << readBp[seed][k].seqpos << '\t' << readBp[seed][k].forward << '\t' << readBp[seed][k].qual << std::endl;
	    candidates.insert(std::make_pair(seed, telrec));
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }


  template<typename TConfig, typename TMappings, typename TTelomereReads>
  inline void
  outputTelomereComponents(TConfig const& c, std::map<std::size_t, uint32_t>& comp, TMappings& mp, TTelomereReads& telreads) {
    typedef std::map<std::size_t, uint32_t> TComponentMap;
    
    // Open file handles
    samFile* samfile = sam_open(c.sample.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.sample.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Pre-compute support
    std::vector<uint32_t> support(comp.size(), 0);
    for(typename TComponentMap::const_iterator cit = comp.begin(); cit != comp.end(); ++cit) ++support[cit->second];
    
    // Mappings
    std::string filename = c.outprefix + ".mappings.tsv";
    std::ofstream ofile(filename.c_str());    
    ofile << "chr\trefstart\trefend\treadname\treadstart\treadend\tforward\tcomponentid" << std::endl;
    for(uint32_t i = 0; i < mp.size(); ++i) {
      ofile << hdr->target_name[mp[i].tid] << '\t' << mp[i].gstart << '\t' << mp[i].gend << '\t' << telreads[mp[i].seed].qname << '\t' << mp[i].rstart << '\t' << mp[i].rend << '\t' << (int) (mp[i].fwd) << '\t' << comp[mp[i].seed] << std::endl;
      //'\t' << comp[mp[i].seed] << '\t' << support[mp[i].seed] << '\t' << "0" << '\t' << "0" << std::endl;
    }
    ofile.close();


    // Reads
    filename = c.outprefix + ".reads.tsv";
    std::ofstream rfile(filename.c_str());
    rfile << "component\tsupport\treadname\tfwdtelomere\trevtelomere" << std::endl;
    for(typename TTelomereReads::iterator it = telreads.begin(); it != telreads.end(); ++it) {
      rfile << comp[it->first] << '\t' << support[comp[it->first]] << '\t' << it->second.qname << '\t' << it->second.telfwd.count() << '\t' << it->second.telrev.count() << std::endl;
    }
    rfile.close();

    // Fasta files for each component
    for(uint32_t cidx = 0; cidx < support.size(); ++cidx) {
      if (support[cidx] > 0) {
	filename = c.outprefix + ".comp" + boost::lexical_cast<std::string>(cidx) + ".fa";
	std::ofstream ffile(filename.c_str());
	for(typename TTelomereReads::iterator it = telreads.begin(); it != telreads.end(); ++it) {
	  if (comp[it->first] == cidx) {
	    ffile << ">" << it->second.qname << std::endl;
	    ffile << it->second.sequence << std::endl;
	  }
	}
	ffile.close();
      }
    }

    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }

  template<typename TConfig, typename TTelomereReads, typename TMappings>
  inline void
  mappings(TConfig const& c, TTelomereReads const& telreads, TMappings& mp) {
    // Open file handles
    samFile* samfile = sam_open(c.sample.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.sample.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Parse BAM
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing... " << hdr->target_name[refIndex] << std::endl;

      // Iterate alignments 
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY)) continue;
	if (rec->core.qual < c.minMapQual) continue;
	std::size_t seed = hash_lr(rec);
	if (telreads.find(seed) != telreads.end()) {	
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
	  bool dir = true;
	  if (rec->core.flag & BAM_FREVERSE) {
	    dir = false;
	    int32_t seqTmp = seqStart;
	    seqStart = sp - seqEnd;
	    seqEnd = sp - seqTmp;
	  }
	  if (gpStart < gpEnd) mp.push_back(Mapping(rec->core.tid, gpStart, gpEnd, seqStart, seqEnd, dir, rec->core.qual, seed));
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }
    

  template<typename TConfig>
  inline int32_t
  runTelomere(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Breakpoints
    typedef std::vector<Junction> TJunctionVector;
    typedef std::map<std::size_t, TJunctionVector> TReadBp;
    TReadBp readBp;
    findJunctions(c, readBp);

    // Compute connected components
    {
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Connected components for " << readBp.size() << " candidates." << std::endl;
      typedef std::map<std::size_t, uint32_t> TComponentMap;
      TComponentMap comp;
      concomp(c, readBp, comp);
      std::vector<uint32_t> support(comp.size(), 0);
      for(typename TComponentMap::const_iterator cit = comp.begin(); cit != comp.end(); ++cit) ++support[cit->second];
      std::vector<std::size_t> delkeys;
      for(TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
	if (support[comp[it->first]] < c.minSupport) delkeys.push_back(it->first);
      }
      for(uint32_t i = 0; i < delkeys.size(); ++i) readBp.erase(delkeys[i]);
    }
    
    // Candidate tumor reads
    typedef std::map<std::size_t, TelomereRecord> TReadTelomere;
    TReadTelomere telreads;
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Telomere content screening for " << readBp.size() << " reads." << std::endl;
    collectCandidates(c, readBp, telreads);

    // Clean-up junctions
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Intersection of telomere and SV reads" << std::endl;
    {
      std::vector<std::size_t> delkeys;
      for(TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
	if (telreads.find(it->first) == telreads.end()) delkeys.push_back(it->first);
      }
      for(uint32_t i = 0; i < delkeys.size(); ++i) readBp.erase(delkeys[i]);
    }

    // Connected components
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Connected components for " << readBp.size() << " SV reads and " << telreads.size() << " telomere reads." << std::endl;
    typedef std::map<std::size_t, uint32_t> TComponentMap;
    TComponentMap comp;
    concomp(c, readBp, comp);
    std::vector<uint32_t> support(comp.size(), 0);
    for(typename TComponentMap::const_iterator cit = comp.begin(); cit != comp.end(); ++cit) ++support[cit->second];
    std::vector<std::size_t> delkeys;
    for(TReadBp::const_iterator it = readBp.begin(); it != readBp.end(); ++it) {
      if (support[comp[it->first]] < c.minSupport) delkeys.push_back(it->first);
    }
    for(uint32_t i = 0; i < delkeys.size(); ++i) {
      readBp.erase(delkeys[i]);
      telreads.erase(delkeys[i]);
    }

    // Get alignments
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Telomeric break identification using  " << readBp.size() << " SV reads and " << telreads.size() << " telomere reads." << std::endl;
    std::vector<Mapping> tumor_mp;
    mappings(c, telreads, tumor_mp);
    
    // Data out
    outputTelomereComponents(c, comp, tumor_mp, telreads);

#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  
  int telomere(int argc, char** argv) {
    TelomereConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Options");
    generic.add_options()
      ("help,?", "show help message")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(10), "min. mapping quality")
      ("minclip,c", boost::program_options::value<uint32_t>(&c.minClip)->default_value(25), "min. clipping length")
      ("support,s", boost::program_options::value<uint32_t>(&c.minSupport)->default_value(3), "min. telomere read support")
      ("delta,d", boost::program_options::value<uint32_t>(&c.delta)->default_value(250), "max. breakpoint junction offset")
      ("tellen,t", boost::program_options::value<uint32_t>(&c.maxTelLen)->default_value(20000), "max. telomere length")
      ("chrlen,l", boost::program_options::value<uint32_t>(&c.minChrLen)->default_value(40000000), "min. chromosome length")
      ("repeats,r", boost::program_options::value<std::string>(&c.repeatStr)->default_value("TTAGGG,TCAGGG,TGAGGG,TTGGGG"), "repeat units")
      ("period,p", boost::program_options::value<uint32_t>(&c.period)->default_value(3), "repeat period")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("out"), "output prefix")
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
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
      std::cout << "Usage: lorax " << argv[0] << " [OPTIONS] -g <ref.fa> <tumor.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Create telomere motifs
    createRepeatMotifs(c, true);

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "lorax ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runTelomere(c);
  }

}

#endif
