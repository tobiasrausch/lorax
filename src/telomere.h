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
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>


#include "junction.h"

namespace lorax
{

  struct TelomereConfig {
    typedef boost::dynamic_bitset<> TBitSet;
    
    uint32_t minClip;
    uint32_t replen;
    uint32_t minChrLen;
    uint32_t win;
    uint32_t medwin;
    double threshold;
    std::string repeatStr;
    std::string outprefix;
    
    TBitSet fwdm;
    TBitSet revm;
    boost::filesystem::path genome;
    boost::filesystem::path sample;
  };


  struct TelomereRecord {
    bool fwd;
    uint32_t seqsize;
    uint32_t telStart;
    uint32_t telEnd;
    uint32_t telRegions;
    double avgTelSig;
    std::size_t seed;
    std::string qname;

    TelomereRecord() {}
    
    TelomereRecord(bool const f, uint32_t slen, uint32_t const tS, uint32_t const tE, uint32_t const tR, double const avgTS, std::size_t const sd, std::string const& qn) : fwd(f), seqsize(slen), telStart(tS), telEnd(tE), telRegions(tR), avgTelSig(avgTS), seed(sd), qname(qn) {}
  };

  double medianVector(std::vector<double> const& myVector) {
    std::vector<double> myVectorCopy = myVector;
    const auto middleItr = myVectorCopy.begin() + myVectorCopy.size() / 2;
    std::nth_element(myVectorCopy.begin(), middleItr, myVectorCopy.end());
    if (myVectorCopy.size() % 2 == 0) {
      const auto leftMiddleItr = std::max_element(myVectorCopy.begin(), middleItr);
      return (*leftMiddleItr + *middleItr) / 2.0;
    } else {
      return *middleItr;
    }
  }


  template<typename TConfig, typename TBitSet>
  inline bool
  repeatFinder(TConfig const& c, std::string const& qname, std::size_t const& seed, TBitSet const& tel, bool const fwd, std::vector<TelomereRecord>& telreads) {
    bool foundReps = false;
    uint32_t seqsize = tel.size();

    std::vector<double> avg(seqsize);
    uint32_t cumsum = 0;
    uint32_t lagging_cumsum = 0;
    for(uint32_t i = 0; i<c.win; ++i) cumsum += tel[i];
    uint32_t half_win = (c.win - 1) / 2;
    for(uint32_t k = 0; k <= half_win; ++k) avg[k] = cumsum/(double)(c.win);
    for(uint32_t i = c.win; i < seqsize; ++i) {
      cumsum += tel[i];
      lagging_cumsum += tel[i - c.win];
      avg[i - half_win] = (cumsum - lagging_cumsum)/(double)(c.win);
    }
    for(uint32_t k = (seqsize - half_win); k < seqsize; ++k) avg[k] = (cumsum - lagging_cumsum)/(double)(c.win);
    std::vector<double> medavg(seqsize);
    std::vector<double> rolling_median(c.medwin);
    half_win = (c.medwin - 1) / 2;
    for(uint32_t i = 0; i<c.medwin; ++i) rolling_median[i] = avg[i];
    double medInit = medianVector(rolling_median);
    for(uint32_t k = 0; k <= half_win; ++k) medavg[k] = medInit;
    for(uint32_t i = c.medwin; i < seqsize; ++i) {
      rolling_median[i % c.medwin] = avg[i];
      medavg[i - half_win] = medianVector(rolling_median);
    }
    medInit = medianVector(rolling_median);
    for(uint32_t k = (seqsize - half_win); k < seqsize; ++k) medavg[k] = medInit;
    // Telomere runs
    std::vector<uint32_t> beginTel;
    std::vector<uint32_t> endTel;
    if (medavg[0] > c.threshold) beginTel.push_back(0);
    for(uint32_t i = 0; i < (seqsize - 1); ++i) {
      if ((medavg[i] <= c.threshold) && (medavg[i+1] > c.threshold)) beginTel.push_back(i+1);
      if ((medavg[i] > c.threshold) && (medavg[i+1] <= c.threshold)) endTel.push_back(i+1);
    }
    if (medavg[seqsize - 1] > c.threshold) endTel.push_back(seqsize);
    if (beginTel.size()) {
      //for(uint32_t i = 0; i < seqsize; ++i) {
      //std::cerr << i << ',' << tel[i] << ',' << avg[i] << ',' << medavg[i] << ',' << std::endl;
      //}

      for(uint32_t i = 0; i < beginTel.size(); ++i) {
	//std::cerr << beginTel[i] << ',' << endTel[i] << std::endl;
	uint32_t telLen = endTel[i] - beginTel[i];
	uint32_t telRegions = beginTel.size();
	double avgTelSig = 0;
	for(uint32_t k = beginTel[i]; k < endTel[i]; ++k) avgTelSig+=medavg[k];
	avgTelSig /= (double) telLen;

	//std::cerr << "Telomere\t" << qname << '\t' << (int) fwd << '\t' << beginTel[i] << '\t' << endTel[i] << '\t' << telLen << '\t' << telRegions << '\t' << avgTelSig << '\t' << seqsize << std::endl;
	telreads.push_back(TelomereRecord(fwd, seqsize, beginTel[i], endTel[i], telRegions, avgTelSig, seed, qname));
	foundReps = true;
      }
    }
    return foundReps;
  }


  template<typename TConfig, typename TBitSet>
  inline void
  markRepeats(TConfig const& c, std::string const& sequence, TBitSet& telfwd, TBitSet& telrev) {
    uint32_t seqsize = sequence.size();
    telfwd.resize(seqsize, 0);
    telrev.resize(seqsize, 0);
    telfwd.reset();
    telrev.reset();
	  
    // Mark telomere hits
    std::string m = sequence.substr(0, c.replen);
    uint32_t ct = 0;
    for(uint32_t i = 0; i < m.size(); ++i) ct += mp(m[i]) * std::pow(4, c.replen - i - 1);
    if (c.fwdm[ct]) telfwd[0] = 1;
    if (c.revm[ct]) telrev[0] = 1;
    for(uint32_t i = c.replen; i < seqsize; ++i) {
      ct -= mp(sequence[i-c.replen]) * std::pow(4, c.replen - 1);
      ct *= 4;
      ct += mp(sequence[i]);
      if (c.fwdm[ct]) telfwd[i - c.replen + 1] = 1;
      if (c.revm[ct]) telrev[i - c.replen + 1] = 1;
    }
  }
      
  
  template<typename TConfig>
  inline void
  telomericReads(TConfig const& c, std::vector<Junction>& junctions, std::vector<TelomereRecord>& telreads) {
    typedef std::vector<Junction> TJunction;
    
    // Open file handles
    samFile* samfile = sam_open(c.sample.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.sample.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Parse BAM
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing... " << hdr->target_name[refIndex] << std::endl;

      // Bitsets
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet telfwd;
      TBitSet telrev;

      // Iterate alignments 
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	std::size_t seed = hash_lr(rec);

	// Load sequence
	std::string sequence(rec->core.l_qseq, 'N');
	uint8_t* seqptr = bam_get_seq(rec);
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	
	// Search telomeric repeats
	if (rec->core.l_qseq > (int32_t) c.medwin) {
	  
	  // Mark telomere hits
	  markRepeats(c, sequence, telfwd, telrev);

	  // Identify telomeric repeats
	  if ((telfwd.count() > c.win) || (telrev.count() > c.win)) {
	    std::string qn = bam_get_qname(rec);
	    bool fwdFound = repeatFinder(c, qn, seed, telfwd, 1, telreads);
	    bool revFound = repeatFinder(c, qn, seed, telrev, 0, telreads);

	    // Check all junctions
	    if ((fwdFound) || (revFound)) {
	      std::vector<TelomereRecord> dummy;
	      typename TJunction::iterator iter = std::lower_bound(junctions.begin(), junctions.end(), Junction(seed));
	      for(;(iter < junctions.end()) && (iter->seed == seed); ++iter) {
		if (iter->seqpos >= (int) c.medwin) {
		  std::string subseq = sequence.substr((uint32_t) (iter->seqpos - c.medwin), c.medwin);
		  markRepeats(c, subseq, telfwd, telrev);
		  if ((repeatFinder(c, qn, seed, telfwd, 1, dummy)) || (repeatFinder(c, qn, seed, telrev, 0, dummy))) iter->telLeft = true;
		}
		if ( (int) (iter->seqpos + c.medwin) <= (int) rec->core.l_qseq) {
		  std::string subseq = sequence.substr((uint32_t) iter->seqpos, c.medwin);
		  markRepeats(c, subseq, telfwd, telrev);
		  if ((repeatFinder(c, qn, seed, telfwd, 1, dummy)) || (repeatFinder(c, qn, seed, telrev, 0, dummy))) iter->telRight = true;
		}
	      }
	    }
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


  template<typename TConfig>
  inline void
  outputTelomereSVs(TConfig const& c, std::vector<Junction>& junctions, std::vector<TelomereRecord>& telreads) {
    // Open file handles
    samFile* samfile = sam_open(c.sample.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Map seeds to names
    std::map<std::size_t, std::string> msr;
    
    // Output telomere reads
    std::string filename = c.outprefix + ".reads.tsv";
    std::ofstream rfile(filename.c_str());
    rfile << "read\treadTelStart\treadTelEnd\treadTelSize\treadTelRegion\treadTelSignal\treadSize\ttelMotifStrand\ttelClass" << std::endl;
    for(uint32_t i = 0; i < telreads.size(); ++i) {
      msr.insert(std::make_pair(telreads[i].seed, telreads[i].qname));
      std::string orientation("reverse");
      if (telreads[i].fwd) orientation = "forward";
      uint32_t size = telreads[i].telEnd - telreads[i].telStart;
      std::string label("intra_telomeric");
      if (telreads[i].telStart < c.win) label = "left_telomeric";
      else if (telreads[i].telEnd + c.win > telreads[i].seqsize) label = "right_telomeric";
      rfile << telreads[i].qname << '\t' << (telreads[i].telStart + 1) << '\t' << (telreads[i].telEnd + 1) << '\t' << size << '\t' << telreads[i].telRegions << '\t' << telreads[i].avgTelSig << '\t' << telreads[i].seqsize << '\t' << orientation << '\t' << label << std::endl;
    }
    rfile.close();

    // Output telomeric breakpoints
    filename = c.outprefix + ".breakpoints.tsv";
    std::ofstream bfile(filename.c_str());
    bfile << "chr\tpos\tmapq\tstrand\tclipped\tclass\tread\tseqpos" << std::endl;
    for(uint32_t i = 0; i < junctions.size(); ++i) {
      if ((junctions[i].telLeft) || (junctions[i].telRight)) {
	std::string strand = "reverse";
	if (junctions[i].forward) strand = "forward";
	std::string clipped = "clipped_right";
	if (junctions[i].scleft) clipped = "clipped_left";
	std::string label = "intra_telomeric";
	if (junctions[i].telLeft) label = "left_telomeric";
	else if (junctions[i].telRight) label = "right_telomeric";
	if (junctions[i].scleft) clipped = "clipped_left";
	bfile << hdr->target_name[junctions[i].refidx] << '\t' << junctions[i].refpos << '\t' << (int) junctions[i].qual << '\t' << strand << '\t' << clipped << '\t' << label << '\t' << msr[junctions[i].seed] << '\t' << junctions[i].seqpos << std::endl;
      }
    }
    bfile.close();

    // Clean-up
    bam_hdr_destroy(hdr);
    sam_close(samfile);
  }

  template<typename TConfig>
  inline int32_t
  runTelomere(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Telomere reads
    typedef std::vector<TelomereRecord> TReadTelomere;
    TReadTelomere telreads;

    // Junctions
    typedef std::vector<Junction> TJunctionVector;
    TJunctionVector junctions;

    if (junctions.empty()) {
      // Parse SV breakpoints
      TJunctionVector readBp;
      findJunctions(c, readBp);

      // Telomere reads
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Telomere content screening" << std::endl;
      telomericReads(c, readBp, telreads);

      // Final junctions
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Subset SV reads to telomeric reads" << std::endl;
      std::set<std::size_t> telSeeds;
      for(uint32_t i = 0; i < telreads.size(); ++i) telSeeds.insert(telreads[i].seed);
      for(uint32_t i = 0; i < readBp.size(); ++i) {
	if (telSeeds.find(readBp[i].seed) != telSeeds.end()) junctions.push_back(readBp[i]);
      }
    }

    // Data out
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Output telomere reads" << std::endl;
    outputTelomereSVs(c, junctions, telreads);

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
      ("minclip,c", boost::program_options::value<uint32_t>(&c.minClip)->default_value(18), "min. clipping length")
      ("movavg,m", boost::program_options::value<uint32_t>(&c.win)->default_value(51), "rolling average window")
      ("medsize,s", boost::program_options::value<uint32_t>(&c.medwin)->default_value(501), "rolling median window")
      ("thres,t", boost::program_options::value<double>(&c.threshold)->default_value(0.35), "repeat threshold")
      ("chrlen,l", boost::program_options::value<uint32_t>(&c.minChrLen)->default_value(40000000), "min. chromosome length")
      ("repeats,r", boost::program_options::value<std::string>(&c.repeatStr)->default_value("TTAGGG,TCAGGG,TGAGGG,TTGGGG"), "repeat units")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("out"), "output file prefix")
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
    createRepeatMotifs(c);
    
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
