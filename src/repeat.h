#ifndef REPEAT_H
#define REPEAT_H

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

  struct RepeatConfig {
    typedef std::set<std::string> TStringSet;
    
    uint32_t period;
    uint32_t replen;
    uint32_t minChrLen;
    uint32_t window;
    uint32_t chrEnd;
    std::string repeatStr;
    TStringSet fwdmotif;
    TStringSet revmotif;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
  };

  template<typename TConfig>
  inline int32_t
  runRepeatFinder(TConfig const& c) {

#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Output file
    std::ofstream ofile(c.outfile.string().c_str());
    ofile << "chr\tstart\tend\tN\tGC\tfwd\trev\trepeat" << std::endl;
    
    // Parse chr-by-chr
    char* seq = NULL;
    faidx_t* fai = fai_load(c.genome.string().c_str());
    uint32_t motiflen = c.period * c.replen;
    uint64_t allbases = 0;
    uint64_t repbases = 0;
    for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
      // Check sequence length
      std::string tname(faidx_iseq(fai, refIndex));
      uint32_t seqlen = faidx_seq_len(fai, tname.c_str());
      if (seqlen < c.minChrLen) continue;

      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Processing... " << tname << std::endl;

      // Counters
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet fwdhit(seqlen);
      TBitSet revhit(seqlen);
      TBitSet nref(seqlen);
      TBitSet gcref(seqlen);

      // Load sequence
      int32_t sl = -1;
      seq = faidx_fetch_seq(fai, tname.c_str(), 0, seqlen, &sl);

      // Iterate chromosome
      for(uint32_t pos = 0; pos < seqlen; ++pos) {
	if ((seq[pos] == 'c') || (seq[pos] == 'C') || (seq[pos] == 'g') || (seq[pos] == 'G')) gcref[pos] = 1;
	if ((seq[pos] == 'n') || (seq[pos] == 'N')) nref[pos] = 1;
	if (pos + motiflen < seqlen) {
	  if (c.fwdmotif.find(boost::to_upper_copy(std::string(seq + pos, seq + pos + motiflen))) != c.fwdmotif.end()) fwdhit[pos] = 1;
	  if (c.revmotif.find(boost::to_upper_copy(std::string(seq + pos, seq + pos + motiflen))) != c.revmotif.end()) revhit[pos] = 1;
	}
      }
      if (seqlen >= motiflen + nref.count()) {
	allbases += (seqlen - motiflen - nref.count());
	repbases += fwdhit.count() + revhit.count();
      }

      // Summarize
      for(uint32_t pos = 0; pos + c.window < seqlen; pos += c.window) {
	uint32_t nsum = 0;
	uint32_t gcsum = 0;
	uint32_t fwdsum = 0;
	uint32_t revsum = 0;
	for(uint32_t k = pos; k < pos + c.window; ++k) {
	  nsum += nref[k];
	  gcsum += gcref[k];
	  fwdsum += fwdhit[k];
	  revsum += revhit[k];
	}
	if ((!c.chrEnd) || ((pos < c.chrEnd) || (pos + c.chrEnd > seqlen))) {
	  ofile << tname << "\t" << pos << "\t" << (pos + c.window) << "\t" << nsum << "\t" << (double) gcsum / (double) c.window << "\t" << fwdsum << "\t" << revsum << "\t" << c.repeatStr << std::endl;
	}
      }  
	  

      if (seq != NULL) free(seq);
    }
    // Close output file
    ofile.close();
      
#ifdef PROFILE
    ProfilerStop();
#endif

    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Repeat nucleotides: " << repbases << ", All nucleotides: " << allbases << std::endl;
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    
    return 0;
  }
  

  int repeat(int argc, char** argv) {
    RepeatConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Options");
    generic.add_options()
      ("help,?", "show help message")
      ("chrlen,l", boost::program_options::value<uint32_t>(&c.minChrLen)->default_value(40000000), "min. chromosome length")
      ("repeats,r", boost::program_options::value<std::string>(&c.repeatStr)->default_value("TTAGGG,TCAGGG,TGAGGG,TTGGGG"), "repeat units")
      ("period,p", boost::program_options::value<uint32_t>(&c.period)->default_value(3), "repeat period")
      ("window,w", boost::program_options::value<uint32_t>(&c.window)->default_value(1000), "window length")
      ("chrend,e", boost::program_options::value<uint32_t>(&c.chrEnd)->default_value(0), "chromosome end window [0: deactivate]")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.tsv"), "output file")
      ("nomix,n", "do not mix repeat units")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.genome), "input file")
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
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << "Usage: lorax " << argv[0] << " [OPTIONS] <ref.fa>" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Generate motifs
    bool mix = true;
    if (vm.count("nomix")) mix = false;
    createRepeatMotifs(c, mix);
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "lorax ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runRepeatFinder(c);
  }

}

#endif
