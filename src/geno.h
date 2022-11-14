#ifndef GENO_H
#define GENO_H

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

  struct GenoConfig {
    float pct;
    uint32_t len;
    uint32_t delcut;
    std::string outprefix;
    boost::filesystem::path genome;
    boost::filesystem::path gfafile;
    boost::filesystem::path seqfile;
    boost::filesystem::path sample;
    boost::filesystem::path readsfile;
  };


  template<typename TConfig>
  inline bool
  plotGraphAlignments(TConfig const& c, std::vector<AlignRecord> const& aln) {
    typedef std::vector<AlignRecord> TAlignRecords;
    
    faidx_t* fai = fai_load(c.readsfile.string().c_str());
    faidx_t* sai = fai_load(c.seqfile.string().c_str());
    for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
      std::string qname(faidx_iseq(fai, refIndex));
      int32_t seqlen = 0;
      char* query = faidx_fetch_seq(fai, qname.c_str(), 0, faidx_seq_len(fai, qname.c_str()), &seqlen);
      std::size_t seed = hash_lr(qname);

      // Find alignment records
      typename TAlignRecords::const_iterator iter = std::lower_bound(aln.begin(), aln.end(), AlignRecord(0, seed), SortAlignRecord<AlignRecord>());
      for(; ((iter != aln.end()) && (iter->seed == seed)); ++iter) {
	// Query slice
	std::string qslice = std::string(query + iter->qstart, query + iter->qend);
	if (iter->strand == '-') reverseComplement(qslice);
	std::string refslice;
	for(uint32_t i = 0; i < iter->path.size(); ++i) {
	  int32_t segment = iter->path[i];
	  bool forward = true;
	  if (segment < 0) {
	    forward = false;
	    segment *= -1;
	  }
	  std::string segname =  boost::lexical_cast<std::string>(segment);
	  seqlen = 0;
	  char* ref = faidx_fetch_seq(sai, segname.c_str(), 0, faidx_seq_len(sai, segname.c_str()), &seqlen);
	  if (forward) refslice += std::string(ref);
	  else {
	    revcomplement(ref);
	    refslice += std::string(ref);
	  }
	  free(ref);
	}
	if (refslice.size() != (uint32_t) iter->plen) {
	  std::cerr << "Path length does not match!" << std::endl;
	  return false;
	}
	refslice = refslice.substr(iter->pstart, iter->pend - iter->pstart);
	// Show base-level alignments
	uint32_t rp = 0;
	uint32_t sp = 0;
	uint32_t al = 0;
	std::string refalign;
	std::string spacer;
	std::string qalign;
	for (uint32_t i = 0; i < iter->cigarop.size(); ++i) {
	  if (iter->cigarop[i] == BAM_CEQUAL) {
	    for(uint32_t k = 0; k < iter->cigarlen[i]; ++k, ++al, ++sp, ++rp) {
	      refalign += refslice[rp];
	      spacer += "|";
	      qalign += qslice[sp];
	    }
	  }
	  else if (iter->cigarop[i] == BAM_CDIFF) {
	    for(uint32_t k = 0; k < iter->cigarlen[i]; ++k, ++al, ++sp, ++rp) {
	      refalign += refslice[rp];
	      spacer += " ";
	      qalign += qslice[sp];
	    }
	  }
	  else if (iter->cigarop[i] == BAM_CDEL) {
	    for(uint32_t k = 0; k < iter->cigarlen[i]; ++k, ++al, ++rp) {
	      refalign += refslice[rp];
	      spacer += " ";
	      qalign += "-";
	    }
	  }
	  else if (iter->cigarop[i] == BAM_CINS) {
	    for(uint32_t k = 0; k < iter->cigarlen[i]; ++k, ++al, ++sp) {
	      refalign += "-";
	      spacer += " ";
	      qalign += qslice[sp];
	    }
	  }
	  else if (iter->cigarop[i] == BAM_CSOFT_CLIP) {
	    std::cerr << "Soft-clips!" << std::endl;
	    return false;
	  }	    
	  else if (iter->cigarop[i] == BAM_CREF_SKIP) {
	    std::cerr << "Reference skip!" << std::endl;
	    return false;
	  }
	  else if (iter->cigarop[i] == BAM_CHARD_CLIP) {
	    std::cerr << "Hard-clips!" << std::endl;
	    return false;
	  }
	  else {
	    std::cerr << "Warning: Unknown Cigar option " << iter->cigarop[i] << std::endl;
	    return false;
	  }
	}
	std::cerr << ">" << qname << std::endl;
	std::cerr << refalign << std::endl;
	std::cerr << spacer << std::endl;
	std::cerr << qalign << std::endl;
	if (refalign.size() != (uint32_t) iter->alignlen) {
	  std::cerr << "Alignment length does not match!" << std::endl;
	  return false;
	}
      }
      free(query);
    }
    fai_destroy(sai);
    fai_destroy(fai);
    return true;
  }


  template<typename TConfig>
  inline int32_t
  runGeno(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Load pan-genome graph
    //std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Load pan-genome graph" << std::endl;
    //Graph g;
    //parseGfa(c, g);
    
    // Parse alignments
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse alignments" << std::endl;
    //pctGaf(c);
    
    // Parse reads
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse reads" << std::endl;

    // Plot pair-wise graph alignments
    //if (!plotGraphAlignments(c, aln)) return -1;
    
    // Write pan-genome graph
    //writeGfa(c, g);

#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  
  int geno(int argc, char** argv) {
    GenoConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("graph,g", boost::program_options::value<boost::filesystem::path>(&c.gfafile), "GFA pan-genome graph")
      ("reads,e", boost::program_options::value<boost::filesystem::path>(&c.readsfile), "reads in FASTA format")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("out"), "output tprefix")
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
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << "Usage:" << std::endl;
      std::cout << "lorax " << argv[0] << " [OPTIONS] -g <pangenome.hg38.gfa.gz> -e <reads.fasta> <sample.gaf>" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }


    // Random name for temporary file
    boost::uuids::uuid uuid = boost::uuids::random_generator()();
    c.seqfile = c.outprefix + "." + boost::lexical_cast<std::string>(uuid) + ".fa";

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "lorax ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runGeno(c);
  }

}

#endif
