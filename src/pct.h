#ifndef PCT_H
#define PCT_H

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

  struct PctConfig {
    bool gfaMode;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
    boost::filesystem::path sample;
  };

  template<typename TConfig>
  inline void
  pctGaf(TConfig const& c) {
    // Output file
    std::ofstream ofile;
    ofile.open(c.outfile.string().c_str(), std::ofstream::out | std::ofstream::trunc);
    ofile << "qname\tqlen\tqsublen\tpctidentity\tlargestdel\tlargestins\tmapped\tmatches\tmismatches\tdeletions\tdelsize\tinsertions\tinssize\tsoftclips\tsoftclipsize\thardclips\thardclipsize" << std::endl;
    
    // Parse GAF
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
      AlignRecord ar;
		   
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep("\t");
      Tokenizer tokens(gline, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter != tokens.end()) {
	std::string qname = *tokIter;
	ar.seed = hash_lr(qname);
	++tokIter;
	if (tokIter != tokens.end()) {
	  ar.qlen = boost::lexical_cast<int32_t>(*tokIter);
	  ++tokIter;
	  if (tokIter != tokens.end()) {
	    ar.qstart = boost::lexical_cast<int32_t>(*tokIter);
	    ++tokIter;
	    if (tokIter != tokens.end()) {
	      ar.qend = boost::lexical_cast<int32_t>(*tokIter);
	      ++tokIter;
	      if (tokIter != tokens.end()) {
		ar.strand = boost::lexical_cast<char>(*tokIter);
		++tokIter;
		if (tokIter != tokens.end()) {
		  // Skip path parsing
		  ++tokIter;
		  if (tokIter != tokens.end()) {
		    ar.plen = boost::lexical_cast<int32_t>(*tokIter);
		    ++tokIter;
		    if (tokIter != tokens.end()) {
		      ar.pstart = boost::lexical_cast<int32_t>(*tokIter);
		      ++tokIter;
		      if (tokIter != tokens.end()) {
			ar.pend = boost::lexical_cast<int32_t>(*tokIter);
			++tokIter;
			if (tokIter != tokens.end()) {
			  ar.matches = boost::lexical_cast<int32_t>(*tokIter);
			  ++tokIter;
			  if (tokIter != tokens.end()) {
			    ar.alignlen = boost::lexical_cast<int32_t>(*tokIter);
			    ++tokIter;
			    if (tokIter != tokens.end()) {
			      ar.mapq = boost::lexical_cast<int32_t>(*tokIter);
			      ++tokIter;
			      for(; tokIter != tokens.end(); ++tokIter) {
				// Optional fields
				boost::char_separator<char> kvsep(":");
				Tokenizer tokens(*tokIter, kvsep);
				Tokenizer::iterator tikv = tokens.begin();
				if (*tikv == "cg") {
				  ++tikv; ++tikv;
				  parseGafCigar(*tikv, ar);
				}
			      }
			      
			      
			      // Evaluate alignment record
			      uint32_t largestdel = 0;
			      uint32_t largestins = 0;
			      uint32_t mismatch = 0;
			      uint32_t match = 0;
			      uint32_t del = 0;
			      uint32_t delsize = 0;
			      uint32_t ins = 0;
			      uint32_t inssize = 0;
			      uint32_t sc = 0;
			      uint32_t scsize = 0;
			      uint32_t hc = 0;
			      uint32_t hcsize = 0;
			      for (uint32_t i = 0; i < ar.cigarop.size(); ++i) {
				if (ar.cigarop[i] == BAM_CEQUAL) match += ar.cigarlen[i];
				else if (ar.cigarop[i] == BAM_CDIFF) mismatch += ar.cigarlen[i];
				else if (ar.cigarop[i] == BAM_CDEL) {
				  if (ar.cigarlen[i] > largestdel) largestdel = ar.cigarlen[i];
				  ++del;
				  delsize += ar.cigarlen[i];
				}
				else if (ar.cigarop[i] == BAM_CINS) {
				  if (ar.cigarlen[i] > largestins) largestins = ar.cigarlen[i];
				  ++ins;
				  inssize += ar.cigarlen[i];
				}
				else {
				  std::cerr << "Warning: Unknown Cigar option " << ar.cigarop[i] << std::endl;
				  exit(-1);
				}
			      }
			      
			      // Percent identity
			      int32_t qsublen = ar.qend - ar.qstart;
			      double pctval = (double) (match) / (double) qsublen;
			      ofile << qname << '\t' << ar.qlen << '\t' << qsublen << '\t' << pctval << '\t' << largestdel << '\t' << largestins << "\taligned\t" << match << '\t' << mismatch << '\t' << del << '\t' << delsize << '\t' << ins << '\t' << inssize << '\t' << sc << '\t' << scsize << '\t' << hc << '\t' << hcsize << std::endl;
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
    dataIn.pop();
    if (is_gz(c.sample)) dataIn.pop();
    gafFile.close();
  }


  
  template<typename TConfig>
  inline void
  pctBam(TConfig const& c) {
    // Open file handles
    samFile* samfile = sam_open(c.sample.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Output file
    std::ofstream ofile;
    ofile.open(c.outfile.string().c_str(), std::ofstream::out | std::ofstream::trunc);
    ofile << "qname\tqlen\tqsublen\tpctidentity\tlargestdel\tlargestins\tmapped\tmatches\tmismatches\tdeletions\tdelsize\tinsertions\tinssize\tsoftclips\tsoftclipsize\thardclips\thardclipsize" << std::endl;

    // Parse BAM
    int32_t refIndex = -1;
    char* seq = NULL;
    faidx_t* fai = fai_load(c.genome.string().c_str());
    bam1_t* rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      // New chromosome?
      if ((!(rec->core.flag & BAM_FUNMAP)) && (rec->core.tid != refIndex)) {
	// Load new chromosome
	if (seq != NULL) free(seq);
	refIndex = rec->core.tid;
	int32_t seqlen = -1;
	std::string tname(hdr->target_name[refIndex]);
	seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);

	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing... " << hdr->target_name[refIndex] << std::endl;
      }
	
      // Parse only primary and supplementary alignments
      if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY)) continue;

      // Get the read sequence
      std::string sequence;
      sequence.resize(rec->core.l_qseq);
      uint8_t* seqptr = bam_get_seq(rec);
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

      
      // Get read qualities
      typedef std::vector<uint8_t> TQuality;
      TQuality quality;
      quality.resize(rec->core.l_qseq, (uint8_t) 33);
      uint8_t* qualptr = bam_get_qual(rec);
      if (!(qualptr[0] == 0xff)) {
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
      }

      // Unmapped read
      if (rec->core.flag & (BAM_FUNMAP)) {
	ofile << bam_get_qname(rec) << '\t' << sequence.size() << "\t0\t0\t0\t0\tunmapped" << std::endl;
	continue;
      }
	
      // Parse cigar
      uint32_t largestdel = 0;
      uint32_t largestins = 0;
      uint32_t rp = rec->core.pos; // reference pointer
      uint32_t sp = 0; // sequence pointer
      uint32_t mismatch = 0;
      uint32_t match = 0;
      uint32_t del = 0;
      uint32_t delsize = 0;
      uint32_t ins = 0;
      uint32_t inssize = 0;
      uint32_t sc = 0;
      uint32_t scsize = 0;
      uint32_t hc = 0;
      uint32_t hcsize = 0;
      uint32_t* cigar = bam_get_cigar(rec);
      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	  for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
	    if (sequence[sp] == upper(seq[rp])) ++match;
	    else ++mismatch;
	    ++sp;
	    ++rp;
	  }
	} else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	  if (bam_cigar_oplen(cigar[i]) > largestdel) largestdel = bam_cigar_oplen(cigar[i]);
	  ++del;
	  delsize += bam_cigar_oplen(cigar[i]);
	  rp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	  if (bam_cigar_oplen(cigar[i]) > largestins) largestins = bam_cigar_oplen(cigar[i]);
	  ++ins;
	  inssize += bam_cigar_oplen(cigar[i]);
	  sp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	  ++sc;
	  scsize += bam_cigar_oplen(cigar[i]);
	  sp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	  rp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	  ++hc;
	  hcsize += bam_cigar_oplen(cigar[i]);
	} else {
	  std::cerr << "Warning: Unknown Cigar options!" << std::endl;
	}
      }

      // Percent identity
      uint32_t qsublen = querySubLength(rec);
      uint32_t qlen = sequenceLength(rec);
      double pctval = (double) (match) / (double) qsublen;
      ofile << bam_get_qname(rec) << '\t' << qlen << '\t' << qsublen << '\t' << pctval << '\t' << largestdel << '\t' << largestins << "\taligned\t" << match << '\t' << mismatch << '\t' << del << '\t' << delsize << '\t' << ins << '\t' << inssize << '\t' << sc << '\t' << scsize << '\t' << hc << '\t' << hcsize << std::endl;
    }
    if (seq != NULL) free(seq);
    ofile.close();
    
    // Clean-up
    bam_destroy1(rec);
    fai_destroy(fai);
    bam_hdr_destroy(hdr);
    sam_close(samfile);
  }
    

  template<typename TConfig>
  inline int32_t
  runPct(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Parse alignments
    if (c.gfaMode) {
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse GAF alignments" << std::endl;
      pctGaf(c);
    } else {
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse BAM alignments" << std::endl;
      pctBam(c);
    }
    

#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  
  int pct(int argc, char** argv) {
    PctConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.tsv"), "output statistics")
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
      std::cout << "Linear reference genome: lorax " << argv[0] << " [OPTIONS] -r <ref.fa> <sample.bam>" << std::endl;
      std::cout << "Pan-genome graph: lorax " << argv[0] << " [OPTIONS] <sample.gaf.gz>" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Pan-genome or linear reference
    if (!vm.count("reference")) c.gfaMode = true;
    else c.gfaMode = false;

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "lorax ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runPct(c);
  }

}

#endif
