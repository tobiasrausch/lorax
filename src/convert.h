#ifndef CONVERT_H
#define CONVERT_H

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

  struct ConvertConfig {
    bool seqCoords;
    std::string outprefix;
    boost::filesystem::path gfafile;
    boost::filesystem::path seqfile;
    boost::filesystem::path sample;
    boost::filesystem::path genome;
    boost::filesystem::path readsfile;
  };


  template<typename TConfig>
  inline bool
  plotGraphAlignments(TConfig const& c, Graph const& g, std::vector<AlignRecord> const& aln) {
    typedef std::vector<AlignRecord> TAlignRecords;

    // Vertex map
    std::vector<std::string> idSegment(g.smap.size());
    for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;

    // Load bam file
    samFile* samfile = sam_open(c.readsfile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Load graph sequences
    faidx_t* fai = fai_load(c.seqfile.string().c_str());
    bam1_t* rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
      std::string qname(bam_get_qname(rec));
      std::size_t seed = hash_lr(qname);

      // Find alignment records
      typename TAlignRecords::const_iterator iter = std::lower_bound(aln.begin(), aln.end(), AlignRecord(0, seed), SortAlignRecord<AlignRecord>());
      if ((iter == aln.end()) || (iter->seed != seed)) continue;
	    
      // Load sequence
      std::string sequence;
      sequence.resize(rec->core.l_qseq, 'N');
      uint8_t* seqptr = bam_get_seq(rec);
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
      if (rec->core.flag & BAM_FREVERSE) reverseComplement(sequence);

      // Plot alignment
      for(; ((iter != aln.end()) && (iter->seed == seed)); ++iter) {
	std::string qslice = sequence.substr(iter->qstart, (iter->qend - iter->qstart));
	if (iter->strand == '-') reverseComplement(qslice);
	std::string refslice;
	for(uint32_t i = 0; i < iter->path.size(); ++i) {
	  // Vertex coordinates
	  std::string seqname = idSegment[iter->path[i].tid];
	  int32_t seqlen = 0;
	  char* ref = faidx_fetch_seq(fai, seqname.c_str(), 0, faidx_seq_len(fai, seqname.c_str()), &seqlen);
	  if (!iter->path[i].forward) revcomplement(ref);
	  // Alternate between upper and lower to see vertex breaks
	  if (i%2 == 0) refslice += boost::to_upper_copy(std::string(ref));
	  else refslice += boost::to_lower_copy(std::string(ref));
	  if (i == 0) refslice = refslice.substr(iter->pstart);
	  if (i + 1 == iter->path.size()) refslice = refslice.substr(0, iter->pend - iter->pstart);
	  free(ref);
	}
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
    }
    fai_destroy(fai);

    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(samfile);

    return true;
  }


  template<typename TConfig>
  inline bool
  convertToBam(TConfig const& c, Graph const& g, std::vector<AlignRecord> const& aln) {
    typedef std::vector<AlignRecord> TAlignRecords;

    // Vertex map
    std::vector<std::string> idSegment(g.smap.size());
    for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;
    
    // Load bam file
    samFile* samfile = sam_open(c.readsfile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Output file
    std::ofstream sfile;
    std::string filen = c.outprefix + ".sam";
    sfile.open(filen.c_str());

    // Load graph sequences
    faidx_t* fai = fai_load(c.seqfile.string().c_str());

    // Write header
    for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
      std::string qname(faidx_iseq(fai, refIndex));
      int32_t seqlen = faidx_seq_len(fai, qname.c_str());
      sfile << "@SQ\tSN:" << qname << "\tLN:" << seqlen << std::endl;
    }
    
    // Convert to BAM
    bam1_t* rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
      std::string qname(bam_get_qname(rec));
      std::size_t seed = hash_lr(qname);

      // Find alignment records
      typename TAlignRecords::const_iterator iter = std::lower_bound(aln.begin(), aln.end(), AlignRecord(0, seed), SortAlignRecord<AlignRecord>());
      if ((iter == aln.end()) || (iter->seed != seed)) continue;

      // Load sequence
      std::string sequence(rec->core.l_qseq, 'N');
      std::string quality(rec->core.l_qseq, 'B');
      uint8_t* seqptr = bam_get_seq(rec);
      uint8_t* qualptr = bam_get_qual(rec);
      bool hasQual = true;
      if ((int32_t) qualptr[0] == 255) hasQual = false;
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) {
	if (hasQual) quality[i] = boost::lexical_cast<char>((uint8_t) (qualptr[i] + 33));
	sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
      }
      if (rec->core.flag & BAM_FREVERSE) reverseComplement(sequence);
      
      // Iterate all graph alignments of this read
      for(; ((iter != aln.end()) && (iter->seed == seed)); ++iter) {
	std::string qslice = sequence.substr(iter->qstart, (iter->qend - iter->qstart));
	std::string qualsl = quality.substr(iter->qstart, (iter->qend - iter->qstart));
	
	if (iter->strand == '-') reverseComplement(qslice);
	uint32_t refstart = 0;
	for(uint32_t i = 0; i < iter->path.size(); ++i) {
	  // Vertex coordinates
	  std::string seqname = idSegment[iter->path[i].tid];
	  int32_t seqlen = faidx_seq_len(fai, seqname.c_str());
	  sfile << qname;
	  if (!iter->path[i].forward) sfile << "\t272";
	  else sfile << "\t256";
	  sfile << "\t" << seqname;
	  uint32_t pstart = 0;
	  uint32_t plen = seqlen;
	  if (i == 0) {
	    plen -= iter->pstart;
	    if (iter->path[i].forward) pstart = iter->pstart;
	  }
	  if (i + 1 == iter->path.size()) {
	    plen = iter->pend - iter->pstart - refstart;
	    if (!iter->path[i].forward) pstart = iter->pstart + refstart + seqlen - iter->pend;
	  }
	  sfile << "\t" << pstart + 1;
	  sfile << "\t" << iter->mapq;

	  // Build CIGAR
	  uint32_t refend = refstart + plen;
	  //std::cerr << refstart << ',' << refend << ':' << iter->pstart << ',' << iter->pend << ';' << seqlen << ',' << refstart << std::endl;
	  uint32_t rp = 0;
	  uint32_t sp = 0;
	  std::string cigout = "";
	  std::string qalign = "";
	  std::string qstr = "";
	  if (!hasQual) qstr = "*";
	  for (uint32_t i = 0; i < iter->cigarop.size(); ++i) {
	    if (iter->cigarop[i] == BAM_CEQUAL) {
	      for(uint32_t k = 0; k < iter->cigarlen[i]; ++k, ++sp, ++rp) {
		if ((rp >= refstart) && (rp < refend)) {
		  cigout += "=";
		  qalign += qslice[sp];
		  if (hasQual) qstr += qualsl[sp];
		}
	      }
	    }
	    else if (iter->cigarop[i] == BAM_CDIFF) {
	      for(uint32_t k = 0; k < iter->cigarlen[i]; ++k, ++sp, ++rp) {
		if ((rp >= refstart) && (rp < refend)) {
		  cigout += "X";
		  qalign += qslice[sp];
		  if (hasQual) qstr += qualsl[sp];
		}
	      }
	    }
	    else if (iter->cigarop[i] == BAM_CDEL) {
	      for(uint32_t k = 0; k < iter->cigarlen[i]; ++k, ++rp) {
		if ((rp >= refstart) && (rp < refend)) cigout += "D";
	      }
	    }
	    else if (iter->cigarop[i] == BAM_CINS) {
	      for(uint32_t k = 0; k < iter->cigarlen[i]; ++k, ++sp) {
		if ((rp >= refstart) && (rp < refend)) {
		  cigout += "I";
		  qalign += qslice[sp];
		  if (hasQual) qstr += qualsl[sp];
		}
	      }
	    }
	    else {
	      std::cerr << "Warning: Unknown Cigar option " << iter->cigarop[i] << std::endl;
	      return false;
	    }
	  }
	  // Reverse path?
	  if (!iter->path[i].forward) {
	    reverseComplement(qalign);
	    std::reverse(qstr.begin(), qstr.end());
	    std::reverse(cigout.begin(), cigout.end());
	  }
	  // Cigar run length encoding
	  sfile << '\t';
	  char old = cigout[0];
	  uint32_t clen = 1;
	  for(uint32_t k = 1; k < cigout.size(); ++k) {
	    if (old == cigout[k]) ++clen;
	    else {
	      sfile << clen << old;
	      old = cigout[k];
	      clen = 1;
	    }
	  }
	  sfile << clen << old << "\t*\t0\t0\t" << qalign << "\t" << qstr << std::endl;
	  
	  // Next segment
	  refstart = refend;
	}
      }
    }
    fai_destroy(fai);
    sfile.close();

    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(samfile);

    return true;
  }

  template<typename TConfig>
  inline int32_t
  runConvert(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Load pan-genome graph
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Load pan-genome graph" << std::endl;
    Graph g;
    parseGfa(c, g);
    
    // Parse alignments
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse alignments" << std::endl;
    std::vector<AlignRecord> aln;
    parseGaf(c, g, aln);

    // Plot pair-wise graph alignments
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parse reads" << std::endl;
    //if (!plotGraphAlignments(c, g, aln)) return -1;
    if (!convertToBam(c, g, aln)) return -1;
    
    // Write pan-genome graph
    //std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Write GFA" << std::endl;
    //writeGfa(c, g);

#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  
  int convert(int argc, char** argv) {
    ConvertConfig c;
    c.seqCoords = false;
    
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("graph,g", boost::program_options::value<boost::filesystem::path>(&c.gfafile), "GFA pan-genome graph")
      ("sequences,s", boost::program_options::value<boost::filesystem::path>(&c.seqfile), "stable sequences")
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "FASTA reference")
      ("align,a", boost::program_options::value<boost::filesystem::path>(&c.readsfile), "BAM file")
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
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("graph"))) {
      std::cout << "Usage:" << std::endl;
      std::cout << "lorax " << argv[0] << " [OPTIONS] -g <pangenome.hg38.gfa.gz> -r <genome.fasta> -a <align.bam> <sample.gaf.gz>" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Sequence file for vertex segments
    c.seqfile = c.outprefix + ".fa";
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "lorax ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runConvert(c);
  }

}

#endif
