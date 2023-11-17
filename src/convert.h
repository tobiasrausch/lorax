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
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include "gfa.h"
#include "gaf.h"

namespace lorax
{

  struct ConvertConfig {
    bool hasFastq;
    uint32_t chunk;
    boost::filesystem::path outfile;
    boost::filesystem::path gfafile;
    boost::filesystem::path seqfile;
    boost::filesystem::path sample;
    boost::filesystem::path genome;
    boost::filesystem::path fastqfile;
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
	  std::string seqname = idSegment[iter->path[i].second];
	  int32_t seqlen = 0;
	  char* ref = faidx_fetch_seq(fai, seqname.c_str(), 0, faidx_seq_len(fai, seqname.c_str()), &seqlen);
	  //std::cerr << "CP\t" << seqlen << '\t' << g.segments[iter->path[i].second].len << std::endl;
	  if (!iter->path[i].first) revcomplement(ref);
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


  template<typename TAlignIter>
  inline bool
  convertToBam(std::vector<std::string> const& idSegment, std::size_t const seed, faidx_t* fai, std::ostream& sfile, std::string const& qname, std::string const& sequence, std::string const& quality, std::vector<AlignRecord> const& aln, TAlignIter& iter) {
    bool hasQual = true;
    if ((quality.size() == 1) && (quality == "*")) hasQual = false;

    // Identify longest alignment to flag as primary
    int32_t qbestsize = 0;
    for(TAlignIter tmpIter = iter; ((tmpIter != aln.end()) && (tmpIter->seed == seed)); ++tmpIter) {
      if ((tmpIter->qend - tmpIter->qstart) > qbestsize) qbestsize = tmpIter->qend - tmpIter->qstart;
    }
    
    // Iterate all graph alignments of this read
    bool primdone = false;
    for(; ((iter != aln.end()) && (iter->seed == seed)); ++iter) {
      std::string qslice = sequence.substr(iter->qstart, (iter->qend - iter->qstart));
      std::string qualsl;
      if (hasQual) qualsl = quality.substr(iter->qstart, (iter->qend - iter->qstart));
	
      if (iter->strand == '-') reverseComplement(qslice);
      uint32_t refstart = 0;
      for(uint32_t i = 0; i < iter->path.size(); ++i) {
	bool primrec = false; // Primary alignment record
	
	// Vertex coordinates
	std::string seqname = idSegment[iter->path[i].second];
	int32_t seqlen = faidx_seq_len(fai, seqname.c_str());
	std::string outstr = qname;
	if (!iter->path[i].first) {
	  if ((!primdone) && ((iter->qend - iter->qstart) == qbestsize)) {
	    primrec = true;
	    outstr += "\t16";
	  } else outstr += "\t272";
	} else {
	  if ((!primdone) && ((iter->qend - iter->qstart) == qbestsize)) {
	    primrec = true;
	    outstr += "\t0";
	  } else outstr += "\t256";
	}
	outstr += "\t";
	outstr += seqname;
	uint32_t pstart = 0;
	uint32_t plen = seqlen;
	if (i == 0) {
	  plen -= iter->pstart;
	  if (iter->path[i].first) pstart = iter->pstart;
	}
	if (i + 1 == iter->path.size()) {
	  plen = iter->pend - iter->pstart - refstart;
	  if (!iter->path[i].first) {
	    if (i == 0) pstart = seqlen - iter->pend;
	    else pstart = iter->pstart + refstart + seqlen - iter->pend;
	  }
	}
	
	// Build CIGAR
	uint32_t refend = refstart + plen;
	//std::cerr << refstart << ',' << refend << ':' << iter->pstart << ',' << iter->pend << ';' << seqlen << ',' << refstart << std::endl;
	uint32_t rp = 0;
	uint32_t sp = 0;
	std::string cigout = "";
	for(int32_t clip = 0; clip < iter->qstart; ++clip) {
	  if (primrec) cigout += "S";
	  else cigout += "H";
	}
	std::string qalign = "";
	std::string qstr = "";
	if (!hasQual) qstr = "*";
	for (uint32_t ci = 0; ci < iter->cigarop.size(); ++ci) {
	  if (iter->cigarop[ci] == BAM_CEQUAL) {
	    for(uint32_t k = 0; k < iter->cigarlen[ci]; ++k, ++sp, ++rp) {
	      if ((rp >= refstart) && (rp < refend)) {
		cigout += "=";
		qalign += qslice[sp];
		if (hasQual) qstr += qualsl[sp];
	      } else {
		if (primrec) cigout += "S";
		else cigout += "H";
	      }
	    }
	  }
	  else if (iter->cigarop[ci] == BAM_CDIFF) {
	    for(uint32_t k = 0; k < iter->cigarlen[ci]; ++k, ++sp, ++rp) {
	      if ((rp >= refstart) && (rp < refend)) {
		cigout += "X";
		qalign += qslice[sp];
		if (hasQual) qstr += qualsl[sp];
	      } else {
		if (primrec) cigout += "S";
		else cigout += "H";
	      }
	    }
	  }
	  else if (iter->cigarop[ci] == BAM_CDEL) {
	    for(uint32_t k = 0; k < iter->cigarlen[ci]; ++k, ++rp) {
	      if ((rp >= refstart) && (rp < refend)) cigout += "D";
	    }
	  }
	  else if (iter->cigarop[ci] == BAM_CINS) {
	    for(uint32_t k = 0; k < iter->cigarlen[ci]; ++k, ++sp) {
	      if ((rp >= refstart) && (rp < refend)) {
		cigout += "I";
		qalign += qslice[sp];
		if (hasQual) qstr += qualsl[sp];
	      } else {
		if (primrec) cigout += "S";
		else cigout += "H";
	      }
	    }
	  }
	  else {
	    std::cerr << "Warning: Unknown Cigar option " << iter->cigarop[ci] << std::endl;
	    return false;
	  }
	}
	// End clipping
	for(int32_t clip = iter->qend; clip < (int32_t) sequence.size(); ++clip) {
	  if (primrec) cigout += "S";
	  else cigout += "H";
	}
	if (primrec) {
	  qalign = sequence; // Assign full sequence to primary record
	  if (hasQual) qstr = quality;
	}
	// Reverse path?
	if (!iter->path[i].first) {
	  reverseComplement(qalign);
	  std::reverse(qstr.begin(), qstr.end());
	  std::reverse(cigout.begin(), cigout.end());
	}

	// Any alignments to this node?
	if (qalign.size()) {
	  sfile << outstr;
	  sfile << "\t" << pstart + 1;
	  sfile << "\t" << iter->mapq;
	
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
	}

	// Primary alignment set?
	if (primrec) primdone = true;
	
	// Next segment
	refstart = refend;
      }
    }
    return true;
  }
    
  
  template<typename TConfig>
  inline bool
  convertToBamViaFASTQ(TConfig const& c, Graph const& g, std::vector<AlignRecord> const& aln, bool const firstPass) {
    typedef std::vector<AlignRecord> TAlignRecords;

    // Vertex map
    std::vector<std::string> idSegment(g.smap.size());
    for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;

    // Output file
    std::streambuf* buf;
    std::ofstream of;
    if (c.outfile != "-") {
      of.open(c.outfile.string().c_str());
      buf = of.rdbuf();
    } else {
      buf = std::cout.rdbuf();
    }
    std::ostream sfile(buf);    

    // Load graph sequences
    faidx_t* fai = fai_load(c.seqfile.string().c_str());

    // Write header
    if (firstPass) {
      for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
	std::string qname(faidx_iseq(fai, refIndex));
	int32_t seqlen = faidx_seq_len(fai, qname.c_str());
	sfile << "@SQ\tSN:" << qname << "\tLN:" << seqlen << std::endl;
      }
    }

    // Load FASTA/FASTQ
    std::ifstream fqfile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    if (is_gz(c.fastqfile)) {
      fqfile.open(c.fastqfile.string().c_str(), std::ios_base::in | std::ios_base::binary);
      dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
    } else fqfile.open(c.fastqfile.string().c_str(), std::ios_base::in);
    dataIn.push(fqfile);
    std::istream instream(&dataIn);
    std::string gline;
    uint64_t lnum = 0;
    std::string qname;
    bool validRec = true;
    while(std::getline(instream, gline)) {
      if (lnum % 2 == 0) {
	// FASTA or FASTQ
	if ((gline[0] == '>') || (gline[0] == '@')) {
	  validRec = true;
	  qname = gline.substr(1);
	  qname = qname.substr(0, qname.find(' '));
	  qname = qname.substr(0, qname.find('\t'));
	} else validRec = false;
      } else if (lnum % 2 == 1) {
	if (validRec) {
	  std::size_t seed = hash_lr(qname);
	
	  // Find alignment records
	  typename TAlignRecords::const_iterator iter = std::lower_bound(aln.begin(), aln.end(), AlignRecord(0, seed), SortAlignRecord<AlignRecord>());
	  if ((iter != aln.end()) || (iter->seed == seed)) {
	    // Convert alignments
	    if (!convertToBam(idSegment, seed, fai, sfile, qname, gline, "*", aln, iter)) return false;
	  }
	}
      }
      ++lnum;
    }
    // Clean-up
    dataIn.pop();
    if (is_gz(c.fastqfile)) dataIn.pop();
    fqfile.close();
    
    return true;
  }

  template<typename TConfig>
  inline bool
  convertToBamViaCRAM(TConfig const& c, Graph const& g, std::vector<AlignRecord> const& aln, bool const firstPass) {
    typedef std::vector<AlignRecord> TAlignRecords;

    // Vertex map
    std::vector<std::string> idSegment(g.smap.size());
    for(typename Graph::TSegmentIdMap::const_iterator it = g.smap.begin(); it != g.smap.end(); ++it) idSegment[it->second] = it->first;
    
    // Load bam file
    samFile* samfile = sam_open(c.readsfile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Output file
    std::streambuf* buf;
    std::ofstream of;
    if (c.outfile != "-") {
      of.open(c.outfile.string().c_str());
      buf = of.rdbuf();
    } else {
      buf = std::cout.rdbuf();
    }
    std::ostream sfile(buf);    

    // Load graph sequences
    faidx_t* fai = fai_load(c.seqfile.string().c_str());

    // Write header
    if (firstPass) {
      for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
	std::string qname(faidx_iseq(fai, refIndex));
	int32_t seqlen = faidx_seq_len(fai, qname.c_str());
	sfile << "@SQ\tSN:" << qname << "\tLN:" << seqlen << std::endl;
      }
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
      std::string quality;
      uint8_t* seqptr = bam_get_seq(rec);
      uint8_t* qualptr = bam_get_qual(rec);
      bool hasQual = true;
      if ((int32_t) qualptr[0] == 255) {
	hasQual = false;
	quality = "*";
      } else quality.resize(rec->core.l_qseq, 'B');
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) {
	if (hasQual) quality[i] = boost::lexical_cast<char>((uint8_t) (qualptr[i] + 33));
	sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
      }
      if (rec->core.flag & BAM_FREVERSE) {
	reverseComplement(sequence);
	std::reverse(quality.begin(), quality.end());
      }
      
      // Convert alignments
      if (!convertToBam(idSegment, seed, fai, sfile, qname, sequence, quality, aln, iter)) return false;
    }
    fai_destroy(fai);
    if (c.outfile != "-") of.close();

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
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Load pan-genome graph" << std::endl;
    Graph g;
    parseGfa(c, g, true);

    // Parse alignments
    
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Convert alignments" << std::endl;
    std::vector<AlignRecord> aln;
    bool firstPass = true;
    std::set<std::size_t> processed;
    bool gafdone = true;
    do {
      gafdone = parseGaf(c, g, aln, c.chunk, processed);
      std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parsed " << aln.size() << " alignments" << std::endl;
      //if (!plotGraphAlignments(c, g, aln)) return -1;
      if (c.hasFastq) {
	if (!convertToBamViaFASTQ(c, g, aln, firstPass)) {
	  std::cerr << "Error converting to BAM!" << std::endl;
	  return -1;
	}
      } else {
	if (!convertToBamViaCRAM(c, g, aln, firstPass)) {
	  std::cerr << "Error converting to BAM!" << std::endl;
	  return -1;
	}
      }

      // Clean-up
      if (c.chunk) std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Processed " << processed.size() << " reads" << std::endl;
      firstPass = false;
      aln.clear();
    } while (!gafdone);

#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }


  
  int convert(int argc, char** argv) {
    ConvertConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("chunk,c", boost::program_options::value<uint32_t>(&c.chunk)->default_value(500000), "chunk size [0: all at once]")
      ("graph,g", boost::program_options::value<boost::filesystem::path>(&c.gfafile), "GFA pan-genome graph")
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "FASTA reference")
      ("align,a", boost::program_options::value<boost::filesystem::path>(&c.readsfile), "BAM/CRAM file")
      ("fastq,f", boost::program_options::value<boost::filesystem::path>(&c.fastqfile), "FASTA/FASTQ file")
      ("sequences,s", boost::program_options::value<boost::filesystem::path>(&c.seqfile)->default_value("out.fa"), "output sequences")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output alignments")
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
      std::cerr << "Usage:" << std::endl;
      std::cerr << "Using BAM/CRAM: lorax " << argv[0] << " [OPTIONS] -g <pangenome.gfa.gz> -r <genome.fasta> -a <align.bam> <sample.gaf.gz>" << std::endl;
      std::cerr << "Using FASTQ: lorax " << argv[0] << " [OPTIONS] -g <pangenome.gfa.gz> -f <reads.fa.gz> <sample.gaf.gz>" << std::endl;
      std::cerr << visible_options << "\n";
      return -1;
    }

    // Check outfile
    if (!vm.count("outfile")) c.outfile = "-";
    else {
      if (c.outfile.string() != "-") {
	if (!_outfileValid(c.outfile)) return 1;
      }
    }

    // FASTQ or BAM/CRAM
    if (vm.count("fastq")) c.hasFastq = true;
    else {
      if ((vm.count("reference")) && (vm.count("align"))) c.hasFastq = false;
      else {
	std::cerr << "You either need to provide a FASTQ file or a BAM/CRAM file with the associated linear reference!" << std::endl;
	return -1;
      }
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cerr << "lorax ";
    for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
    std::cerr << std::endl;
    
    return runConvert(c);
  }

}

#endif
